#!/usr/bin/env python3
"""Vendor the immutable primary sources for the 2026 SIDIS tranches.

The command is intentionally separate from ordinary ``fetch-data`` validation:
it is a maintainer tool for creating a new pinned source snapshot.  Runtime
validation never rewrites tracked inputs.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
import shutil
import tempfile
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Mapping


ROOT = Path(__file__).resolve().parents[1]
RETRIEVED = "2026-08-28"
HEPDATA = {
    "COMPASS_2026_I3096394": {
        "inspire": 3096394,
        "record": 169859,
        "doi": "10.17182/hepdata.169859.v1",
        "tables": 3,
    },
    "COMPASS_2025_I2840545": {
        "inspire": 2840545,
        "record": 159544,
        "doi": "10.17182/hepdata.159544.v1",
        "tables": 3,
    },
    "HERMES_2013_I1208547": {
        "inspire": 1208547,
        "record": 62097,
        "doi": "10.17182/hepdata.62097.v1",
        "tables": 64,
    },
    "COMPASS_2018_I1624692": {
        "inspire": 1624692,
        "record": 83542,
        "doi": "10.17182/hepdata.83542.v1",
        "tables": 162,
    },
}
HERMES_ARCHIVE_SHA256 = (
    "e54ae9e96fb417f43c7845e11319977c21b4c0f7349f00ca987658e00b181e1b"
)
HERMES_ARCHIVE_SIZE = 25_937_678


class VendorError(RuntimeError):
    """A source could not be pinned reproducibly."""


def sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(prefix=path.name + ".", dir=path.parent,
                                     delete=False) as stream:
        temporary = Path(stream.name)
        stream.write(payload)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def atomic_json(path: Path, payload: Mapping[str, Any]) -> None:
    atomic_write(
        path,
        (json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False)
         + "\n").encode("utf-8"),
    )


def download(url: str) -> bytes:
    for attempt in range(8):
        request = urllib.request.Request(
            url, headers={"User-Agent": "HerwigPol-reference-data/1.0"}
        )
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                return response.read()
        except urllib.error.HTTPError as exc:
            if exc.code not in {429, 500, 502, 503, 504} or attempt == 7:
                raise
            retry_after = exc.headers.get("Retry-After")
            delay = float(retry_after) if retry_after else min(30.0, 2.0 ** attempt)
            time.sleep(delay)
        except OSError:
            if attempt == 7:
                raise
            time.sleep(min(30.0, 2.0 ** attempt))
    raise VendorError(f"Unreachable retry state for {url}")


def load_record(cache: Path, inspire: int) -> tuple[bytes, dict[str, Any]]:
    cached = cache / f"hepdata-ins{inspire}-v1-record.json"
    payload = cached.read_bytes() if cached.is_file() else download(
        f"https://www.hepdata.net/record/ins{inspire}?format=json"
    )
    try:
        record = json.loads(payload)
    except json.JSONDecodeError as exc:
        raise VendorError(f"Invalid HEPData record for ins{inspire}: {exc}") from exc
    return payload, record


def vendor_hepdata(measurement: str, specification: Mapping[str, Any],
                   cache: Path, workers: int) -> None:
    inspire = int(specification["inspire"])
    record_payload, record = load_record(cache, inspire)
    metadata = record.get("record", {})
    tables = list(record.get("data_tables", []))
    if (
        metadata.get("hepdata_doi") != specification["doi"]
        or int(metadata.get("version", -1)) != 1
        or int(metadata.get("inspire_id", -1)) != inspire
        or len(tables) != int(specification["tables"])
    ):
        raise VendorError(f"HEPData record identity/inventory changed for {measurement}")

    base = ROOT / "data" / "phenomenology" / measurement
    raw = base / "raw"
    atomic_write(raw / "record.json", record_payload)

    def fetch_table(item: tuple[int, Mapping[str, Any]]) -> tuple[int, bytes, str]:
        number, table = item
        cached = cache / f"hepdata-ins{inspire}-t{number:03d}.json"
        url = (
            f"https://www.hepdata.net/record/data/{int(specification['record'])}/"
            f"{int(table['id'])}/1/"
        )
        payload = cached.read_bytes() if cached.is_file() else download(url)
        if not cached.is_file():
            atomic_write(cached, payload)
        return number, payload, url

    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        fetched = list(executor.map(fetch_table, enumerate(tables, start=1)))

    entries: list[dict[str, Any]] = []
    for number, payload, url in sorted(fetched):
        table = json.loads(payload)
        expected = tables[number - 1]
        if (
            table.get("name") != expected.get("name")
            or table.get("doi") != expected.get("doi")
        ):
            raise VendorError(f"HEPData table identity changed for {measurement}/t{number}")
        path = raw / f"t{number:03d}.json"
        atomic_write(path, payload)
        entries.append({
            "number": number,
            "table_id": int(expected["id"]),
            "name": str(table["name"]),
            "doi": str(table["doi"]),
            "path": str(path.relative_to(ROOT)),
            "url": url,
            "sha256": sha256(payload),
            "headers": [str(header.get("name")) for header in table.get("headers", [])],
            "rows": len(table.get("values", [])),
        })
    atomic_json(base / "source-manifest.json", {
        "measurement": measurement,
        "retrieved": RETRIEVED,
        "record_doi": specification["doi"],
        "record_version": 1,
        "record_path": str((raw / "record.json").relative_to(ROOT)),
        "record_sha256": sha256(record_payload),
        "tables": entries,
    })


def copy_checked(source: Path, destination: Path, expected_sha256: str | None = None,
                 expected_size: int | None = None) -> dict[str, Any]:
    payload = source.read_bytes()
    digest = sha256(payload)
    if expected_sha256 is not None and digest != expected_sha256:
        raise VendorError(f"Checksum mismatch for {source}: {digest}")
    if expected_size is not None and len(payload) != expected_size:
        raise VendorError(f"Size mismatch for {source}: {len(payload)}")
    atomic_write(destination, payload)
    return {"path": str(destination.relative_to(ROOT)), "bytes": len(payload),
            "sha256": digest}


def vendor_auxiliary_sources(args: argparse.Namespace) -> None:
    hermes_base = ROOT / "data/phenomenology/HERMES_2013_I1208547"
    archive = copy_checked(
        args.hermes_archive,
        hermes_base / "raw/HERMES-multiplicities-20161105.tar.gz",
        HERMES_ARCHIVE_SHA256,
        HERMES_ARCHIVE_SIZE,
    )
    manifest_path = hermes_base / "source-manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["full_archive"] = {
        **archive,
        "url": (
            "https://web.archive.org/web/20161105011149id_/"
            "http://www-hermes.desy.de/multiplicities/database/data/"
            "HERMES-multiplicities.tar.gz"
        ),
        "retrieved": RETRIEVED,
        "role": "VM-subtracted five-binning 3D numerical authority",
    }
    atomic_json(manifest_path, manifest)

    compass_base = ROOT / "data/phenomenology/COMPASS_2010_I862410"
    source = copy_checked(
        args.compass_2010_source,
        compass_base / "raw/arXiv-1007.4061-source.tar.gz",
    )
    pdf = copy_checked(
        args.compass_2010_pdf,
        compass_base / "raw/arXiv-1007.4061.pdf",
    )
    atomic_json(compass_base / "source-manifest.json", {
        "measurement": "COMPASS_2010_I862410",
        "retrieved": RETRIEVED,
        "numerical_authority": "arXiv:1007.4061 source/PDF",
        "tex_member": "cern-ph-ep_2010-023.tex",
        "source_archive": {**source, "url": "https://arxiv.org/e-print/1007.4061"},
        "paper_pdf": {**pdf, "url": "https://arxiv.org/pdf/1007.4061"},
    })

    ratio_base = ROOT / "data/phenomenology/COMPASS_2020_I1788430"
    ratio_source = copy_checked(
        args.compass_2020_source,
        ratio_base / "raw/arXiv-2003.11791-source.tar.gz",
        "1be86deab9f235c9ee316e43fcb67a8d003082b2903432e46b00197a9bb17aa3",
        94_008,
    )
    ratio_pdf = copy_checked(
        args.compass_2020_pdf,
        ratio_base / "raw/arXiv-2003.11791.pdf",
        "421123edc14883620d2a70b7d9da5306a5c048b9c5f9f4a46503d414a10a48e1",
        437_294,
    )
    inherited = copy_checked(
        args.compass_2018_ratio_source,
        ratio_base / "raw/arXiv-1802.00584-source.tar.gz",
        "f4e639143dfc09b6a4ac894a2abc07033f15c8b0c9cc15b21c75d2f16d961303",
        810_394,
    )
    atomic_json(ratio_base / "source-manifest.json", {
        "measurement": "COMPASS_2020_I1788430",
        "retrieved": RETRIEVED,
        "numerical_authority": "arXiv:2003.11791 paper TeX tables 1--3",
        "source_archive": {
            **ratio_source, "url": "https://arxiv.org/e-print/2003.11791",
            "tex_member": "main.tex",
        },
        "paper_pdf": {
            **ratio_pdf, "url": "https://arxiv.org/pdf/2003.11791",
        },
        "inherited_kaon_selection_source": {
            **inherited, "url": "https://arxiv.org/e-print/1802.00584",
            "tex_member": "paper_v9.tex",
            "role": "source for the kaon cuts explicitly inherited by arXiv:2003.11791",
        },
        "hepdata_audit": {
            "record_url": "https://www.hepdata.net/record/ins1788430",
            "record_http_status": 404,
            "searches": [
                "inspire id 1788430", "arXiv 2003.11791",
                "DOI 10.1016/j.physletb.2020.135600", "exact paper title",
            ],
            "result": "no official HEPData submission found",
            "checked": RETRIEVED,
        },
    })


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-cache", type=Path, required=True)
    parser.add_argument("--hermes-archive", type=Path, required=True)
    parser.add_argument("--compass-2010-source", type=Path, required=True)
    parser.add_argument("--compass-2010-pdf", type=Path, required=True)
    parser.add_argument("--compass-2020-source", type=Path, required=True)
    parser.add_argument("--compass-2020-pdf", type=Path, required=True)
    parser.add_argument("--compass-2018-ratio-source", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=12)
    args = parser.parse_args()
    if args.workers <= 0:
        raise VendorError("--workers must be positive")
    for measurement, specification in HEPDATA.items():
        vendor_hepdata(measurement, specification, args.source_cache, args.workers)
        print(f"vendored {measurement}")
    vendor_auxiliary_sources(args)
    print("vendored HERMES archive and COMPASS 2010/2020 paper sources")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
