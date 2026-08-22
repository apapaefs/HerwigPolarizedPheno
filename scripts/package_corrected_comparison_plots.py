#!/usr/bin/env python3
"""Package corrected data overlays as vector PDF and 600-dpi PNG."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import struct
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence


REPOSITORY = Path(__file__).resolve().parents[1]


class PackageError(RuntimeError):
    """Raised when a curated plot cannot be identified or verified."""


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise PackageError(f"Could not read {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise PackageError(f"Expected a JSON object in {path}")
    return value


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def campaign_dir(spec: Mapping[str, Any]) -> Path:
    return (
        REPOSITORY / "campaigns" / str(spec["storage"])
        / str(spec["measurement"]) / str(spec["tag"])
    )


def one_pdf(root: Path, stem: str) -> Path:
    matches = sorted(path for path in root.rglob(f"{stem}.pdf") if path.is_file())
    if len(matches) != 1:
        raise PackageError(
            f"Expected exactly one {stem}.pdf below {root}, found {len(matches)}"
        )
    if matches[0].stat().st_size == 0:
        raise PackageError(f"Empty source PDF {matches[0]}")
    return matches[0]


def png_dimensions(path: Path) -> tuple[int, int]:
    with path.open("rb") as stream:
        header = stream.read(24)
    if len(header) != 24 or header[:8] != b"\x89PNG\r\n\x1a\n":
        raise PackageError(f"Not a valid PNG header: {path}")
    return struct.unpack(">II", header[16:24])


def render_png(pdf: Path, destination: Path) -> tuple[int, int]:
    executable = shutil.which("pdftocairo")
    if not executable:
        raise PackageError("pdftocairo is required for the 600-dpi export")
    output_stem = destination.with_suffix("")
    completed = subprocess.run(
        [executable, "-png", "-r", "600", "-singlefile", str(pdf), str(output_stem)],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )
    if completed.returncode != 0:
        raise PackageError(f"pdftocairo failed for {pdf}: {completed.stdout}")
    if not destination.is_file() or destination.stat().st_size == 0:
        raise PackageError(f"No rendered PNG was produced for {pdf}")
    width, height = png_dimensions(destination)
    if min(width, height) < 2000:
        raise PackageError(
            f"600-dpi render is unexpectedly small ({width}x{height}): {destination}"
        )
    return width, height


def fixed_plot_specs(config: Mapping[str, Any]) -> list[dict[str, Any]]:
    names = {
        "HERMES_2007_I726689": "hermes-a1",
        "COMPASS_2010_I843494": "compass-2010-a1",
        "COMPASS_2016_I1357198": "compass-2016-a1",
        "COMPASS_2017_I1501480": "compass-2017-a1",
    }
    output = []
    for spec in config["fixed"]:
        measurement = str(spec["measurement"])
        if measurement not in names:
            continue
        snapshot = load_json(
            REPOSITORY / "data" / "experimental" / measurement / "reference.json"
        )
        output.append(
            {
                "spec": spec, "stem": Path(str(snapshot["rivet_path"])).name,
                "destination": Path("dis") / names[measurement],
            }
        )
    return output


def sidis_plot_specs(config: Mapping[str, Any]) -> list[dict[str, Any]]:
    spec = next(
        item for item in config["fixed"]
        if item["measurement"] == "HERMES_2019_I1698889"
    )
    snapshot = load_json(
        REPOSITORY / "data" / "phenomenology" / "HERMES_2019_I1698889"
        / "reference.json"
    )
    output = []
    for identifier, dataset in sorted(snapshot["datasets"].items()):
        if dataset.get("observable") != "A_parallel":
            continue
        output.append(
            {
                "spec": spec, "stem": Path(str(dataset["rivet_path"])).name,
                "destination": Path("hermes-sidis") / identifier,
            }
        )
    if len(output) != 24:
        raise PackageError(f"Expected 24 HERMES SIDIS primary overlays, found {len(output)}")
    return output


def star_plot_specs(config: Mapping[str, Any]) -> list[dict[str, Any]]:
    spec = config["star510"][0]
    snapshot = load_json(
        REPOSITORY / "data" / "phenomenology" / "STAR_2022_I1949588"
        / "reference.json"
    )
    output = [
        {
            "spec": spec, "stem": Path(str(dataset["rivet_path"])).name,
            "destination": Path("star-jets-510") / f"{identifier}_ALL",
        }
        for identifier, dataset in sorted(snapshot["datasets"].items())
    ]
    if len(output) != 5:
        raise PackageError(f"Expected five STAR 510 primary overlays, found {len(output)}")
    return output


def package(
    config_path: Path,
    output_root: Path,
    source_commit: str,
    *,
    selection: str = "all",
) -> Path:
    config = load_json(config_path)
    if selection == "all":
        entries = (
            fixed_plot_specs(config)
            + sidis_plot_specs(config)
            + star_plot_specs(config)
        )
        expected = 33
        destination_name = (
            f"compatibility-corrected-20260819-{source_commit[:12]}-600dpi"
        )
        selection_record = {
            "inclusive_fixed_target": 4,
            "hermes_sidis": 24,
            "star_510": 5,
            "total": 33,
            "withheld": [
                "HERMES_2007_I726689_LEGACY",
                "PHENIX_2023_I2033856",
            ],
        }
        description = (
            "four corrected inclusive fixed-target overlays, 24 HERMES SIDIS "
            "overlays, and five STAR 510 GeV overlays"
        )
    elif selection == "star510":
        entries = star_plot_specs(config)
        expected = 5
        destination_name = (
            f"star510-pt13p1-bloch-20260822-{source_commit[:12]}-600dpi"
        )
        selection_record = {
            "inclusive_fixed_target": 0,
            "hermes_sidis": 0,
            "star_510": 5,
            "total": 5,
            "inclusive_primary_bins": "5-14 (analysis pT >= 13.1 GeV)",
            "inclusive_diagnostic_only_bins": "1-4",
            "dijet_policy": "all published bins remain primary",
            "spin_density_policy": "radial_bloch_ball_projection",
            "withheld": [
                "HERMES_2007_I726689_LEGACY",
                "PHENIX_2023_I2033856",
            ],
        }
        description = (
            "five STAR 510 GeV overlays using inclusive analysis bins 5--14 "
            "(pT >= 13.1 GeV) and all dijet bins"
        )
    else:
        raise PackageError(f"Unknown plot selection {selection!r}")
    if len(entries) != expected:
        raise PackageError(
            f"Expected {expected} primary overlays for {selection}, found {len(entries)}"
        )
    destination_root = output_root / destination_name
    if destination_root.exists():
        raise PackageError(
            f"Refusing to overwrite existing curated package {destination_root}"
        )
    destination_root.mkdir(parents=True)
    provenance_entries = []
    campaign_records: dict[str, Any] = {}
    for entry in entries:
        spec = entry["spec"]
        root = campaign_dir(spec)
        manifest_path = root / "manifest.json"
        manifest = load_json(manifest_path)
        if manifest.get("status") != "complete" or manifest.get("plots") is None:
            raise PackageError(f"Campaign is not completely plotted: {root}")
        source = one_pdf(root / "plots", str(entry["stem"]))
        destination_base = destination_root / entry["destination"]
        destination_base.parent.mkdir(parents=True, exist_ok=True)
        pdf = destination_base.with_suffix(".pdf")
        png = destination_base.with_suffix(".png")
        shutil.copy2(source, pdf)
        width, height = render_png(pdf, png)
        key = f"{spec['measurement']}/{spec['tag']}"
        campaign_records.setdefault(
            key,
            {
                "manifest": str(manifest_path.resolve()),
                "manifest_sha256": sha256_file(manifest_path),
                "summary_sha256": sha256_file(root / "postprocess" / "summary.json"),
                "source_commit": (
                    manifest.get("runtime", {}).get("provenance", {})
                    .get("source_control", {}).get("commit")
                ),
                "plugin_sha256": manifest.get("plugin_provenance", {}).get("sha256"),
            },
        )
        provenance_entries.append(
            {
                "measurement": spec["measurement"], "tag": spec["tag"],
                "source_pdf": str(source.resolve()),
                "pdf": str(pdf.relative_to(destination_root)),
                "pdf_sha256": sha256_file(pdf),
                "png": str(png.relative_to(destination_root)),
                "png_sha256": sha256_file(png),
                "png_pixels": [width, height], "png_dpi": 600,
            }
        )
    if any(record["source_commit"] != source_commit for record in campaign_records.values()):
        raise PackageError("A packaged campaign was not generated from the requested commit")
    payload = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "source_commit": source_commit,
        "selection": selection_record,
        "campaigns": campaign_records,
        "plots": provenance_entries,
    }
    (destination_root / "PROVENANCE.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (destination_root / "README.md").write_text(
        "# Corrected experimental comparisons\n\n"
        f"This bundle contains {description}. Each vector "
        "PDF has a 600-dpi PNG sibling. HERMES low-Q2 and PHENIX are deliberately "
        "withheld. Exact campaign, runtime, plugin, and file hashes are recorded "
        "in `PROVENANCE.json`.\n",
        encoding="utf-8",
    )
    print(f"Packaged {expected} corrected overlays at {destination_root}")
    return destination_root


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign-config", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--selection", choices=("all", "star510"), default="all")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = make_parser().parse_args(argv)
    try:
        package(
            args.campaign_config,
            args.output_root,
            args.source_commit,
            selection=args.selection,
        )
    except PackageError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
