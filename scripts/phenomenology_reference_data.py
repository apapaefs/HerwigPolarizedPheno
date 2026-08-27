#!/usr/bin/env python3
"""Pinned reference-data validation and YODA generation for phenomenology.

The normalized JSON files are the offline campaign inputs.  Refreshing data is
deliberately conservative: every downloaded byte stream, table identity,
header, row count, and normalized value must agree with the vendored snapshot.
An upstream revision therefore fails loudly instead of changing a prediction
comparison silently.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import os
import shutil
import tempfile
import urllib.request
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


SCRIPT_PATH = Path(__file__).resolve()
DISPOL_ROOT = SCRIPT_PATH.parents[1]


class ReferenceDataError(RuntimeError):
    """A pinned source or normalized reference failed validation."""


SOURCES: dict[str, list[dict[str, Any]]] = {
    "HERMES_2007_I726689_LEGACY": [
        {
            "url": "https://www.hepdata.net/record/data/11211/10813/1",
            "sha256": "faa0f3fc562ebd8fc96148a57f43b95367ff0b37a1b747c1f740f6b8613412ae",
            "raw": "data/experimental/HERMES_2007_I726689_LEGACY/raw-hepdata-table.json",
            "kind": "hepdata",
            "name": "Table 7",
            "doi": "10.17182/hepdata.11211.v1/t7",
            "headers": ["X", "Y", "Q**2 [GEV**2]", "ASYM(NAME=PARALLEL,C=MEASURED)", "ASYM(NAME=PARALLEL,C=BORN)"],
            "rows": 45,
        }
    ],
    "COMPASS_2010_I843494": [
        {
            "url": "https://www.hepdata.net/record/data/61588/102092/2",
            "sha256": "97514a9cac2987377ff3531533dee5fe6ce58c87f9ffef5ff8abf47ebba06e7a",
            "raw": "data/experimental/COMPASS_2010_I843494/raw-hepdata-table.json",
            "kind": "hepdata",
            "name": "Table 1",
            "doi": "10.17182/hepdata.61588.v2/t1",
            "headers": ["X", "Q**2 [GEV**2]", "A1", " G1"],
            "rows": 15,
        }
    ],
    "STAR_2019_I1708793": [
        {
            "url": "https://www.hepdata.net/record/data/105912/1134045/1",
            "sha256": "fb38647137874114997af0ba4303be0431e92c9a83f4495676b561596cf00818",
            "raw": "data/phenomenology/STAR_2019_I1708793/raw-hepdata-table17.json",
            "kind": "hepdata",
            "name": "Figure 5, $A_L$ for $W^+ \\rightarrow e^+$, combined data samples",
            "doi": "10.17182/hepdata.105912.v1/t17",
            "headers": ["$\\eta_e$", "$A_L$"],
            "rows": 6,
        },
        {
            "url": "https://www.hepdata.net/record/data/105912/1134046/1",
            "sha256": "388aadd7d2eef07b1fbd7c3160904cfd62e335afca8683ca845cf2157e515dd1",
            "raw": "data/phenomenology/STAR_2019_I1708793/raw-hepdata-table18.json",
            "kind": "hepdata",
            "name": "Figure 5, $A_L$ for $W^-\\rightarrow e^-$, combined data samples",
            "doi": "10.17182/hepdata.105912.v1/t18",
            "headers": ["$\\eta_e$", "$A_L$"],
            "rows": 6,
        },
        {
            "url": "https://arxiv.org/pdf/1812.04817",
            "sha256": "d7fde389ec21927b4869f9b8b772b01326284cd33255b0e84ef5e0511e058a7c",
            "kind": "pdf",
        },
    ],
    "PHENIX_2023_I2033856": [
        {
            "url": "https://www.hepdata.net/record/data/129088/1314457/1",
            "sha256": "f5977a346f9c144bc0f96ab9c196424bbbc0e61bb85715e8888d8cbda37657b4",
            "raw": "data/phenomenology/PHENIX_2023_I2033856/raw-hepdata-table1.json",
            "kind": "hepdata",
            "name": "Figure 1",
            "doi": "10.17182/hepdata.129088.v1/t1",
            "headers": ["$p_T$ (GeV/c)", "$Ed^{3}\\sigma^{inc}/dp^{3}$ (pb GeV$^{-2}$ c$^{3}$)", "$Ed^{3}\\sigma^{iso}/dp^{3}$ (pb GeV$^{-2}$ c$^{3}$)"],
            "rows": 18,
        },
        {
            "url": "https://www.hepdata.net/record/data/129088/1314458/1",
            "sha256": "187b02458b4c07c6741de52b2e7c25f9c30a71d401c54ef49369082d75bca57d",
            "raw": "data/phenomenology/PHENIX_2023_I2033856/raw-hepdata-table2.json",
            "kind": "hepdata",
            "name": "Figure 2",
            "doi": "10.17182/hepdata.129088.v1/t2",
            "headers": ["$p_T$ (GeV/c)", "$A_{LL}$"],
            "rows": 7,
        },
    ],
}


REFERENCE_PATHS = {
    "HERMES_2007_I726689_LEGACY": "data/experimental/HERMES_2007_I726689_LEGACY/reference.json",
    "COMPASS_2010_I843494": "data/experimental/COMPASS_2010_I843494/reference.json",
    "STAR_2019_I1708793": "data/phenomenology/STAR_2019_I1708793/reference.json",
    "PHENIX_2023_I2033856": "data/phenomenology/PHENIX_2023_I2033856/reference.json",
}


REFERENCE_YODA_PATHS = {
    "HERMES_2007_I726689_LEGACY": "analyses/rivet/dis/HERMES_2007_I726689_LEGACY.yoda.gz",
    "COMPASS_2010_I843494": "analyses/rivet/dis/COMPASS_2010_I843494.yoda.gz",
    "STAR_2019_I1708793": "analyses/rivet/pp/STAR_2019_I1708793.yoda.gz",
    "PHENIX_2023_I2033856": "analyses/rivet/pp/PHENIX_2023_I2033856.yoda.gz",
}


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ReferenceDataError(f"Could not read {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ReferenceDataError(f"Expected a JSON object in {path}")
    return value


def _table(payload: bytes, source: Mapping[str, Any]) -> dict[str, Any]:
    try:
        table = json.loads(payload.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ReferenceDataError(f"Invalid HEPData JSON from {source['url']}: {exc}") from exc
    headers = [str(item.get("name")) for item in table.get("headers", [])]
    if str(table.get("name")) != source["name"] or str(table.get("doi")) != source["doi"]:
        raise ReferenceDataError(f"HEPData table identity changed at {source['url']}")
    if headers != list(source["headers"]):
        raise ReferenceDataError(f"HEPData schema changed at {source['url']}: {headers}")
    if len(table.get("values", [])) != int(source["rows"]):
        raise ReferenceDataError(f"HEPData row count changed at {source['url']}")
    return table


def _error(cell: Mapping[str, Any], *labels: str) -> float:
    for item in cell.get("errors", []):
        if str(item.get("label")) in labels and "symerror" in item:
            return float(item["symerror"])
    raise ReferenceDataError(f"Missing symmetric uncertainty {labels}")


def _dependent(row: Mapping[str, Any], group: int) -> Mapping[str, Any]:
    for cell in row.get("y", []):
        if int(cell.get("group", -1)) == group:
            return cell
    raise ReferenceDataError(f"Missing dependent-variable group {group}")


def normalize_compass(table: Mapping[str, Any]) -> dict[str, Any]:
    points: list[dict[str, Any]] = []
    edges: list[float] = []
    for index, row in enumerate(table["values"], start=1):
        xcell, qcell = row["x"]
        value = _dependent(row, 0)
        low, high = float(xcell["low"]), float(xcell["high"])
        if not edges:
            edges.append(low)
        edges.append(high)
        systematic = _error(value, "sys", "syst")
        points.append({
            "bin": index, "x_low": low, "x_high": high,
            "x_mean": float(xcell.get("value", 0.5*(low+high))),
            "q2_mean": float(qcell["value"]), "value": float(value["value"]),
            "stat": _error(value, "stat"), "systematic_combined": systematic,
            "systematics": {"published_total": systematic},
        })
    return {
        "schema_version": 2, "measurement": "COMPASS_2010_I843494",
        "observable": "A1p", "selection": "Q2 > 1 GeV2",
        "bin_edges": edges, "points": points,
        "rivet_path": "/COMPASS_2010_I843494/d01-x01-y01",
        "provenance": {
            "source": "HEPData Table 1 v2",
            "source_url": SOURCES["COMPASS_2010_I843494"][0]["url"],
            "source_sha256": SOURCES["COMPASS_2010_I843494"][0]["sha256"],
            "table_doi": "10.17182/hepdata.61588.v2/t1",
            "record_doi": "10.17182/hepdata.61588.v2",
            "retrieved": "2026-07-21",
        },
    }


HERMES_Q2_THRESHOLDS = [
    [], [], [], [], [1.0], [1.0], [1.0], [1.0],
    [1.505, 2.27], [1.62, 2.62], [1.74, 3.02], [1.88, 3.49],
    [2.06, 4.02], [2.23, 4.61], [2.66, 5.5], [3.3, 6.625],
    [4.09, 7.98], [5.04, 9.455], [7.645, 12.5],
]


def normalize_hermes(table: Mapping[str, Any]) -> dict[str, Any]:
    groups: list[dict[str, Any]] = []
    for row in table["values"]:
        xcell, ycell, qcell = row["x"]
        low, high = float(xcell["low"]), float(xcell["high"])
        if not groups or groups[-1]["x_low"] != low:
            groups.append({"x_low": low, "x_high": high, "points": []})
        value = _dependent(row, 1)
        systematic = _error(value, "sys", "syst")
        q2 = float(qcell["value"])
        groups[-1]["points"].append({
            "bin": len(groups[-1]["points"])+1,
            "x_mean": float(xcell.get("value", 0.5*(low+high))),
            "x_mean_source": "published" if "value" in xcell else "bin midpoint",
            "y_mean": float(ycell["value"]), "q2_mean": q2,
            "value": float(value["value"]), "stat": _error(value, "stat"),
            "systematic_combined": systematic,
            "systematics": {"published_total_including_5p2pct_normalization": systematic},
            "nonperturbative_extrapolation": q2 < 1.0,
        })
    for index, group in enumerate(groups):
        edges = [0.18, *HERMES_Q2_THRESHOLDS[index], 20.0]
        if len(group["points"]) != len(edges)-1:
            raise ReferenceDataError(f"HERMES x slice {index+1} changed multiplicity")
        group.update({
            "slice": index+1, "q2_edges": edges,
            "rivet_path": f"/HERMES_2007_I726689_LEGACY/d07-x{index+1:02d}-y02",
        })
        for bin_index, point in enumerate(group["points"]):
            point["q2_low"], point["q2_high"] = edges[bin_index:bin_index+2]
    return {
        "schema_version": 3, "measurement": "HERMES_2007_I726689_LEGACY",
        "observable": "Aparallel Born", "selection": "0.18 < Q2 < 20 GeV2",
        "x_edges": [0.0041, 0.0073, 0.0118, 0.0168, 0.0212, 0.0295, 0.0362,
                    0.0444, 0.0568, 0.0727, 0.0929, 0.119, 0.152, 0.194,
                    0.249, 0.318, 0.406, 0.520, 0.665, 0.9],
        "datasets": groups, "point_count": sum(len(g["points"]) for g in groups),
        "scale_floor_gev": 1.0, "quantitative_mask": "q2_mean < 1 GeV2",
        "provenance": {
            "source": "HEPData Table 7", "source_url": SOURCES["HERMES_2007_I726689_LEGACY"][0]["url"],
            "source_sha256": SOURCES["HERMES_2007_I726689_LEGACY"][0]["sha256"],
            "table_doi": "10.17182/hepdata.11211.v1/t7",
            "publication": "Phys. Rev. D 75 (2007) 012007", "retrieved": "2026-07-21",
            "normalization_uncertainty_percent": 5.2,
            "normalization_included_in_published_systematic": True,
        },
    }


def normalize_star(table17: Mapping[str, Any], table18: Mapping[str, Any]) -> dict[str, Any]:
    channels: dict[str, Any] = {}
    for charge, table, doi in (
        ("Wplus", table17, "10.17182/hepdata.105912.v1/t17"),
        ("Wminus", table18, "10.17182/hepdata.105912.v1/t18"),
    ):
        points = []
        for index, row in enumerate(table["values"], start=1):
            value = _dependent(row, 0)
            systematic = _error(value, "sys", "syst")
            points.append({
                "point": index, "eta": float(row["x"][0]["value"]),
                "value": float(value["value"]), "stat": _error(value, "stat"),
                "systematic_combined": systematic,
                "systematics": {"published_total": systematic},
            })
        channels[charge] = {
            "observable": "AL", "rivet_path": f"/STAR_2019_I1708793/{charge}_AL",
            "eta_bin_edges": [-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5],
            "points": points, "hepdata_table_doi": doi,
        }
    publication = {
        "Wplus": [(0.25, 0.016, 0.042, 0.011), (0.72, 0.072, 0.054, 0.011), (1.24, 0.000, 0.262, 0.028)],
        "Wminus": [(0.27, -0.012, 0.101, 0.019), (0.74, -0.028, 0.092, 0.020), (1.27, -0.147, 0.260, 0.038)],
    }
    for charge, values in publication.items():
        channels[charge]["all"] = {
            "observable": "ALL", "rivet_path": f"/STAR_2019_I1708793/{charge}_ALL",
            "abs_eta_bin_edges": [0.0, 0.5, 1.0, 1.5],
            "points": [{
                "point": index+1, "abs_eta": x, "value": value, "stat": stat,
                "systematic_combined": systematic,
                "systematics": {"published_total": systematic},
            } for index, (x, value, stat, systematic) in enumerate(values)],
        }
    channels["Zgamma"] = {
        "observable": "AL", "rivet_path": "/STAR_2019_I1708793/Zgamma_AL",
        "points": [{"point": 1, "coordinate": 0.0, "value": -0.04, "stat": 0.07,
                    "systematic_combined": 0.0,
                    "systematics": {"published_negligible": 0.0}}],
    }
    return {
        "schema_version": 3, "measurement": "STAR_2019_I1708793", "sqrt_s_gev": 510.0,
        "channels": channels,
        "provenance": {
            "hepdata_record_doi": "10.17182/hepdata.105912.v1",
            "table17_url": SOURCES["STAR_2019_I1708793"][0]["url"],
            "table17_sha256": SOURCES["STAR_2019_I1708793"][0]["sha256"],
            "table18_url": SOURCES["STAR_2019_I1708793"][1]["url"],
            "table18_sha256": SOURCES["STAR_2019_I1708793"][1]["sha256"],
            "publication_arxiv": "1812.04817", "publication_pdf_url": SOURCES["STAR_2019_I1708793"][2]["url"],
            "publication_pdf_sha256": SOURCES["STAR_2019_I1708793"][2]["sha256"],
            "publication_page": 7, "publication_table": "Table I", "z_result_page": 7,
            "retrieved": "2026-07-21",
        },
    }


def normalize_phenix(table1: Mapping[str, Any], table2: Mapping[str, Any]) -> dict[str, Any]:
    edges: list[float] = []
    inclusive: list[dict[str, Any]] = []
    isolated: list[dict[str, Any]] = []
    for index, row in enumerate(table1["values"], start=1):
        xcell = row["x"][0]
        low, high = float(xcell["low"]), float(xcell["high"])
        if not edges:
            edges.append(low)
        edges.append(high)
        for output, group in ((inclusive, 0), (isolated, 1)):
            value = _dependent(row, group)
            systematic = _error(value, "sys", "syst")
            output.append({
                "bin": index, "pt_low": low, "pt_high": high, "pt_mean": 0.5*(low+high),
                "value": float(value["value"]), "stat": _error(value, "stat"),
                "systematic_combined": systematic,
                "systematics": {"point_to_point_total": systematic},
            })
    all_edges: list[float] = []
    all_points: list[dict[str, Any]] = []
    for index, row in enumerate(table2["values"], start=1):
        xcell, value = row["x"][0], _dependent(row, 0)
        low, high = float(xcell["low"]), float(xcell["high"])
        if not all_edges:
            all_edges.append(low)
        all_edges.append(high)
        systematic = _error(value, "sys", "syst")
        all_points.append({
            "bin": index, "pt_low": low, "pt_high": high, "pt_mean": 0.5*(low+high),
            "value": float(value["value"]), "stat": _error(value, "stat"),
            "systematic_combined": systematic,
            "systematics": {"point_to_point_total": systematic},
        })
    return {
        "schema_version": 3, "measurement": "PHENIX_2023_I2033856", "sqrt_s_gev": 510.0,
        "channels": {
            "inclusive_cross_section": {"observable": "E d3sigma/dp3", "units": "pb GeV^-2", "rivet_path": "/PHENIX_2023_I2033856/d01-x01-y01", "bin_edges": edges, "points": inclusive},
            "isolated_cross_section": {"observable": "E d3sigma/dp3", "units": "pb GeV^-2", "rivet_path": "/PHENIX_2023_I2033856/d01-x01-y02", "bin_edges": edges, "points": isolated},
            "isolated_all": {"observable": "ALL", "rivet_path": "/PHENIX_2023_I2033856/d02-x01-y01", "bin_edges": all_edges, "points": all_points},
        },
        "global_uncertainties": {"cross_section_luminosity_relative": 0.10, "all_relative_luminosity_absolute": 3.9e-4, "all_polarization_relative": 0.066},
        "provenance": {
            "hepdata_record_doi": "10.17182/hepdata.129088.v1",
            "table1_url": SOURCES["PHENIX_2023_I2033856"][0]["url"], "table1_sha256": SOURCES["PHENIX_2023_I2033856"][0]["sha256"],
            "table2_url": SOURCES["PHENIX_2023_I2033856"][1]["url"], "table2_sha256": SOURCES["PHENIX_2023_I2033856"][1]["sha256"],
            "retrieved": "2026-07-21",
        },
    }


def normalized_from_raw(measurement: str) -> dict[str, Any]:
    sources = SOURCES[measurement]
    tables = [
        _table((DISPOL_ROOT/source["raw"]).read_bytes(), source)
        for source in sources if source["kind"] == "hepdata"
    ]
    if measurement == "COMPASS_2010_I843494":
        return normalize_compass(tables[0])
    if measurement == "HERMES_2007_I726689_LEGACY":
        return normalize_hermes(tables[0])
    if measurement == "STAR_2019_I1708793":
        return normalize_star(tables[0], tables[1])
    if measurement == "PHENIX_2023_I2033856":
        return normalize_phenix(tables[0], tables[1])
    raise ReferenceDataError(f"Unknown measurement {measurement}")


def validate_vendored(measurement: str) -> dict[str, Any]:
    if measurement not in SOURCES:
        raise ReferenceDataError(f"Unknown measurement {measurement}")
    for source in SOURCES[measurement]:
        raw = source.get("raw")
        if raw:
            payload = (DISPOL_ROOT/raw).read_bytes()
            if sha256_bytes(payload) != source["sha256"]:
                raise ReferenceDataError(f"Vendored raw checksum mismatch for {raw}")
            _table(payload, source)
    generated = normalized_from_raw(measurement)
    snapshot = load_json(DISPOL_ROOT/REFERENCE_PATHS[measurement])
    if generated != snapshot:
        raise ReferenceDataError(f"Normalized snapshot is stale for {measurement}")
    return snapshot


def fetch_and_validate(measurement: str, cache_directory: Path | None = None) -> list[Path]:
    snapshot = validate_vendored(measurement)
    cache = cache_directory or (DISPOL_ROOT/"campaigns"/"phenomenology"/"_data_cache"/measurement)
    cache.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    for index, source in enumerate(SOURCES[measurement], start=1):
        request = urllib.request.Request(str(source["url"]), headers={"User-Agent": "HerwigPol-phenomenology-reference/1.0"})
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                payload = response.read()
        except OSError as exc:
            raise ReferenceDataError(f"Could not download {source['url']}: {exc}") from exc
        actual = sha256_bytes(payload)
        if actual != source["sha256"]:
            raise ReferenceDataError(f"Pinned checksum changed at {source['url']}: expected {source['sha256']}, got {actual}")
        if source["kind"] == "hepdata":
            _table(payload, source)
            raw = DISPOL_ROOT/source["raw"]
            if raw.read_bytes() != payload:
                raise ReferenceDataError(f"Checksum-valid source differs from vendored bytes: {raw}")
            suffix = ".json"
        else:
            suffix = ".pdf"
        destination = cache/f"source-{index:02d}{suffix}"
        destination.write_bytes(payload)
        outputs.append(destination)
    write_reference_yoda(measurement, snapshot)
    return outputs


def _import_yoda() -> Any:
    try:
        import yoda  # type: ignore
    except (ImportError, OSError) as exc:
        raise ReferenceDataError("YODA Python bindings are unavailable; load herwig/pol") from exc
    return yoda


def _write_yoda(yoda: Any, objects: Sequence[Any], destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="phenodata-", dir=destination.parent) as temporary:
        plain = Path(temporary)/destination.with_suffix("").name
        yoda.write(list(objects), str(plain))
        if destination.suffix == ".gz":
            compressed = Path(temporary)/destination.name
            with plain.open("rb") as source, compressed.open("wb") as target:
                with gzip.GzipFile(filename="", mode="wb", fileobj=target, mtime=0) as stream:
                    shutil.copyfileobj(source, stream)
            os.replace(compressed, destination)
        else:
            os.replace(plain, destination)


def _reference_path(path: str) -> str:
    """Return Rivet's canonical reference-data namespace for an AO path."""
    return path if path.startswith("/REF/") else "/REF" + path


def _estimate(yoda: Any, edges: Sequence[float], points: Sequence[Mapping[str, Any]], path: str, annotations: Mapping[str, Any] | None = None) -> Any:
    result = yoda.BinnedEstimate1D(
        [float(edge) for edge in edges], _reference_path(path)
    )
    result.setAnnotation("IsRef", 1)
    for key, value in (annotations or {}).items():
        result.setAnnotation(str(key), value)
    for index, point in enumerate(points, start=1):
        bin_object = result.bin(index)
        bin_object.setVal(float(point["value"]))
        bin_object.setErr(-float(point["stat"]), float(point["stat"]), "stat")
        for label, error in point.get("systematics", {}).items():
            bin_object.setErr(-float(error), float(error), str(label))
    return result


def _scatter(yoda: Any, points: Sequence[Mapping[str, Any]], coordinate: str, path: str, total: bool) -> Any:
    result = yoda.Scatter2D(_reference_path(path))
    result.setAnnotation("IsRef", 1)
    result.setAnnotation("ErrorDisplay", "total" if total else "statistical only")
    for point in points:
        error = float(point["stat"])
        if total:
            error = math.hypot(error, float(point["systematic_combined"]))
        result.addPoint(float(point[coordinate]), float(point["value"]), 0.0, error)
    return result


def write_reference_yoda(measurement: str, snapshot: Mapping[str, Any] | None = None) -> Path:
    snapshot = snapshot or validate_vendored(measurement)
    yoda = _import_yoda()
    objects: list[Any] = []
    if measurement == "HERMES_2007_I726689_LEGACY":
        for dataset in snapshot["datasets"]:
            objects.append(_estimate(yoda, dataset["q2_edges"], dataset["points"], dataset["rivet_path"], {
                "Observable": "Born A_parallel", "XLow": dataset["x_low"], "XHigh": dataset["x_high"],
                "NonperturbativeMask": json.dumps([bool(point["nonperturbative_extrapolation"]) for point in dataset["points"]]),
            }))
    elif measurement == "COMPASS_2010_I843494":
        objects.append(_estimate(
            yoda,
            snapshot["bin_edges"],
            snapshot["points"],
            snapshot["rivet_path"],
            {
                "Observable": "A1p",
                "A2G2Assumption": "eta*A2 = 0; g2 is not modeled",
            },
        ))
    elif measurement == "STAR_2019_I1708793":
        for charge in ("Wplus", "Wminus"):
            channel = snapshot["channels"][charge]
            objects.append(_scatter(yoda, channel["points"], "eta", channel["rivet_path"], True))
            objects.append(_scatter(yoda, channel["points"], "eta", channel["rivet_path"]+"_StatOnly", False))
            objects.append(_scatter(yoda, channel["all"]["points"], "abs_eta", channel["all"]["rivet_path"], True))
            objects.append(_scatter(yoda, channel["all"]["points"], "abs_eta", channel["all"]["rivet_path"]+"_StatOnly", False))
        z = snapshot["channels"]["Zgamma"]
        objects.append(_scatter(yoda, z["points"], "coordinate", z["rivet_path"], True))
        objects.append(_scatter(yoda, z["points"], "coordinate", z["rivet_path"]+"_StatOnly", False))
    elif measurement == "PHENIX_2023_I2033856":
        for channel in snapshot["channels"].values():
            objects.append(_estimate(yoda, channel["bin_edges"], channel["points"], channel["rivet_path"], {
                "Observable": channel["observable"], "GlobalUncertainties": json.dumps(snapshot["global_uncertainties"], sort_keys=True),
            }))
    else:
        raise ReferenceDataError(f"No YODA writer for {measurement}")
    destination = DISPOL_ROOT/REFERENCE_YODA_PATHS[measurement]
    _write_yoda(yoda, objects, destination)
    return destination


# The STAR jet records and the APS HERMES SIDIS archive use a richer source
# model (cross-observable covariance blocks and a version-pinned ZIP,
# respectively).  Keep their parser isolated while exposing one stable public
# API to the campaign runner and tests.
import polarized_jet_sidis_reference_data as _jet_sidis_reference
import compass_sidis_reference_data as _compass_sidis_reference

for _measurement in _jet_sidis_reference.MEASUREMENTS:
    SOURCES.setdefault(_measurement, [])
for _measurement in _compass_sidis_reference.MEASUREMENTS:
    SOURCES.setdefault(_measurement, [])
REFERENCE_PATHS.update(_jet_sidis_reference.REFERENCE_PATHS)
REFERENCE_YODA_PATHS.update(_jet_sidis_reference.REFERENCE_YODA_PATHS)
REFERENCE_PATHS.update(_compass_sidis_reference.REFERENCE_PATHS)
REFERENCE_YODA_PATHS.update(_compass_sidis_reference.REFERENCE_YODA_PATHS)

_legacy_normalized_from_raw = normalized_from_raw
_legacy_validate_vendored = validate_vendored
_legacy_fetch_and_validate = fetch_and_validate
_legacy_write_reference_yoda = write_reference_yoda


def _translate_new_reference_error(function: Any, *args: Any, **kwargs: Any) -> Any:
    try:
        return function(*args, **kwargs)
    except _jet_sidis_reference.NewReferenceDataError as exc:
        raise ReferenceDataError(str(exc)) from exc


def _translate_compass_reference_error(
    function: Any, *args: Any, **kwargs: Any
) -> Any:
    try:
        return function(*args, **kwargs)
    except _compass_sidis_reference.CompassReferenceDataError as exc:
        raise ReferenceDataError(str(exc)) from exc


def normalized_from_raw(measurement: str) -> dict[str, Any]:
    if measurement in _compass_sidis_reference.MEASUREMENTS:
        return _translate_compass_reference_error(
            _compass_sidis_reference.normalized_from_raw, measurement
        )
    if measurement in _jet_sidis_reference.MEASUREMENTS:
        return _translate_new_reference_error(
            _jet_sidis_reference.normalized_from_raw, measurement
        )
    return _legacy_normalized_from_raw(measurement)


def validate_vendored(measurement: str) -> dict[str, Any]:
    if measurement in _compass_sidis_reference.MEASUREMENTS:
        return _translate_compass_reference_error(
            _compass_sidis_reference.validate_vendored, measurement
        )
    if measurement in _jet_sidis_reference.MEASUREMENTS:
        return _translate_new_reference_error(
            _jet_sidis_reference.validate_vendored, measurement
        )
    return _legacy_validate_vendored(measurement)


def fetch_and_validate(
    measurement: str,
    cache_directory: Path | None = None,
    source_file: Path | None = None,
) -> list[Path]:
    if measurement in _compass_sidis_reference.MEASUREMENTS:
        if source_file is not None:
            raise ReferenceDataError(
                "--source-file is not used for version-pinned HEPData records"
            )
        return _translate_compass_reference_error(
            _compass_sidis_reference.fetch_and_validate,
            measurement,
            cache_directory,
        )
    if measurement in _jet_sidis_reference.MEASUREMENTS:
        return _translate_new_reference_error(
            _jet_sidis_reference.fetch_and_validate,
            measurement,
            cache_directory,
            source_file,
        )
    if source_file is not None:
        raise ReferenceDataError(
            "--source-file is only supported for the APS HERMES SIDIS archive"
        )
    return _legacy_fetch_and_validate(measurement, cache_directory)


def write_reference_yoda(
    measurement: str, snapshot: Mapping[str, Any] | None = None
) -> Path:
    if measurement in _compass_sidis_reference.MEASUREMENTS:
        return _translate_compass_reference_error(
            _compass_sidis_reference.write_reference_yoda,
            measurement,
            snapshot,
        )
    if measurement in _jet_sidis_reference.MEASUREMENTS:
        return _translate_new_reference_error(
            _jet_sidis_reference.write_reference_yoda, measurement, snapshot
        )
    return _legacy_write_reference_yoda(measurement, snapshot)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("measurement", choices=sorted(SOURCES))
    parser.add_argument("--fetch", action="store_true")
    parser.add_argument("--write-yoda", action="store_true")
    parser.add_argument("--source-file", type=Path)
    args = parser.parse_args()
    validated = validate_vendored(args.measurement)
    if args.fetch:
        fetch_and_validate(args.measurement, source_file=args.source_file)
    elif args.write_yoda:
        write_reference_yoda(args.measurement, validated)
    print(f"validated {args.measurement}")
