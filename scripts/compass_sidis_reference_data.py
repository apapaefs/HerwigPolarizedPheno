#!/usr/bin/env python3
"""Version-pinned COMPASS SIDIS reference data and Rivet exports.

HEPData v1 is the numerical authority for all three measurements.  The 2009
publication table is retained only as a frozen, checksum-attributed audit; no
paper value is substituted into the normalized reference snapshot.
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
ROOT = SCRIPT_PATH.parents[1]


class CompassReferenceDataError(RuntimeError):
    """A COMPASS source, audit, or normalized snapshot failed validation."""


MEASUREMENTS = {
    "COMPASS_2009_I820721",
    "COMPASS_2017_I1444985",
    "COMPASS_2017_I1483098",
}

REFERENCE_PATHS = {
    identifier: f"data/phenomenology/{identifier}/reference.json"
    for identifier in MEASUREMENTS
}
SOURCE_MANIFEST_PATHS = {
    identifier: f"data/phenomenology/{identifier}/source-manifest.json"
    for identifier in MEASUREMENTS
}
REFERENCE_YODA_PATHS = {
    identifier: f"analyses/rivet/dis/{identifier}.yoda.gz"
    for identifier in MEASUREMENTS
}

AUDIT_PATH = (
    "data/phenomenology/COMPASS_2009_I820721/"
    "paper-hepdata-discrepancy-audit.json"
)

SPECIFICATIONS: dict[str, dict[str, Any]] = {
    "COMPASS_2009_I820721": {
        "record_doi": "10.17182/hepdata.55300.v1",
        "record_id": 55300,
        "inspire_id": 820721,
        "table_count": 4,
        "species": {
            "piplus": {"table": 2, "group": 0, "pid": 211},
            "piminus": {"table": 2, "group": 1, "pid": -211},
            "kplus": {"table": 2, "group": 2, "pid": 321},
            "kminus": {"table": 2, "group": 3, "pid": -321},
        },
    },
    "COMPASS_2017_I1444985": {
        "record_doi": "10.17182/hepdata.76800.v1",
        "record_id": 76800,
        "inspire_id": 1444985,
        "table_count": 4,
        "species": {
            "piplus": {"table": 1, "pid": 211, "pid_mode": "identified"},
            "piminus": {"table": 2, "pid": -211, "pid_mode": "identified"},
            "hplus": {"table": 3, "charge": 1, "pid_mode": "stable_charged_hadron"},
            "hminus": {"table": 4, "charge": -1, "pid_mode": "stable_charged_hadron"},
        },
    },
    "COMPASS_2017_I1483098": {
        "record_doi": "10.17182/hepdata.77892.v1",
        "record_id": 77892,
        "inspire_id": 1483098,
        "table_count": 2,
        "species": {
            "kplus": {"table": 1, "pid": 321, "pid_mode": "identified"},
            "kminus": {"table": 2, "pid": -321, "pid_mode": "identified"},
        },
    },
}

X_EDGES_2009 = [0.004, 0.006, 0.01, 0.02, 0.03, 0.04, 0.06, 0.1, 0.15, 0.2, 0.3]
X_EDGES_2017 = [0.004, 0.01, 0.02, 0.03, 0.04, 0.06, 0.1, 0.14, 0.18, 0.4]
Y_EDGES_2017 = [0.1, 0.15, 0.2, 0.3, 0.5, 0.7]
Z_EDGES_2017 = [0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.85]
SPECIES_ORDER_2009 = ["piplus", "piminus", "kplus", "kminus"]


# Literal transcription of tables_withRC2.tex, table tab:SIDIS_asym, from the
# checksum-pinned arXiv source archive.  Values remain separate from HEPData.
PAPER_2009_ROWS: list[dict[str, Any]] = [
    {"x_mean": 0.0052, "q2_mean": 1.16, "piplus": [0.006, 0.014, 0.006], "piminus": [0.009, 0.014, 0.006], "kplus": [-0.018, 0.029, 0.011], "kminus": [-0.084, 0.030, 0.014]},
    {"x_mean": 0.0079, "q2_mean": 1.42, "piplus": [-0.003, 0.008, 0.003], "piminus": [-0.002, 0.008, 0.003], "kplus": [-0.017, 0.017, 0.007], "kminus": [0.014, 0.018, 0.007]},
    {"x_mean": 0.0141, "q2_mean": 2.03, "piplus": [-0.003, 0.007, 0.003], "piminus": [-0.007, 0.007, 0.003], "kplus": [-0.039, 0.014, 0.006], "kminus": [-0.004, 0.016, 0.006]},
    {"x_mean": 0.0244, "q2_mean": 3.19, "piplus": [-0.001, 0.011, 0.004], "piminus": [0.006, 0.012, 0.005], "kplus": [0.020, 0.022, 0.009], "kminus": [0.025, 0.026, 0.010]},
    {"x_mean": 0.0346, "q2_mean": 4.43, "piplus": [0.026, 0.015, 0.006], "piminus": [0.004, 0.016, 0.006], "kplus": [0.021, 0.030, 0.012], "kminus": [0.012, 0.035, 0.014]},
    {"x_mean": 0.0487, "q2_mean": 6.10, "piplus": [0.016, 0.015, 0.006], "piminus": [0.037, 0.016, 0.007], "kplus": [0.066, 0.029, 0.011], "kminus": [-0.058, 0.035, 0.015]},
    {"x_mean": 0.0763, "q2_mean": 9.26, "piplus": [0.046, 0.017, 0.008], "piminus": [0.018, 0.018, 0.007], "kplus": [0.064, 0.032, 0.013], "kminus": [0.015, 0.041, 0.017]},
    {"x_mean": 0.121, "q2_mean": 14.9, "piplus": [0.094, 0.025, 0.012], "piminus": [0.087, 0.028, 0.013], "kplus": [0.117, 0.046, 0.019], "kminus": [-0.007, 0.065, 0.026]},
    {"x_mean": 0.171, "q2_mean": 22.4, "piplus": [0.102, 0.039, 0.017], "piminus": [0.132, 0.044, 0.020], "kplus": [0.116, 0.070, 0.028], "kminus": [0.002, 0.103, 0.041]},
    {"x_mean": 0.240, "q2_mean": 32.8, "piplus": [0.218, 0.044, 0.024], "piminus": [0.147, 0.051, 0.023], "kplus": [0.208, 0.078, 0.031], "kminus": [-0.018, 0.118, 0.047]},
]


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise CompassReferenceDataError(f"Could not read {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise CompassReferenceDataError(f"Expected a JSON object in {path}")
    return value


def _manifest(measurement: str) -> dict[str, Any]:
    path = ROOT / SOURCE_MANIFEST_PATHS[measurement]
    value = _json(path)
    specification = SPECIFICATIONS[measurement]
    if value.get("measurement") != measurement:
        raise CompassReferenceDataError(f"Source manifest identity mismatch in {path}")
    if value.get("record_doi") != specification["record_doi"] or int(value.get("record_version", -1)) != 1:
        raise CompassReferenceDataError(f"Pinned HEPData version changed for {measurement}")
    return value


def _validate_table(payload: bytes, expected: Mapping[str, Any]) -> dict[str, Any]:
    if sha256_bytes(payload) != expected["sha256"]:
        raise CompassReferenceDataError(f"Vendored checksum mismatch for {expected['path']}")
    try:
        table = json.loads(payload.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise CompassReferenceDataError(f"Invalid HEPData table JSON: {exc}") from exc
    headers = [str(item.get("name")) for item in table.get("headers", [])]
    if str(table.get("name")) != expected["name"] or str(table.get("doi")) != expected["doi"]:
        raise CompassReferenceDataError(f"HEPData table identity changed for {expected['path']}")
    if headers != list(expected["headers"]):
        raise CompassReferenceDataError(f"HEPData headers changed for {expected['path']}")
    if len(table.get("values", [])) != int(expected["rows"]):
        raise CompassReferenceDataError(f"HEPData row count changed for {expected['path']}")
    return table


def _tables(measurement: str) -> dict[int, dict[str, Any]]:
    manifest = _manifest(measurement)
    specification = SPECIFICATIONS[measurement]
    record_path = ROOT / str(manifest["record_path"])
    record_payload = record_path.read_bytes()
    if sha256_bytes(record_payload) != str(manifest["record_sha256"]):
        raise CompassReferenceDataError(f"Vendored HEPData record checksum mismatch for {measurement}")
    record = json.loads(record_payload)
    metadata = record.get("record", {})
    if metadata.get("hepdata_doi") != specification["record_doi"] or int(metadata.get("version", -1)) != 1:
        raise CompassReferenceDataError(f"HEPData record identity changed for {measurement}")
    entries = manifest.get("tables", [])
    if len(entries) != int(specification["table_count"]):
        raise CompassReferenceDataError(f"Expected complete table inventory for {measurement}")
    expected_numbers = set(range(1, int(specification["table_count"]) + 1))
    if {int(entry["number"]) for entry in entries} != expected_numbers:
        raise CompassReferenceDataError(f"HEPData table inventory changed for {measurement}")
    if {str(item.get("doi")) for item in record.get("data_tables", [])} != {
        str(entry["doi"]) for entry in entries
    }:
        raise CompassReferenceDataError(f"HEPData record/table DOI inventory changed for {measurement}")
    return {
        int(entry["number"]): _validate_table(
            (ROOT / str(entry["path"])).read_bytes(), entry
        )
        for entry in entries
    }


def _dependent(row: Mapping[str, Any], group: int) -> Mapping[str, Any]:
    for item in row.get("y", []):
        if int(item.get("group", -1)) == group:
            return item
    raise CompassReferenceDataError(f"Missing dependent-variable group {group}")


def _error(cell: Mapping[str, Any], label: str) -> float:
    for item in cell.get("errors", []):
        if str(item.get("label", "")).lower() in {label, "syst" if label == "sys" else label}:
            return float(item["symerror"])
    raise CompassReferenceDataError(f"Missing {label} uncertainty")


def _bin_index(low: float, high: float, edges: Sequence[float]) -> int:
    for index, (left, right) in enumerate(zip(edges, edges[1:])):
        if math.isclose(low, left, rel_tol=0.0, abs_tol=1.0e-12) and math.isclose(high, right, rel_tol=0.0, abs_tol=1.0e-12):
            return index
    raise CompassReferenceDataError(f"Cell [{low}, {high}] is outside the pinned binning")


def paper_hepdata_audit(table: Mapping[str, Any]) -> dict[str, Any]:
    hepdata_rows: list[dict[str, Any]] = []
    for row in table["values"]:
        entry: dict[str, Any] = {
            "x_mean": float(row["x"][0]["value"]),
            "q2_mean": float(row["x"][1]["value"]),
        }
        for species, group in zip(SPECIES_ORDER_2009, range(4)):
            cell = _dependent(row, group)
            entry[species] = [
                float(cell["value"]), _error(cell, "stat"), _error(cell, "sys")
            ]
        hepdata_rows.append(entry)
    if len(hepdata_rows) != len(PAPER_2009_ROWS):
        raise CompassReferenceDataError("2009 paper/HEPData audit row count changed")
    discrepancies: list[dict[str, Any]] = []
    labels = ["a1", "stat", "systematic"]
    for index, (paper, hepdata) in enumerate(zip(PAPER_2009_ROWS, hepdata_rows), start=1):
        for field in ("x_mean", "q2_mean"):
            if float(paper[field]) != float(hepdata[field]):
                discrepancies.append({
                    "bin": index, "field": field,
                    "paper": float(paper[field]), "hepdata": float(hepdata[field]),
                    "delta_hepdata_minus_paper": float(hepdata[field]) - float(paper[field]),
                })
        for species in SPECIES_ORDER_2009:
            for offset, label in enumerate(labels):
                paper_value = float(paper[species][offset])
                hepdata_value = float(hepdata[species][offset])
                if paper_value != hepdata_value:
                    discrepancies.append({
                        "bin": index, "field": f"{species}.{label}",
                        "paper": paper_value, "hepdata": hepdata_value,
                        "delta_hepdata_minus_paper": hepdata_value - paper_value,
                    })
    manifest = _manifest("COMPASS_2009_I820721")
    return {
        "schema_version": 1,
        "measurement": "COMPASS_2009_I820721",
        "numerical_authority": "HEPData v1",
        "policy": "Differences are audited; paper values are never substituted.",
        "paper_source": manifest["paper_audit_source"],
        "paper_rows": PAPER_2009_ROWS,
        "hepdata_rows": hepdata_rows,
        "compared_fields": [
            "x_mean", "q2_mean",
            *[
                f"{species}.{field}"
                for species in SPECIES_ORDER_2009
                for field in labels
            ],
        ],
        "discrepancy_count": len(discrepancies),
        "discrepancies": discrepancies,
    }


def _correlation_blocks(table: Mapping[str, Any], datasets: Mapping[str, Any]) -> tuple[list[dict[str, Any]], list[list[float]], list[str]]:
    rows = table["values"]
    if len(rows) != 100:
        raise CompassReferenceDataError("COMPASS 2009 correlation table must contain 100 rows")
    pair_blocks = {
        2: ("piminus", "piplus"),
        4: ("kplus", "piplus"),
        5: ("kplus", "piminus"),
        7: ("kminus", "piplus"),
        8: ("kminus", "piminus"),
        9: ("kminus", "kplus"),
    }
    x_means = [float(point["x_mean"]) for point in datasets["piplus"]["points"]]
    blocks: list[dict[str, Any]] = []
    order: list[str] = []
    full = [[0.0 for _ in range(40)] for _ in range(40)]
    for bin_index in range(10):
        correlation = [
            [1.0 if first == second else 0.0 for second in SPECIES_ORDER_2009]
            for first in SPECIES_ORDER_2009
        ]
        for block_index, (first, second) in pair_blocks.items():
            row = rows[block_index * 10 + bin_index]
            if not math.isclose(float(row["x"][0]["value"]), x_means[bin_index], rel_tol=0.0, abs_tol=5.0e-4):
                raise CompassReferenceDataError("COMPASS 2009 correlation-row order changed")
            value = float(_dependent(row, 0)["value"])
            i, j = SPECIES_ORDER_2009.index(first), SPECIES_ORDER_2009.index(second)
            correlation[i][j] = correlation[j][i] = value
        stat = [float(datasets[species]["points"][bin_index]["stat"]) for species in SPECIES_ORDER_2009]
        covariance = [
            [correlation[i][j] * stat[i] * stat[j] for j in range(4)]
            for i in range(4)
        ]
        for i, species in enumerate(SPECIES_ORDER_2009):
            order.append(f"{species}:x{bin_index + 1:02d}")
            for j in range(4):
                full[4 * bin_index + i][4 * bin_index + j] = covariance[i][j]
        blocks.append({
            "x_bin": bin_index + 1,
            "x_low": X_EDGES_2009[bin_index],
            "x_high": X_EDGES_2009[bin_index + 1],
            "order": SPECIES_ORDER_2009,
            "correlation": correlation,
            "covariance": covariance,
        })
    return blocks, full, order


def normalize_2009() -> dict[str, Any]:
    tables = _tables("COMPASS_2009_I820721")
    table = tables[2]
    datasets: dict[str, Any] = {}
    for species, metadata in SPECIFICATIONS["COMPASS_2009_I820721"]["species"].items():
        points: list[dict[str, Any]] = []
        for index, row in enumerate(table["values"]):
            xcell, qcell = row["x"]
            cell = _dependent(row, int(metadata["group"]))
            systematic = _error(cell, "sys")
            value = float(cell["value"])
            correlated = 0.08 * value
            uncorrelated = math.sqrt(max(0.0, systematic * systematic - correlated * correlated))
            points.append({
                "bin": index + 1,
                "flat_bin": index,
                "x_low": float(xcell["low"]),
                "x_high": float(xcell["high"]),
                "x_mean": float(xcell["value"]),
                "q2_mean": float(qcell["value"]),
                "a1": value,
                "stat": _error(cell, "stat"),
                "systematic": systematic,
                "systematic_correlated_8pct": correlated,
                "systematic_uncorrelated": uncorrelated,
            })
        datasets[species] = {
            "id": species,
            "pid": metadata["pid"],
            "observable": "A1d",
            "rivet_path": f"/COMPASS_2009_I820721/A1_{species}",
            "raw_objects": {
                "ordinary": f"Yield_{species}_x",
                "inverse_depolarization": f"YieldOverD_{species}_x",
                "covariance": f"CovarianceProxy_{species}_x",
            },
            "points": points,
        }
    if any(
        [point["x_low"] for point in dataset["points"]] != X_EDGES_2009[:-1]
        or [point["x_high"] for point in dataset["points"]] != X_EDGES_2009[1:]
        for dataset in datasets.values()
    ):
        raise CompassReferenceDataError("COMPASS 2009 x-bin edges changed")
    blocks, statistical_covariance, point_order = _correlation_blocks(tables[3], datasets)
    audit = paper_hepdata_audit(table)
    return {
        "schema_version": 1,
        "measurement": "COMPASS_2009_I820721",
        "observable": "identified-hadron A1d(x)",
        "beam_energy_gev": 160.0,
        "binning": {"x": X_EDGES_2009},
        "selection": {
            "q2_min_gev2": 1.0, "x": [0.004, 0.3], "y": [0.1, 0.9],
            "z": [0.2, 0.85], "hadron_momentum_gev": [10.0, 50.0],
        },
        "datasets": datasets,
        "point_order": point_order,
        "statistical_correlation_blocks": blocks,
        "statistical_covariance": statistical_covariance,
        "systematics": {
            "correlated_multiplicative_fraction": 0.08,
            "correlated_scope": "one nuisance shared by all 40 SIDIS points",
            "remaining_variance": "diagonal max(sigma_syst^2-(0.08*A1)^2,0)",
        },
        "target_model": {
            "sigma_uu": "(p+n)/2",
            "sigma_ll": "0.925*(p+n)/2",
            "deuteron_d_state_factor": 0.925,
        },
        "depolarization": {
            "factor": "finite-muon-mass COMPASS D with R1998",
            "a2_assumption": "eta*A2=0",
        },
        "paper_hepdata_audit": {
            "path": AUDIT_PATH,
            "discrepancy_count": audit["discrepancy_count"],
            "numerical_authority": "HEPData v1",
        },
        "provenance": {
            "record_doi": SPECIFICATIONS["COMPASS_2009_I820721"]["record_doi"],
            "source_manifest": SOURCE_MANIFEST_PATHS["COMPASS_2009_I820721"],
            "tables": [1, 2, 3, 4],
            "retrieved": _manifest("COMPASS_2009_I820721")["retrieved"],
        },
    }


def _multiplicity_cell(row: Mapping[str, Any], index: int) -> dict[str, Any]:
    xcell, ycell, qcell, zcell = row["x"]
    value = _dependent(row, 0)
    systematic = _error(value, "sys")
    return {
        "cell": index + 1,
        "flat_bin": index,
        "x_low": float(xcell["low"]), "x_high": float(xcell["high"]),
        "x_mean": float(xcell.get("value", 0.5 * (float(xcell["low"]) + float(xcell["high"])))),
        "y_low": float(ycell["low"]), "y_high": float(ycell["high"]),
        "y_mean": float(ycell.get("value", 0.5 * (float(ycell["low"]) + float(ycell["high"])))),
        "q2_mean": float(qcell["value"]),
        "z_low": float(zcell["low"]), "z_high": float(zcell["high"]),
        "z_mean": float(zcell.get("value", 0.5 * (float(zcell["low"]) + float(zcell["high"])))),
        "value": float(value["value"]),
        "stat": _error(value, "stat"),
        "systematic": systematic,
        "systematic_correlated_80pct": 0.8 * systematic,
        "systematic_uncorrelated_60pct": 0.6 * systematic,
        "corrections": {
            "dvm_hadron": float(_dependent(row, 1)["value"]),
            "dvm_dis": float(_dependent(row, 2)["value"]),
            "radiative_hadron": float(_dependent(row, 3)["value"]),
            "radiative_dis": float(_dependent(row, 4)["value"]),
        },
    }


def _cell_signature(point: Mapping[str, Any]) -> tuple[float, ...]:
    return tuple(float(point[key]) for key in (
        "x_low", "x_high", "y_low", "y_high", "z_low", "z_high"
    ))


def _slices(measurement: str, species: str, points: list[dict[str, Any]]) -> list[dict[str, Any]]:
    slices: list[dict[str, Any]] = []
    for point in points:
        key = (point["x_low"], point["x_high"], point["y_low"], point["y_high"])
        if not slices or tuple(slices[-1][name] for name in ("x_low", "x_high", "y_low", "y_high")) != key:
            ix = _bin_index(point["x_low"], point["x_high"], X_EDGES_2017)
            iy = _bin_index(point["y_low"], point["y_high"], Y_EDGES_2017)
            slices.append({
                "slice": len(slices) + 1,
                "x_bin": ix + 1, "y_bin": iy + 1,
                "x_low": point["x_low"], "x_high": point["x_high"],
                "y_low": point["y_low"], "y_high": point["y_high"],
                "flat_bins": [], "z_edges": [],
                "rivet_path": f"/{measurement}/{species}_x{ix + 1:02d}_y{iy + 1:02d}",
            })
        current = slices[-1]
        if not current["z_edges"]:
            current["z_edges"].append(point["z_low"])
        if not math.isclose(current["z_edges"][-1], point["z_low"], rel_tol=0.0, abs_tol=1.0e-12):
            raise CompassReferenceDataError(f"Non-contiguous z cells in {measurement}/{species}")
        current["z_edges"].append(point["z_high"])
        current["flat_bins"].append(point["flat_bin"])
        point["slice"] = current["slice"]
        point["slice_bin"] = len(current["flat_bins"])
    return slices


def normalize_multiplicity(measurement: str) -> dict[str, Any]:
    tables = _tables(measurement)
    specification = SPECIFICATIONS[measurement]
    datasets: dict[str, Any] = {}
    signatures: list[list[tuple[float, ...]]] = []
    for species, metadata in specification["species"].items():
        table_number = int(metadata["table"])
        points = [
            _multiplicity_cell(row, index)
            for index, row in enumerate(tables[table_number]["values"])
        ]
        signatures.append([_cell_signature(point) for point in points])
        slices = _slices(measurement, species, points)
        datasets[species] = {
            "id": species,
            "observable": "dM/dz",
            "pid_mode": metadata["pid_mode"],
            **({"pid": metadata["pid"]} if "pid" in metadata else {"charge": metadata["charge"]}),
            "table": table_number,
            "table_doi": tables[table_number]["doi"],
            "flat_rivet_path": f"/{measurement}/Multiplicity_{species}_cells",
            "raw_objects": {
                "numerator": f"HadronNumerator_{species}_cells",
                "denominator": f"DISDenominator_{species}_cells",
                "covariance": f"CovarianceProxy_{species}_cells",
            },
            "points": points,
            "slices": slices,
        }
    if any(signature != signatures[0] for signature in signatures[1:]):
        raise CompassReferenceDataError(f"Species valid-cell maps differ in {measurement}")
    if len(datasets[next(iter(datasets))]["slices"]) != 38:
        raise CompassReferenceDataError(f"Expected 38 published (x,y) slices for {measurement}")
    cell_count = len(signatures[0])
    expected = 311 if measurement == "COMPASS_2017_I1444985" else 309
    if cell_count != expected:
        raise CompassReferenceDataError(f"Expected {expected} valid cells for {measurement}")
    return {
        "schema_version": 1,
        "measurement": measurement,
        "observable": "charged-hadron multiplicity dM(x,y,z)/dz",
        "beam_energy_gev": 160.0,
        "binning": {"x": X_EDGES_2017, "y": Y_EDGES_2017, "z": Z_EDGES_2017},
        "selection": {
            "q2_min_gev2": 1.0, "w_min_gev": 5.0,
            "x": [0.004, 0.4], "y": [0.1, 0.7], "z": [0.2, 0.85],
            "hadron_momentum_gev": [12.0, 40.0],
            "hadron_lab_angle_mrad": [10.0, 120.0],
        },
        "cell_count": cell_count,
        "slice_count": 38,
        "valid_cell_map": [list(item) for item in signatures[0]],
        "datasets": datasets,
        "target_model": {
            "isoscalar_numerator": "(p+n)/2",
            "isoscalar_dis_denominator": "(p+n)/2",
            "ratio_order": "sum P/N before numerator/DIS ratio and divide by dz",
        },
        "systematics": {
            "correlated_fraction_of_published_systematic": 0.8,
            "uncorrelated_fraction_of_published_systematic": 0.6,
            "correlated_scope": "one nuisance shared by every table and cell in the record",
            "statistical_covariance": "diagonal; none published",
        },
        "corrections": {
            "reference": "final radiative- and diffractive-vector-meson-corrected multiplicity",
            "herwig": "no detector, radiative, or experimental correction applied",
            "provenance_fields": ["dvm_hadron", "dvm_dis", "radiative_hadron", "radiative_dis"],
        },
        "unidentified_hadron_z": (
            "pion-mass assumption"
            if measurement == "COMPASS_2017_I1444985" else "not applicable"
        ),
        "provenance": {
            "record_doi": specification["record_doi"],
            "source_manifest": SOURCE_MANIFEST_PATHS[measurement],
            "tables": list(range(1, int(specification["table_count"]) + 1)),
            "retrieved": _manifest(measurement)["retrieved"],
        },
    }


def normalized_from_raw(measurement: str) -> dict[str, Any]:
    if measurement == "COMPASS_2009_I820721":
        return normalize_2009()
    if measurement in {"COMPASS_2017_I1444985", "COMPASS_2017_I1483098"}:
        return normalize_multiplicity(measurement)
    raise CompassReferenceDataError(f"Unknown COMPASS SIDIS measurement {measurement}")


def _assert_equal(generated: Any, snapshot: Any, path: str = "$") -> None:
    if type(generated) is not type(snapshot):
        raise CompassReferenceDataError(f"Snapshot type mismatch at {path}")
    if isinstance(generated, Mapping):
        if set(generated) != set(snapshot):
            raise CompassReferenceDataError(f"Snapshot key mismatch at {path}: {sorted(set(generated) ^ set(snapshot))}")
        for key in generated:
            _assert_equal(generated[key], snapshot[key], f"{path}.{key}")
        return
    if isinstance(generated, list):
        if len(generated) != len(snapshot):
            raise CompassReferenceDataError(f"Snapshot length mismatch at {path}")
        for index, (actual, expected) in enumerate(zip(generated, snapshot)):
            _assert_equal(actual, expected, f"{path}[{index}]")
        return
    if isinstance(generated, float):
        if not math.isclose(generated, snapshot, rel_tol=1.0e-15, abs_tol=1.0e-24):
            raise CompassReferenceDataError(f"Snapshot float mismatch at {path}: {generated!r} != {snapshot!r}")
        return
    if generated != snapshot:
        raise CompassReferenceDataError(f"Snapshot mismatch at {path}: {generated!r} != {snapshot!r}")


def _validate_psd(matrix: Sequence[Sequence[float]], label: str) -> None:
    try:
        import numpy as np
    except ImportError as exc:
        raise CompassReferenceDataError("NumPy is required for covariance validation") from exc
    array = np.asarray(matrix, dtype=float)
    if not np.allclose(array, array.T, rtol=0.0, atol=1.0e-12):
        raise CompassReferenceDataError(f"{label} is not symmetric")
    values = np.linalg.eigvalsh(array)
    tolerance = max(1.0, float(np.max(np.abs(values)))) * 1.0e-10
    if float(values[0]) < -tolerance:
        raise CompassReferenceDataError(f"{label} is not positive semidefinite")


def validate_vendored(measurement: str) -> dict[str, Any]:
    if measurement not in MEASUREMENTS:
        raise CompassReferenceDataError(f"Unknown COMPASS SIDIS measurement {measurement}")
    generated = normalized_from_raw(measurement)
    snapshot = _json(ROOT / REFERENCE_PATHS[measurement])
    _assert_equal(generated, snapshot)
    if measurement == "COMPASS_2009_I820721":
        audit = paper_hepdata_audit(_tables(measurement)[2])
        _assert_equal(audit, _json(ROOT / AUDIT_PATH))
        _validate_psd(snapshot["statistical_covariance"], "COMPASS 2009 statistical covariance")
        for block in snapshot["statistical_correlation_blocks"]:
            _validate_psd(block["correlation"], f"COMPASS 2009 x bin {block['x_bin']} correlation")
            _validate_psd(block["covariance"], f"COMPASS 2009 x bin {block['x_bin']} covariance")
    return snapshot


def fetch_and_validate(measurement: str, cache_directory: Path | None = None) -> list[Path]:
    validate_vendored(measurement)
    cache = cache_directory or ROOT / "campaigns" / "phenomenology" / "_data_cache" / measurement
    cache.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    for entry in _manifest(measurement)["tables"]:
        request = urllib.request.Request(str(entry["url"]), headers={"User-Agent": "HerwigPol-reference-data/1.0"})
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                payload = response.read()
        except OSError as exc:
            raise CompassReferenceDataError(f"Could not download {entry['url']}: {exc}") from exc
        _validate_table(payload, entry)
        if payload != (ROOT / str(entry["path"])).read_bytes():
            raise CompassReferenceDataError(f"Checksum-valid HEPData table differs from {entry['path']}")
        destination = cache / f"table-{int(entry['number']):02d}.json"
        destination.write_bytes(payload)
        outputs.append(destination)
    write_reference_yoda(measurement)
    return outputs


def _import_yoda() -> Any:
    try:
        import yoda  # type: ignore
    except (ImportError, OSError) as exc:
        raise CompassReferenceDataError("YODA Python bindings are unavailable; load herwig/pol") from exc
    return yoda


def _estimate(yoda: Any, edges: Sequence[float], points: Sequence[Mapping[str, Any]], path: str, value_key: str) -> Any:
    result = yoda.BinnedEstimate1D([float(edge) for edge in edges], "/REF" + path)
    result.setAnnotation("IsRef", 1)
    for index, point in enumerate(points, start=1):
        target = result.bin(index)
        target.setVal(float(point[value_key]))
        target.setErr(-float(point["stat"]), float(point["stat"]), "stat")
        correlated_key = "systematic_correlated_8pct" if value_key == "a1" else "systematic_correlated_80pct"
        uncorrelated_key = "systematic_uncorrelated" if value_key == "a1" else "systematic_uncorrelated_60pct"
        target.setErr(-abs(float(point[correlated_key])), abs(float(point[correlated_key])), "syst_correlated")
        target.setErr(-float(point[uncorrelated_key]), float(point[uncorrelated_key]), "syst_uncorrelated")
    return result


def _write_yoda(yoda: Any, objects: Iterable[Any], destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="compass-sidis-", dir=destination.parent) as temporary:
        plain = Path(temporary) / destination.with_suffix("").name
        yoda.write(list(objects), str(plain))
        compressed = Path(temporary) / destination.name
        with plain.open("rb") as source, compressed.open("wb") as target:
            with gzip.GzipFile(filename="", mode="wb", fileobj=target, mtime=0) as stream:
                shutil.copyfileobj(source, stream)
        os.replace(compressed, destination)


def write_reference_yoda(measurement: str, snapshot: Mapping[str, Any] | None = None) -> Path:
    snapshot = snapshot or validate_vendored(measurement)
    yoda = _import_yoda()
    objects: list[Any] = []
    if measurement == "COMPASS_2009_I820721":
        for dataset in snapshot["datasets"].values():
            item = _estimate(yoda, snapshot["binning"]["x"], dataset["points"], dataset["rivet_path"], "a1")
            item.setAnnotation("Observable", "A1d")
            item.setAnnotation("NumericalAuthority", "HEPData v1")
            objects.append(item)
    else:
        for dataset in snapshot["datasets"].values():
            flat = _estimate(
                yoda, list(range(len(dataset["points"]) + 1)),
                dataset["points"], dataset["flat_rivet_path"], "value"
            )
            flat.setAnnotation("Observable", "dM/dz")
            flat.setAnnotation("Layout", "published sparse flattened cells")
            objects.append(flat)
            for item in dataset["slices"]:
                points = [dataset["points"][int(index)] for index in item["flat_bins"]]
                sliced = _estimate(yoda, item["z_edges"], points, item["rivet_path"], "value")
                sliced.setAnnotation("Observable", "dM/dz")
                sliced.setAnnotation("XBin", item["x_bin"])
                sliced.setAnnotation("YBin", item["y_bin"])
                objects.append(sliced)
    destination = ROOT / REFERENCE_YODA_PATHS[measurement]
    _write_yoda(yoda, objects, destination)
    return destination


def _atomic_json(path: Path, payload: Mapping[str, Any]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    data = (json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n").encode("utf-8")
    with tempfile.NamedTemporaryFile(prefix=path.name + ".", dir=path.parent, delete=False) as stream:
        temporary = Path(stream.name)
        stream.write(data)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)
    return path


def write_normalized_snapshot(measurement: str) -> Path:
    snapshot = normalized_from_raw(measurement)
    destination = _atomic_json(ROOT / REFERENCE_PATHS[measurement], snapshot)
    if measurement == "COMPASS_2009_I820721":
        _atomic_json(ROOT / AUDIT_PATH, paper_hepdata_audit(_tables(measurement)[2]))
    return destination


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("measurement", choices=sorted(MEASUREMENTS))
    parser.add_argument("--write-snapshot", action="store_true")
    parser.add_argument("--write-yoda", action="store_true")
    parser.add_argument("--fetch", action="store_true")
    args = parser.parse_args()
    if args.write_snapshot:
        write_normalized_snapshot(args.measurement)
    snapshot = validate_vendored(args.measurement)
    if args.fetch:
        fetch_and_validate(args.measurement)
    elif args.write_yoda:
        write_reference_yoda(args.measurement, snapshot)
    print(f"validated {args.measurement}")
