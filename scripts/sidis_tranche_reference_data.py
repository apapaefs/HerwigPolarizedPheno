#!/usr/bin/env python3
"""Pinned reference data for the 2026 SIDIS implementation tranches.

The module turns four complete HEPData-v1 submissions, the full HERMES
multiplicity archive, and the COMPASS proton-asymmetry paper source into one
deterministic offline contract.  Raw inputs are immutable: normalisation may
write the tracked JSON snapshots, audits and reference YODA, but validation
never repairs or rewrites a source.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import math
import os
import re
import shutil
import tarfile
import tempfile
import urllib.request
from decimal import Decimal
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]


class SIDISTrancheDataError(RuntimeError):
    """A pinned SIDIS source, snapshot, or audit failed validation."""


MEASUREMENTS = {
    "COMPASS_2026_I3096394",
    "COMPASS_2025_I2840545",
    "COMPASS_2010_I862410",
    "HERMES_2013_I1208547",
    "COMPASS_2018_I1624692",
}

HEPDATA_MEASUREMENTS = MEASUREMENTS - {"COMPASS_2010_I862410"}
REFERENCE_PATHS = {
    item: f"data/phenomenology/{item}/reference.json" for item in MEASUREMENTS
}
SOURCE_MANIFEST_PATHS = {
    item: f"data/phenomenology/{item}/source-manifest.json" for item in MEASUREMENTS
}
REFERENCE_YODA_PATHS = {
    item: f"analyses/rivet/dis/{item}.yoda.gz" for item in MEASUREMENTS
}
AUDIT_PATHS = {
    "COMPASS_2010_I862410": (
        "data/phenomenology/COMPASS_2010_I862410/paper-extraction-audit.json"
    ),
    "HERMES_2013_I1208547": (
        "data/phenomenology/HERMES_2013_I1208547/"
        "archive-hepdata-projection-audit.json"
    ),
}

EXPECTED = {
    "COMPASS_2026_I3096394": {
        "doi": "10.17182/hepdata.169859.v1", "tables": 3,
        "rows": [306, 306, 305],
    },
    "COMPASS_2025_I2840545": {
        "doi": "10.17182/hepdata.159544.v1", "tables": 3,
        "rows": [302, 302, 298],
    },
    "HERMES_2013_I1208547": {
        "doi": "10.17182/hepdata.62097.v1", "tables": 64,
        "archive_bytes": 25_937_678,
        "archive_sha256": (
            "e54ae9e96fb417f43c7845e11319977c21b4c0f7349f00ca987658e00b181e1b"
        ),
    },
    "COMPASS_2018_I1624692": {
        "doi": "10.17182/hepdata.83542.v1", "tables": 162,
        # The complete v1 table inventory contains 2 x 2332 numerical rows.
        # The earlier assessment's 4918 figure is retained as a documented
        # source-count discrepancy; no absent cells are fabricated.
        "cells": 4664,
    },
}

SPECIES = ("piplus", "piminus", "kplus", "kminus")
SPECIES_PID = {"piplus": 211, "piminus": -211, "kplus": 321, "kminus": -321}
PAIR_TABLES = {
    1: ("hplus", "hminus"),
    2: ("piplus", "piminus"),
    3: ("kplus", "kminus"),
}

X_EDGES_2010 = [
    0.004, 0.006, 0.010, 0.020, 0.030, 0.040, 0.060,
    0.100, 0.150, 0.200, 0.300, 0.400, 0.700,
]
PAPER_2010_ROWS = [
    [0.0052, 1.16, 0.008, .029, .016, 0.020, .029, .016, 0.078, .067, .038, -0.112, .069, .039],
    [0.0079, 1.46, 0.041, .018, .010, 0.016, .018, .010, 0.126, .036, .021, -0.040, .039, .022],
    [0.0142, 2.12, 0.040, .014, .008, 0.049, .015, .009, 0.046, .028, .016, 0.038, .031, .018],
    [0.0245, 3.22, 0.122, .022, .014, 0.055, .023, .013, 0.117, .041, .024, 0.092, .048, .028],
    [0.0346, 4.36, 0.156, .030, .019, 0.060, .032, .018, 0.196, .054, .033, 0.074, .066, .037],
    [0.0487, 5.97, 0.141, .029, .018, 0.118, .031, .019, 0.174, .051, .031, 0.027, .064, .036],
    [0.0765, 8.96, 0.230, .031, .022, 0.053, .033, .019, 0.215, .054, .033, 0.029, .071, .040],
    [0.121, 13.8, 0.243, .041, .027, 0.096, .047, .027, 0.315, .072, .044, 0.212, .101, .058],
    [0.172, 19.6, 0.392, .058, .040, 0.165, .066, .038, 0.355, .099, .059, 0.195, .147, .083],
    [0.240, 27.6, 0.518, .060, .046, 0.233, .069, .041, 0.450, .101, .063, 0.264, .157, .089],
    [0.341, 40.1, 0.549, .097, .064, 0.134, .113, .064, 0.512, .163, .097, 0.375, .259, .147],
    [0.480, 55.6, 0.871, .122, .086, 0.520, .142, .085, 0.726, .207, .124, 0.654, .339, .194],
]

# The six correlations needed for the 4x4 SIDIS-only block, in the order
# (pi-,pi+), (K+,pi+), (K+,pi-), (K-,pi+), (K-,pi-), (K-,K+).
PAPER_2010_CORRELATIONS = [
    [.12,.15,.17,.16,.15,.16,.16,.15,.16,.16,.19,.20],
    [-.17,-.09,-.04,-.02,-.02,-.01,-.02,-.01,-.01,-.02,-.02,-.01],
    [.03,.04,.04,.05,.05,.05,.05,.06,.05,.05,.04,.03],
    [.03,.03,.04,.04,.04,.04,.03,.04,.02,.03,.02,.05],
    [-.16,-.09,-.05,-.03,-.03,-.03,-.03,-.03,-.02,-.03,-.04,-.02],
    [.05,.08,.10,.10,.10,.11,.11,.12,.11,.11,.13,.16],
]
PAPER_2010_DEUTERON_CORRECTIONS = [
    [.001,0,.001,0,0], [.001,0,.001,0,.001], [.001,.001,.002,0,.001],
    [.002,.001,.002,.001,.001], [.002,.001,.003,.001,.002],
    [.003,.001,.003,.001,.002], [.004,.002,.005,.002,.003],
    [.006,.002,.006,.003,.004], [.008,.003,.008,.004,.006],
    [.011,.004,.010,.005,.008], [.015,.005,.013,.009,.011],
    [.020,.006,.017,.013,.015],
]

HERMES_CONFIGS: dict[str, dict[str, Any]] = {
    "z-3D": {
        "axis": "x", "axis_edges": [.023,.085,.6],
        "z_edges": [.1,.15,.2,.25,.3,.4,.5,.6,.7,.8,1.1],
        "phperp_edges": [0,.1,.3,.45,.6,1.2], "usable": 80,
    },
    "zpt-3D": {
        "axis": "x", "axis_edges": [.023,.085,.6],
        "z_edges": [.1,.2,.3,.4,.6,.8,1.1],
        "phperp_edges": [0,.1,.2,.3,.4,.5,.6,.7,.8,1.2], "usable": 90,
    },
    "zx-3D": {
        "axis": "x", "axis_edges": [.023,.04,.055,.075,.1,.14,.2,.3,.4,.6],
        "z_edges": [.1,.2,.3,.4,.6,.8,1.1],
        "phperp_edges": [0,.3,.5,.7,1.2], "usable": 180,
    },
    "zQ2-3D": {
        "axis": "q2", "axis_edges": [1,1.25,1.5,1.75,2,2.25,2.5,3,5,15],
        "x_edges": [.023,.6], "z_edges": [.1,.2,.3,.4,.6,.8,1.1],
        "phperp_edges": [0,.3,.5,.7,1.2], "usable": 180,
    },
    "zxpt-3D": {
        "axis": "x", "axis_edges": [.023,.047,.075,.12,.2,.35,.6],
        "z_edges": [.1,.2,.25,.3,.375,.475,.6,.8,1.1],
        "phperp_edges": [0,.15,.25,.35,.45,.6,.8,1.2], "usable": 294,
    },
}


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SIDISTrancheDataError(f"Could not read {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise SIDISTrancheDataError(f"Expected a JSON object in {path}")
    return value


def _manifest(measurement: str) -> dict[str, Any]:
    value = _load_json(ROOT / SOURCE_MANIFEST_PATHS[measurement])
    if value.get("measurement") != measurement:
        raise SIDISTrancheDataError(f"Source manifest identity mismatch for {measurement}")
    return value


def _tables(measurement: str) -> dict[int, dict[str, Any]]:
    if measurement not in HEPDATA_MEASUREMENTS:
        raise SIDISTrancheDataError(f"{measurement} is not a HEPData source")
    manifest = _manifest(measurement)
    expected = EXPECTED[measurement]
    if manifest.get("record_doi") != expected["doi"] or int(manifest.get("record_version", -1)) != 1:
        raise SIDISTrancheDataError(f"Pinned HEPData version changed for {measurement}")
    record_payload = (ROOT / str(manifest["record_path"])).read_bytes()
    if _sha256(record_payload) != manifest["record_sha256"]:
        raise SIDISTrancheDataError(f"Record checksum mismatch for {measurement}")
    record = json.loads(record_payload)
    entries = manifest.get("tables", [])
    count = int(expected["tables"])
    if len(entries) != count or {int(e["number"]) for e in entries} != set(range(1, count + 1)):
        raise SIDISTrancheDataError(f"Incomplete HEPData inventory for {measurement}")
    if len(record.get("data_tables", [])) != count:
        raise SIDISTrancheDataError(f"HEPData record table count changed for {measurement}")
    output: dict[int, dict[str, Any]] = {}
    for entry in entries:
        payload = (ROOT / str(entry["path"])).read_bytes()
        if _sha256(payload) != entry["sha256"]:
            raise SIDISTrancheDataError(f"Table checksum mismatch: {entry['path']}")
        table = json.loads(payload)
        headers = [str(item.get("name")) for item in table.get("headers", [])]
        if (table.get("name") != entry["name"] or table.get("doi") != entry["doi"]
                or headers != entry["headers"] or len(table.get("values", [])) != int(entry["rows"])):
            raise SIDISTrancheDataError(f"Table identity/schema changed: {entry['path']}")
        output[int(entry["number"])] = table
    return output


def _dependent(row: Mapping[str, Any], group: int) -> Mapping[str, Any]:
    for item in row.get("y", []):
        if int(item.get("group", -1)) == group:
            return item
    raise SIDISTrancheDataError(f"Missing dependent-variable group {group}")


def _error(cell: Mapping[str, Any], label: str) -> float:
    aliases = {label, "sys" if label == "syst" else label, "syst" if label == "sys" else label}
    for item in cell.get("errors", []):
        if str(item.get("label", "")).lower() in aliases and "symerror" in item:
            return float(item["symerror"])
    raise SIDISTrancheDataError(f"Missing {label} uncertainty")


def _atomic_json(path: Path, payload: Mapping[str, Any]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    data = (json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False,
                       allow_nan=False) + "\n").encode("utf-8")
    with tempfile.NamedTemporaryFile(prefix=path.name + ".", dir=path.parent,
                                     delete=False) as stream:
        temporary = Path(stream.name)
        stream.write(data)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)
    return path


def _cell_signature(point: Mapping[str, Any], dimensions: Sequence[str]) -> list[float]:
    result: list[float] = []
    for dimension in dimensions:
        result.extend([float(point[f"{dimension}_low"]), float(point[f"{dimension}_high"])])
    return result


def _slices(measurement: str, dataset: str, points: list[dict[str, Any]],
            slice_dimensions: Sequence[str], plotted_dimension: str) -> list[dict[str, Any]]:
    groups: list[dict[str, Any]] = []
    previous: tuple[float, ...] | None = None
    for point in points:
        key = tuple(float(point[f"{axis}_{side}"])
                    for axis in slice_dimensions for side in ("low", "high"))
        if key != previous:
            groups.append({
                "slice": len(groups) + 1,
                "slice_dimensions": list(slice_dimensions),
                "bounds": {axis: [point[f"{axis}_low"], point[f"{axis}_high"]]
                           for axis in slice_dimensions},
                "flat_bins": [], "edges": [],
                "rivet_path": f"/{measurement}/{dataset}_slice{len(groups) + 1:04d}",
            })
            previous = key
        group = groups[-1]
        low, high = float(point[f"{plotted_dimension}_low"]), float(point[f"{plotted_dimension}_high"])
        if not group["edges"]:
            group["edges"].append(low)
        if not math.isclose(float(group["edges"][-1]), low, abs_tol=1.e-12, rel_tol=0):
            # Sparse HEPData cells can have gaps.  Start a new readable slice
            # rather than fabricating an intervening prediction bin.
            groups.append({
                "slice": len(groups) + 1,
                "slice_dimensions": list(slice_dimensions),
                "bounds": {axis: [point[f"{axis}_low"], point[f"{axis}_high"]]
                           for axis in slice_dimensions},
                "flat_bins": [], "edges": [low],
                "rivet_path": f"/{measurement}/{dataset}_slice{len(groups) + 1:04d}",
            })
            group = groups[-1]
        group["edges"].append(high)
        group["flat_bins"].append(int(point["flat_bin"]))
        point["slice"] = int(group["slice"])
        point["slice_bin"] = len(group["flat_bins"])
    return groups


def normalize_modern_compass(measurement: str) -> dict[str, Any]:
    tables = _tables(measurement)
    expected_rows = EXPECTED[measurement]["rows"]
    datasets: dict[str, Any] = {}
    for table_number, pair in PAIR_TABLES.items():
        table = tables[table_number]
        if len(table["values"]) != expected_rows[table_number - 1]:
            raise SIDISTrancheDataError(f"Modern COMPASS row count changed in table {table_number}")
        headers = [str(item["name"]) for item in table["headers"]]
        for sign_index, species in enumerate(pair):
            points: list[dict[str, Any]] = []
            for flat_bin, row in enumerate(table["values"]):
                xcell, ycell, zcell = row["x"]
                cell = _dependent(row, sign_index)
                systematic = _error(cell, "sys")
                corrections = {
                    re.sub(r"[^a-z0-9]+", "_", headers[group + 3].lower()).strip("_"):
                    float(_dependent(row, group + 2)["value"])
                    for group in range(5)
                }
                points.append({
                    "cell": flat_bin + 1, "flat_bin": flat_bin,
                    "x_low": float(xcell["low"]), "x_high": float(xcell["high"]),
                    "y_low": float(ycell["low"]), "y_high": float(ycell["high"]),
                    "z_low": float(zcell["low"]), "z_high": float(zcell["high"]),
                    "value": float(cell["value"]), "stat": _error(cell, "stat"),
                    "systematic": systematic,
                    "systematic_correlated_80pct": .8 * systematic,
                    "systematic_uncorrelated_60pct": .6 * systematic,
                    "corrections": corrections,
                })
            slices = _slices(measurement, species, points, ("x", "y"), "z")
            mode: dict[str, Any]
            if species.startswith("h"):
                mode = {"pid_mode": "stable_charged_hadron", "charge": 1 if species.endswith("plus") else -1,
                        "z_mass_assumption": "pion"}
            else:
                mode = {"pid_mode": "identified", "pid": SPECIES_PID[species]}
            datasets[species] = {
                "id": species, "table": table_number, "table_doi": table["doi"],
                "observable": "dM/dz", "density_widths": ["z"], **mode,
                "flat_rivet_path": f"/{measurement}/Multiplicity_{species}_cells",
                "raw_objects": {
                    "numerator": f"HadronNumerator_{species}_cells",
                    "denominator": f"DISDenominator_{species}_cells",
                    "covariance": f"CovarianceProxy_{species}_cells",
                },
                "valid_cell_map": [_cell_signature(point, ("x", "y", "z")) for point in points],
                "points": points, "slices": slices,
            }
    target = "H" if measurement == "COMPASS_2025_I2840545" else "D"
    return {
        "schema_version": 2, "measurement": measurement,
        "observable": "charged-hadron multiplicity dM(x,y,z)/dz",
        "beam": {"pid": -13, "energy_gev": 160.0},
        "selection": {
            "q2_min_gev2": 1.0, "w_min_gev": 5.0, "x": [.004,.4],
            "y": [.1,.7], "z": [.2,.85], "hadron_momentum_gev": [12.,40.],
            "hadron_lab_angle_mrad": [10.,120.],
        },
        "datasets": datasets,
        "target_outputs": {target: {"P": 1.0} if target == "H" else {"P": .5, "N": .5}},
        "systematics": {
            "model": "diag((0.6*syst)^2) + (0.8*syst)(0.8*syst)^T",
            "correlated_scope": "record-wide", "statistical_covariance": "diagonal; not released",
        },
        "corrections": {
            "reference": "final radiative- and diffractive-vector-meson-corrected multiplicity",
            "herwig": "no detector or experimental correction applied",
        },
        "provenance": {"record_doi": EXPECTED[measurement]["doi"],
                       "source_manifest": SOURCE_MANIFEST_PATHS[measurement],
                       "tables": [1,2,3], "retrieved": _manifest(measurement)["retrieved"]},
    }


def _paper_tex() -> tuple[str, bytes]:
    manifest = _manifest("COMPASS_2010_I862410")
    path = ROOT / manifest["source_archive"]["path"]
    payload = path.read_bytes()
    if _sha256(payload) != manifest["source_archive"]["sha256"]:
        raise SIDISTrancheDataError("COMPASS 2010 source archive checksum mismatch")
    member = str(manifest["tex_member"])
    with tarfile.open(fileobj=io.BytesIO(payload), mode="r:gz") as archive:
        stream = archive.extractfile(member)
        if stream is None:
            raise SIDISTrancheDataError(f"Missing {member} in COMPASS 2010 source")
        tex = stream.read()
    return tex.decode("latin-1"), tex


def paper_2010_audit() -> dict[str, Any]:
    tex, tex_payload = _paper_tex()
    # Audit every numeric row by requiring the exact TeX tokens after
    # whitespace normalisation.  The structured transcript below is then the
    # sole source used by the generator and postprocessor.
    compact = re.sub(r"\s+", "", tex)
    missing: list[str] = []
    for row in PAPER_2010_ROWS:
        token = f"{row[0]:.4f}" if row[0] < .1 else f"{row[0]:.3f}"
        if token.rstrip("0") not in compact:
            missing.append(f"asymmetry-x={row[0]}")
    for low, high, correction in zip(X_EDGES_2010[:-1], X_EDGES_2010[1:], PAPER_2010_DEUTERON_CORRECTIONS):
        # Range formatting varies only by trailing zeroes; ensure every row's
        # five numerical values occurs in the table source in sequence.
        values = "&".join("0.~~~" if value == 0 else f"{value:.3f}" for value in correction)
        if values.replace("0.000", "0.~~~") not in compact:
            # Preserve the finding in the audit instead of silently accepting
            # an extraction drift.  Current pinned source yields no entries.
            missing.append(f"deuteron-correction={low}-{high}")
    if missing:
        raise SIDISTrancheDataError(f"COMPASS 2010 TeX extraction audit failed: {missing}")
    manifest = _manifest("COMPASS_2010_I862410")
    return {
        "schema_version": 1, "measurement": "COMPASS_2010_I862410",
        "numerical_authority": "arXiv:1007.4061 paper TeX",
        "source_archive_sha256": manifest["source_archive"]["sha256"],
        "tex_member": manifest["tex_member"], "tex_sha256": _sha256(tex_payload),
        "asymmetry_rows": PAPER_2010_ROWS,
        "correlation_rows": PAPER_2010_CORRELATIONS,
        "deuteron_correction_rows": [
            {"x_low": low, "x_high": high,
             "piplus": row[0], "piminus": row[1], "kplus": row[2],
             "kminus": row[3], "inclusive": row[4]}
            for low, high, row in zip(X_EDGES_2010[:-1], X_EDGES_2010[1:],
                                      PAPER_2010_DEUTERON_CORRECTIONS)
        ],
        "systematic_statements": {
            "multiplicative_fraction": .06,
            "false_asymmetry_upper_bound_stat_fraction": .56,
        },
        "missing_rows": [], "discrepancies": [],
    }


def normalize_compass_2010() -> dict[str, Any]:
    audit = paper_2010_audit()
    datasets: dict[str, Any] = {}
    for species_index, species in enumerate(SPECIES):
        points: list[dict[str, Any]] = []
        offset = 2 + 3 * species_index
        for index, row in enumerate(PAPER_2010_ROWS):
            value, stat, systematic = map(float, row[offset:offset + 3])
            correlated = .06 * value
            points.append({
                "bin": index + 1, "flat_bin": index,
                "x_low": X_EDGES_2010[index], "x_high": X_EDGES_2010[index + 1],
                "x_mean": float(row[0]), "q2_mean": float(row[1]), "a1": value,
                "stat": stat, "systematic": systematic,
                "systematic_correlated_6pct": correlated,
                "systematic_uncorrelated": math.sqrt(max(0., systematic**2 - correlated**2)),
            })
        datasets[species] = {
            "id": species, "pid": SPECIES_PID[species], "observable": "A1p",
            "density_widths": [], "rivet_path": f"/COMPASS_2010_I862410/A1_{species}",
            "raw_objects": {
                "ordinary": f"Yield_{species}_x",
                "inverse_depolarization": f"YieldOverD_{species}_x",
                "covariance": f"CovarianceProxy_{species}_x",
            }, "points": points,
        }
    pairs = [(1,0),(2,0),(2,1),(3,0),(3,1),(3,2)]
    blocks: list[dict[str, Any]] = []
    full = [[0. for _ in range(48)] for _ in range(48)]
    order: list[str] = []
    for xbin in range(12):
        corr = [[1. if i == j else 0. for j in range(4)] for i in range(4)]
        for values, (i,j) in zip(PAPER_2010_CORRELATIONS, pairs):
            corr[i][j] = corr[j][i] = values[xbin]
        stat = [datasets[s]["points"][xbin]["stat"] for s in SPECIES]
        cov = [[corr[i][j] * stat[i] * stat[j] for j in range(4)] for i in range(4)]
        for i, species in enumerate(SPECIES):
            order.append(f"{species}:x{xbin + 1:02d}")
            for j in range(4):
                full[4*xbin+i][4*xbin+j] = cov[i][j]
        blocks.append({"x_bin": xbin + 1, "order": list(SPECIES),
                       "correlation": corr, "covariance": cov})
    return {
        "schema_version": 2, "measurement": "COMPASS_2010_I862410",
        "observable": "identified-hadron A1p(x)",
        "beam": {"pid": -13, "energy_gev": 160.0},
        "binning": {"x": X_EDGES_2010},
        "selection": {"q2_min_gev2": 1., "x": [.004,.7], "y": [.1,.9],
                      "z": [.2,.85], "hadron_momentum_gev": [10.,50.]},
        "datasets": datasets, "point_order": order,
        "statistical_correlation_blocks": blocks, "statistical_covariance": full,
        "target_outputs": {"H": {"P": 1.0}}, "longitudinal_target_scale": 1.0,
        "systematics": {"correlated_multiplicative_fraction": .06,
                        "remaining_variance": "diagonal max(syst^2-(0.06*A1)^2,0)"},
        "depolarization": {"factor": "finite-muon-mass COMPASS D with R1998",
                           "a2_assumption": "eta*A2=0"},
        "published_deuteron_corrections": audit["deuteron_correction_rows"],
        "paper_extraction_audit": AUDIT_PATHS["COMPASS_2010_I862410"],
        "provenance": {"source_manifest": SOURCE_MANIFEST_PATHS["COMPASS_2010_I862410"],
                       "numerical_authority": "paper TeX", "arxiv": "1007.4061"},
    }


def _qualifier_bin(table: Mapping[str, Any], key: str) -> tuple[float, float, float]:
    text = str(table["qualifiers"][key][0]["value"])
    match = re.search(
        r"^\s*([-+0-9.eE]+)\s*\(BIN=\s*([-+0-9.eE]+)\s+TO\s+([-+0-9.eE]+)",
        text,
    )
    if not match:
        raise SIDISTrancheDataError(f"Could not parse qualifier {key}: {text}")
    return tuple(map(float, match.groups()))  # type: ignore[return-value]


def normalize_compass_pt2() -> dict[str, Any]:
    measurement = "COMPASS_2018_I1624692"
    tables = _tables(measurement)
    datasets: dict[str, Any] = {}
    for species, table_range in (("hplus", range(1,82)), ("hminus", range(82,163))):
        points: list[dict[str, Any]] = []
        for table_number in table_range:
            table = tables[table_number]
            xmean, xlow, xhigh = _qualifier_bin(table, "$x$")
            qmean, qlow, qhigh = _qualifier_bin(table, "$Q^2$")
            zmean, zlow, zhigh = _qualifier_bin(table, "$z$")
            ymean = float(str(table["qualifiers"]["$y$"][0]["value"]).split()[0])
            for row in table["values"]:
                pcell = row["x"][0]
                primary = _dependent(row, 0)
                hadron_correction = float(_dependent(row, 1)["value"])
                dis_correction = float(_dependent(row, 2)["value"])
                value = float(primary["value"])
                points.append({
                    "cell": len(points) + 1, "flat_bin": len(points), "table": table_number,
                    "x_low": xlow, "x_high": xhigh, "x_mean": xmean,
                    "q2_low": qlow, "q2_high": qhigh, "q2_mean": qmean,
                    "y_mean": ymean, "z_low": zlow, "z_high": zhigh, "z_mean": zmean,
                    "pt2_low": float(pcell["low"]), "pt2_high": float(pcell["high"]),
                    "pt2_mean": float(pcell.get("value", .5*(float(pcell["low"])+float(pcell["high"])))),
                    "value": value, "stat": _error(primary, "stat"),
                    "systematic": _error(primary, "sys"),
                    "corrections": {"vm_hadron": hadron_correction, "vm_dis": dis_correction},
                    "unsubtracted_value": (
                        value * dis_correction / hadron_correction
                        if hadron_correction > 0.0 else None
                    ),
                })
        if len(points) != 2332:
            raise SIDISTrancheDataError(f"Expected 2332 {species} cells, found {len(points)}")
        slices = _slices(measurement, species, points, ("x","q2","z"), "pt2")
        datasets[species] = {
            "id": species, "pid_mode": "stable_charged_hadron",
            "charge": 1 if species == "hplus" else -1, "z_mass_assumption": "pion",
            "observable": "d2M/(dz dpt2)", "density_widths": ["z","pt2"],
            "flat_rivet_path": f"/{measurement}/Multiplicity_{species}_cells",
            "raw_objects": {
                "numerator": f"HadronNumerator_{species}_cells",
                "denominator": f"DISDenominator_{species}_cells",
                "covariance": f"CovarianceProxy_{species}_cells",
            },
            "valid_cell_map": [_cell_signature(p, ("x","q2","z","pt2")) for p in points],
            "points": points, "slices": slices,
        }
    if sum(len(d["points"]) for d in datasets.values()) != EXPECTED[measurement]["cells"]:
        raise SIDISTrancheDataError("COMPASS 2018 complete cell inventory changed")
    return {
        "schema_version": 2, "measurement": measurement,
        "observable": "isoscalar d2M(h+/-)/(dz dP_hT^2)",
        "beam": {"pid": -13, "energy_gev": 160.},
        "selection": {"q2_min_gev2": 1., "w_min_gev": 5., "x": [.003,.4],
                      "y": [.1,.9], "z": [.2,.8], "pt2_gev2": [.02,3.]},
        "datasets": datasets, "target_outputs": {"D": {"P": .5, "N": .5}},
        "systematics": {"statistical_covariance": "diagonal; not released",
                        "systematic_covariance": "diagonal; not released"},
        "corrections": {"primary": "vector-meson-subtracted", "unsubtracted_retained": True,
                        "herwig": "no detector or experimental correction applied"},
        "source_count_audit": {
            "assessment_claim": 4918, "hepdata_v1_rows": 4664,
            "disposition": "use all released rows; never invent the 254-point difference",
        },
        "provenance": {"record_doi": EXPECTED[measurement]["doi"],
                       "source_manifest": SOURCE_MANIFEST_PATHS[measurement],
                       "tables": list(range(1,163)), "retrieved": _manifest(measurement)["retrieved"]},
    }


def _archive_files() -> tuple[tarfile.TarFile, dict[str, tarfile.TarInfo]]:
    measurement = "HERMES_2013_I1208547"
    manifest = _manifest(measurement)
    entry = manifest["full_archive"]
    path = ROOT / str(entry["path"])
    payload = path.read_bytes()
    expected = EXPECTED[measurement]
    if len(payload) != expected["archive_bytes"] or _sha256(payload) != expected["archive_sha256"]:
        raise SIDISTrancheDataError("HERMES full archive size/checksum mismatch")
    archive = tarfile.open(fileobj=io.BytesIO(payload), mode="r:gz")
    members = {member.name: member for member in archive.getmembers()
               if member.isfile() and "/._" not in member.name and not Path(member.name).name.startswith("._")}
    return archive, members


def _archive_text(archive: tarfile.TarFile, members: Mapping[str, tarfile.TarInfo],
                  suffix: str) -> str:
    matches = [name for name in members if name.endswith(suffix)]
    if len(matches) != 1:
        raise SIDISTrancheDataError(f"Expected one HERMES archive member ending {suffix}; got {matches}")
    stream = archive.extractfile(members[matches[0]])
    if stream is None:
        raise SIDISTrancheDataError(f"Could not extract HERMES member {matches[0]}")
    payload = stream.read()
    if suffix.endswith(".gz"):
        payload = gzip.decompress(payload)
    return payload.decode("utf-8")


def _numeric_rows(text: str) -> list[list[float]]:
    rows: list[list[float]] = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        try:
            rows.append([float(token) for token in stripped.split()])
        except ValueError as exc:
            raise SIDISTrancheDataError(f"Invalid HERMES numeric row: {line}") from exc
    return rows


def _hermes_layout(config: Mapping[str, Any]) -> list[dict[str, float]]:
    axis = str(config["axis"])
    axis_edges = list(config["axis_edges"])
    z_edges = list(config["z_edges"])
    ph_edges = list(config["phperp_edges"])
    result: list[dict[str, float]] = []
    for iaxis in range(len(axis_edges)-1):
        for iz in range(len(z_edges)-1):
            for ip in range(len(ph_edges)-1):
                cell = {
                    f"{axis}_low": axis_edges[iaxis], f"{axis}_high": axis_edges[iaxis+1],
                    "z_low": z_edges[iz], "z_high": z_edges[iz+1],
                    "phperp_low": ph_edges[ip], "phperp_high": ph_edges[ip+1],
                }
                if axis == "q2":
                    cell["x_low"], cell["x_high"] = config["x_edges"]
                result.append(cell)
    return result


def _matrix_sha256(matrix: Sequence[Sequence[float]]) -> str:
    return _sha256(json.dumps(matrix, separators=(",", ":")).encode("utf-8"))


def _normalize_released_covariance(
    matrix: Sequence[Sequence[float]], label: str
) -> tuple[list[list[float]], dict[str, Any]]:
    """Make a released covariance mathematically usable with a frozen audit.

    Several HERMES archive matrices contain transpose-level numerical
    differences and nine one-sided decimal/exponent corruptions.  A covariance
    is symmetric and obeys |Cij| <= sqrt(Cii Cjj).  Ordinary asymmetric pairs
    are averaged.  If exactly one direction violates that necessary bound, the
    valid transpose partner is retained.  No diagonal is changed.
    """
    size = len(matrix)
    if any(len(row) != size for row in matrix):
        raise SIDISTrancheDataError(f"HERMES {label} covariance is not square")
    result = [list(map(float, row)) for row in matrix]
    raw_hash = _matrix_sha256(result)
    asymmetric = 0
    max_asymmetry = 0.0
    replacements: list[dict[str, Any]] = []
    for i in range(size):
        if not math.isfinite(result[i][i]) or result[i][i] < 0.0:
            raise SIDISTrancheDataError(
                f"HERMES {label} covariance has invalid diagonal {i}"
            )
        for j in range(i+1, size):
            forward, reverse = result[i][j], result[j][i]
            difference = abs(forward-reverse)
            if difference:
                asymmetric += 1
                max_asymmetry = max(max_asymmetry, difference)
            bound = math.sqrt(result[i][i]*result[j][j])
            tolerance = bound*1.e-10+1.e-15
            forward_valid = abs(forward) <= bound+tolerance
            reverse_valid = abs(reverse) <= bound+tolerance
            if forward_valid and not reverse_valid:
                chosen, disposition = forward, "used forward; reverse violates bound"
            elif reverse_valid and not forward_valid:
                chosen, disposition = reverse, "used reverse; forward violates bound"
            elif forward_valid and reverse_valid:
                chosen, disposition = .5*(forward+reverse), "averaged transpose pair"
            else:
                raise SIDISTrancheDataError(
                    f"HERMES {label} covariance pair ({i},{j}) violates "
                    "Cauchy-Schwarz in both directions"
                )
            if forward_valid != reverse_valid:
                replacements.append({
                    "i": i, "j": j, "forward": forward, "reverse": reverse,
                    "cauchy_schwarz_bound": bound, "chosen": chosen,
                    "disposition": disposition,
                })
            result[i][j] = result[j][i] = chosen
    return result, {
        "dataset": label,
        "policy": (
            "average ordinary transpose differences; when exactly one direction "
            "violates |Cij|<=sqrt(Cii*Cjj), retain the valid transpose partner"
        ),
        "diagonal_changed": False,
        "raw_sha256": raw_hash,
        "normalized_sha256": _matrix_sha256(result),
        "raw_asymmetric_pair_count": asymmetric,
        "maximum_absolute_transpose_difference": max_asymmetry,
        "cauchy_schwarz_replacement_count": len(replacements),
        "cauchy_schwarz_replacements": replacements,
    }


def hermes_projection_audit(tables: Mapping[int, Mapping[str, Any]],
                            archive: tarfile.TarFile,
                            members: Mapping[str, tarfile.TarInfo]) -> dict[str, Any]:
    """Freeze all 64 official projections and report archive differences.

    The archive carries more significant digits than HEPData.  Values are
    compared after rounding to the decimal precision published by HEPData;
    every non-matching field is retained in ``discrepancies``.
    """
    inventory: list[dict[str, Any]] = []
    for number, table in sorted(tables.items()):
        inventory.append({
            "table": number, "doi": table["doi"], "name": table["name"],
            "description": table.get("description", ""),
            "headers": [str(item["name"]) for item in table["headers"]],
            "rows": len(table["values"]),
            "values": table["values"],
        })
    # The checksum-pinned HEPData serialisation is the official cross-check
    # surface.  Archive projection filenames are inventoried so missing files
    # cannot be hidden; every coordinate, value, statistical error and
    # systematic error is compared at the precision printed by HEPData.
    projection_members = sorted(
        name for name in members
        if any(part in name for part in ("/z-proj/", "/pt-proj/", "/x-proj/", "/Q2-proj/"))
        and name.endswith(".list.gz")
    )
    if not projection_members:
        raise SIDISTrancheDataError("HERMES projection inventory is empty")

    def source_rows(suffix: str) -> list[list[float]]:
        return _numeric_rows(_archive_text(archive, members, suffix))

    def source_specs(number: int) -> list[dict[str, Any]]:
        """Map each HEPData dependent-variable group to an archive projection."""
        if number <= 4:
            species = SPECIES[number-1]
            return [{
                "label": target, "group": group, "coordinate_column": 6,
                "indices": list(range(10)),
                "suffix": (
                    f"/z-3D/z-proj/hermes.{target}.z-3D.z-proj."
                    f"vmsub.mults_{species}.list.gz"
                ),
            } for group, target in enumerate(("proton", "deuteron"))]
        if number <= 8:
            species = SPECIES[number-5]
            return [{
                "label": vm, "group": group, "coordinate_column": 6,
                "indices": list(range(10)),
                "suffix": (
                    f"/z-3D/z-proj/hermes.proton.z-3D.z-proj."
                    f"{vm}.mults_{species}.list.gz"
                ),
            } for group, vm in enumerate(("vmsub", "no-vmsub"))]
        if number <= 12:
            species = SPECIES[number-9]
            return [{
                "label": "target-asymmetry", "group": 0,
                "coordinate_column": 6, "indices": list(range(10)),
                "suffix": (
                    f"/z-3D/z-proj/hermes.z-3D.z-proj."
                    f"vmsub.asymm_{species}.list.gz"
                ),
            }]
        if number <= 16:
            z_index = number-12
            return [{
                "label": species, "group": group, "coordinate_column": 5,
                "indices": [axis*6+z_index for axis in range(9)],
                "suffix": (
                    f"/zx-3D/zx-proj/hermes.zx-3D.zx-proj."
                    f"vmsub.asymm_{species}.list.gz"
                ),
            } for group, species in enumerate(("piplus", "piminus"))]

        if number <= 32:
            offset, config, projection, coordinate = number-17, "zpt-3D", "zpt-proj", 7
        elif number <= 48:
            offset, config, projection, coordinate = number-33, "zx-3D", "zx-proj", 5
        else:
            offset, config, projection, coordinate = number-49, "zQ2-3D", "zQ2-proj", 4
        species = SPECIES[offset//4]
        z_index = offset % 4 + 1
        if config == "zpt-3D":
            indices = list(range(z_index*9, (z_index+1)*9))
        else:
            indices = [axis*6+z_index for axis in range(9)]
        return [{
            "label": target, "group": group, "coordinate_column": coordinate,
            "indices": indices,
            "suffix": (
                f"/{config}/{projection}/hermes.{target}.{config}."
                f"{projection}.vmsub.mults_{species}.list.gz"
            ),
        } for group, target in enumerate(("proton", "deuteron"))]

    def printed_tolerance(token: Any) -> float:
        decimal = Decimal(str(token))
        return float(Decimal("0.5") * (Decimal(10) ** decimal.as_tuple().exponent))

    discrepancies: list[dict[str, Any]] = []
    comparisons = 0
    for number, table in sorted(tables.items()):
        for spec in source_specs(number):
            rows = source_rows(str(spec["suffix"]))
            indices = list(spec["indices"])
            if len(indices) != len(table["values"]):
                raise SIDISTrancheDataError(
                    f"HERMES table {number} projection length changed"
                )
            for point_index, (published, archive_index) in enumerate(
                zip(table["values"], indices), start=1
            ):
                row = rows[int(archive_index)]
                dependent = published["y"][int(spec["group"])]
                errors = {item["label"]: item["symerror"]
                          for item in dependent["errors"]}
                fields = {
                    "coordinate": (
                        row[int(spec["coordinate_column"])],
                        published["x"][int(spec["group"])]["value"],
                    ),
                    "value": (row[1], dependent["value"]),
                    "stat": (row[2], errors["stat"]),
                    "sys": (row[3], errors["sys"]),
                }
                for field, (archive_value, published_token) in fields.items():
                    comparisons += 1
                    tolerance = printed_tolerance(published_token)
                    if abs(float(archive_value)-float(published_token)) > tolerance+1.e-12:
                        discrepancies.append({
                            "table": number, "group": spec["label"],
                            "point": point_index, "field": field,
                            "archive": archive_value,
                            "hepdata": float(published_token),
                            "hepdata_printed_tolerance": tolerance,
                            "archive_member_suffix": spec["suffix"],
                        })
    return {
        "schema_version": 1, "measurement": "HERMES_2013_I1208547",
        "authority": "VM-subtracted full five-binning archive",
        "crosscheck": "complete 64-table HEPData v1 submission",
        "comparison_policy": (
            "HEPData projection values are preserved verbatim; differences from the "
            "higher-precision archive projections are never substituted into 3D data."
        ),
        "hepdata_tables": inventory,
        "archive_projection_member_count": len(projection_members),
        "archive_projection_members": projection_members,
        "numeric_comparison_count": comparisons,
        "discrepancy_count": len(discrepancies),
        "discrepancies": discrepancies,
        "note": (
            "No discrepancies were found at the precision published in HEPData; "
            "unrounded archive values remain the 3D authority."
            if not discrepancies else
            "All archive/HEPData projection discrepancies are listed explicitly; "
            "the unrounded archive remains the 3D authority."
        ),
    }


def _hermes_readable_projections(
    tables: Mapping[int, Mapping[str, Any]],
) -> dict[str, Any]:
    """Return official projections tied to separately integrated raw yields."""

    measurement = "HERMES_2013_I1208547"
    projections: dict[str, Any] = {}

    def points(
        table_number: int,
        group: int,
        edges: Sequence[float],
        axis: str,
        *,
        first_row: int = 0,
        z_interval: tuple[float, float] | None = None,
    ) -> list[dict[str, Any]]:
        rows = list(tables[table_number]["values"])[first_row:]
        if len(rows)+1 != len(edges):
            raise SIDISTrancheDataError(
                f"HERMES projection table {table_number} edge count changed"
            )
        result: list[dict[str, Any]] = []
        for index, row in enumerate(rows):
            dependent = row["y"][group]
            errors = {
                item["label"]: float(item["symerror"])
                for item in dependent["errors"]
            }
            point = {
                "flat_bin": index,
                f"{axis}_low": float(edges[index]),
                f"{axis}_high": float(edges[index+1]),
                f"{axis}_mean": float(row["x"][group]["value"]),
                "value": float(dependent["value"]),
                "stat": errors["stat"],
                "systematic": errors["sys"],
                "hepdata_table": table_number,
                "hepdata_row": index + first_row + 1,
            }
            if z_interval is not None:
                point["z_low"], point["z_high"] = z_interval
            result.append(point)
        return result

    def add(
        target: str,
        species: str,
        table_number: int,
        group: int,
        configuration: str,
        key: str,
        edges: Sequence[float],
        axis: str,
        density_widths: Sequence[str],
        *,
        first_row: int = 0,
        z_interval: tuple[float, float] | None = None,
    ) -> None:
        identifier = f"{target}_{species}_{key}"
        projections[identifier] = {
            "id": identifier,
            "integrated_projection": True,
            "published_target": target,
            "species": species,
            "binning": configuration,
            "projection_key": key,
            "axis": axis,
            "density_widths": list(density_widths),
            "edges": list(map(float, edges)),
            "rivet_path": f"/{measurement}/Projection_{identifier}",
            "raw_objects": {
                "numerator": (
                    f"ProjectionHadronNumerator_{species}_{configuration}_{key}"
                ),
                "denominator": (
                    f"ProjectionDISDenominator_{configuration}_{key}"
                ),
                "covariance": (
                    f"ProjectionCovarianceProxy_{species}_{configuration}_{key}"
                ),
            },
            "source": "official HEPData v1 projection cross-check",
            "points": points(
                table_number, group, edges, axis,
                first_row=first_row, z_interval=z_interval,
            ),
        }

    z_edges = [.2,.25,.3,.4,.5,.6,.7,.8,1.1]
    for species_index, species in enumerate(SPECIES):
        for group, target in enumerate(("H", "D")):
            add(target, species, species_index+1, group, "z-3D", "z",
                z_edges, "z", ["z"], first_row=2)

    fixed_z = ((.2,.3),(.3,.4),(.4,.6),(.6,.8))
    families = (
        (17, "zpt-3D", "pt", "phperp",
         [0.,.1,.2,.3,.4,.5,.6,.7,.8,1.2], ["z","phperp"]),
        (33, "zx-3D", "x", "x",
         [.023,.04,.055,.075,.1,.14,.2,.3,.4,.6], ["z"]),
        (49, "zQ2-3D", "q2", "q2",
         [1.,1.25,1.5,1.75,2.,2.25,2.5,3.,5.,15.], ["z"]),
    )
    for first_table, configuration, prefix, axis, edges, widths in families:
        for species_index, species in enumerate(SPECIES):
            for z_index, interval in enumerate(fixed_z):
                table_number = first_table + 4*species_index + z_index
                token = (
                    f"z{int(round(100*interval[0])):03d}_"
                    f"{int(round(100*interval[1])):03d}"
                )
                key = f"{prefix}_{token}"
                for group, target in enumerate(("H", "D")):
                    add(target, species, table_number, group, configuration,
                        key, edges, axis, widths, z_interval=interval)
    if len(projections) != 104:
        raise SIDISTrancheDataError(
            f"HERMES readable projection inventory changed: {len(projections)}"
        )
    return projections


def normalize_hermes() -> tuple[dict[str, Any], dict[str, Any]]:
    measurement = "HERMES_2013_I1208547"
    tables = _tables(measurement)
    archive, members = _archive_files()
    datasets: dict[str, Any] = {}
    covariance_audits: dict[str, Any] = {}
    try:
        for config_name, config in HERMES_CONFIGS.items():
            layout = _hermes_layout(config)
            usable_indices = [i for i, cell in enumerate(layout) if cell["z_low"] >= .2]
            if len(usable_indices) != int(config["usable"]):
                raise SIDISTrancheDataError(f"HERMES {config_name} valid-cell map changed")
            for target, target_name in (("H","proton"),("D","deuteron")):
                covariance_rows = _numeric_rows(_archive_text(
                    archive, members,
                    f"/full/hermes.{target_name}.{config_name}.vmsub.covmat_mults.list.gz",
                ))
                nfull = len(layout)
                if len(covariance_rows) != nfull*nfull:
                    raise SIDISTrancheDataError(f"HERMES {target}/{config_name} covariance dimensions changed")
                covariance_pairs = {
                    (int(row[0]), int(row[1])) for row in covariance_rows
                    if len(row) >= 2
                    and row[0] == int(row[0]) and row[1] == int(row[1])
                }
                if covariance_pairs != {
                    (i, j) for i in range(nfull) for j in range(nfull)
                } or any(len(row) != 2+len(SPECIES) for row in covariance_rows):
                    raise SIDISTrancheDataError(
                        f"HERMES {target}/{config_name} covariance grid changed"
                    )
                covariance_columns = {species: [[0.]*len(usable_indices) for _ in usable_indices]
                                      for species in SPECIES}
                index_map = {full: reduced for reduced, full in enumerate(usable_indices)}
                for row in covariance_rows:
                    i, j = int(row[0]), int(row[1])
                    if i not in index_map or j not in index_map:
                        continue
                    for column, species in enumerate(SPECIES, start=2):
                        covariance_columns[species][index_map[i]][index_map[j]] = row[column]
                for species in SPECIES:
                    rows = _numeric_rows(_archive_text(
                        archive, members,
                        f"/full/hermes.{target_name}.{config_name}.vmsub.mults_{species}.list.gz",
                    ))
                    if len(rows) != nfull or any(int(row[0]) != i for i,row in enumerate(rows)):
                        raise SIDISTrancheDataError(f"HERMES {target}/{config_name}/{species} row order changed")
                    points: list[dict[str, Any]] = []
                    for flat_bin, full_index in enumerate(usable_indices):
                        row = rows[full_index]
                        point = {
                            "cell": flat_bin + 1, "flat_bin": flat_bin,
                            "archive_index": full_index, **layout[full_index],
                            "value": row[1], "stat": row[2], "systematic": row[3],
                            "q2_mean": row[4], "x_mean": row[5], "z_mean": row[6],
                            "phperp_mean": row[7],
                        }
                        points.append(point)
                    covariance = covariance_columns[species]
                    for i, point in enumerate(points):
                        if not math.isclose(covariance[i][i], point["stat"]**2,
                                            rel_tol=2.e-5, abs_tol=2.e-12):
                            raise SIDISTrancheDataError(
                                f"HERMES covariance/stat mismatch {target}/{config_name}/{species}/{i}"
                            )
                    dataset_id = f"{target}_{species}_{config_name.replace('-', '_')}"
                    covariance, covariance_audits[dataset_id] = (
                        _normalize_released_covariance(covariance, dataset_id)
                    )
                    slices = _slices(measurement, dataset_id, points,
                                     (str(config["axis"]), "z"), "phperp")
                    datasets[dataset_id] = {
                        "id": dataset_id, "published_target": target,
                        "species": species, "pid": SPECIES_PID[species],
                        "binning": config_name, "observable": "d2M/(dz dphperp)",
                        "density_widths": ["z","phperp"],
                        "flat_rivet_path": f"/{measurement}/Multiplicity_{dataset_id}_cells",
                        "raw_objects": {
                            "numerator": f"HadronNumerator_{species}_{config_name}_cells",
                            "denominator": f"DISDenominator_{config_name}_cells",
                            "covariance": f"CovarianceProxy_{species}_{config_name}_cells",
                        },
                        "valid_cell_map": [_cell_signature(p, (str(config["axis"]),"z","phperp")) for p in points],
                        "points": points, "statistical_covariance": covariance,
                        "slices": slices,
                    }
        if sum(len(dataset["points"]) for dataset in datasets.values()) != 6592:
            raise SIDISTrancheDataError("HERMES usable target/species cell count changed")
        audit = hermes_projection_audit(tables, archive, members)
        readable_projections = _hermes_readable_projections(tables)
    finally:
        archive.close()
    return ({
        "schema_version": 2, "measurement": measurement,
        "observable": "H/D identified-hadron multiplicities in five independent 3D binnings",
        "beam": {"pid": "e+/-", "energy_gev": 27.6},
        "selection": {"q2_min_gev2": 1., "w2_min_gev2": 10., "y": [.1,.85],
                      "hadron_momentum_gev": [2.,15.], "z_lower_edge_min": .2,
                      "high_z_bins": "retained and annotated"},
        "binnings": HERMES_CONFIGS, "datasets": datasets,
        "readable_projections": readable_projections,
        "target_outputs": {"H": {"P": 1.0}, "D": {"P": .5, "N": .5}},
        "statistical_covariance": "dense within each target/species/binning",
        "covariance_normalization": {
            "reason": (
                "released matrices contain transpose differences and isolated "
                "one-sided Cauchy-Schwarz violations"
            ),
            "datasets": covariance_audits,
        },
        "systematics": "point-to-point diagonal",
        "fit_policy": "never combine overlapping binnings into one goodness of fit",
        "projection_policy": "integrate numerator and DIS denominator separately before ratio",
        "projection_audit": AUDIT_PATHS[measurement],
        "provenance": {"record_doi": EXPECTED[measurement]["doi"],
                       "source_manifest": SOURCE_MANIFEST_PATHS[measurement],
                       "archive_sha256": EXPECTED[measurement]["archive_sha256"],
                       "hepdata_tables": list(range(1,65))},
    }, audit)


def normalized_from_raw(measurement: str) -> dict[str, Any]:
    if measurement in {"COMPASS_2025_I2840545", "COMPASS_2026_I3096394"}:
        return normalize_modern_compass(measurement)
    if measurement == "COMPASS_2010_I862410":
        return normalize_compass_2010()
    if measurement == "COMPASS_2018_I1624692":
        return normalize_compass_pt2()
    if measurement == "HERMES_2013_I1208547":
        return normalize_hermes()[0]
    raise SIDISTrancheDataError(f"Unknown SIDIS tranche measurement {measurement}")


def _assert_equal(generated: Any, snapshot: Any, path: str = "$") -> None:
    if type(generated) is not type(snapshot):
        raise SIDISTrancheDataError(f"Snapshot type mismatch at {path}")
    if isinstance(generated, Mapping):
        if set(generated) != set(snapshot):
            raise SIDISTrancheDataError(f"Snapshot key mismatch at {path}")
        for key in generated:
            _assert_equal(generated[key], snapshot[key], f"{path}.{key}")
    elif isinstance(generated, list):
        if len(generated) != len(snapshot):
            raise SIDISTrancheDataError(f"Snapshot length mismatch at {path}")
        for index, (left, right) in enumerate(zip(generated, snapshot)):
            _assert_equal(left, right, f"{path}[{index}]")
    elif isinstance(generated, float):
        if not math.isclose(generated, snapshot, rel_tol=1.e-14, abs_tol=1.e-24):
            raise SIDISTrancheDataError(f"Snapshot float mismatch at {path}: {generated} != {snapshot}")
    elif generated != snapshot:
        raise SIDISTrancheDataError(f"Snapshot mismatch at {path}: {generated!r} != {snapshot!r}")


def _validate_psd(matrix: Sequence[Sequence[float]], label: str) -> None:
    try:
        import numpy as np
    except ImportError as exc:
        raise SIDISTrancheDataError("NumPy is required for covariance validation") from exc
    array = np.asarray(matrix, dtype=float)
    if array.shape[0] != array.shape[1] or not np.allclose(array, array.T, rtol=0, atol=1.e-10):
        raise SIDISTrancheDataError(f"{label} is not symmetric")
    eigenvalues = np.linalg.eigvalsh(array)
    tolerance = max(1., float(np.max(np.abs(eigenvalues)))) * 2.e-7
    if float(eigenvalues[0]) < -tolerance:
        raise SIDISTrancheDataError(f"{label} is not positive semidefinite: {eigenvalues[0]}")


def validate_vendored(measurement: str, *, full_covariance: bool = False) -> dict[str, Any]:
    generated = normalized_from_raw(measurement)
    snapshot = _load_json(ROOT / REFERENCE_PATHS[measurement])
    _assert_equal(generated, snapshot)
    if measurement == "COMPASS_2010_I862410":
        _assert_equal(paper_2010_audit(), _load_json(ROOT / AUDIT_PATHS[measurement]))
        _validate_psd(snapshot["statistical_covariance"], "COMPASS 2010 statistical covariance")
    elif measurement == "HERMES_2013_I1208547":
        _, audit = normalize_hermes()
        _assert_equal(audit, _load_json(ROOT / AUDIT_PATHS[measurement]))
        if full_covariance:
            for dataset in snapshot["datasets"].values():
                _validate_psd(dataset["statistical_covariance"], dataset["id"])
    return snapshot


def _import_yoda() -> Any:
    try:
        import yoda  # type: ignore
    except (ImportError, OSError) as exc:
        raise SIDISTrancheDataError("YODA Python bindings are unavailable; load herwig/pol") from exc
    return yoda


def _estimate(yoda: Any, edges: Sequence[float], points: Sequence[Mapping[str, Any]],
              path: str, value_key: str) -> Any:
    result = yoda.BinnedEstimate1D(list(map(float, edges)), "/REF" + path)
    result.setAnnotation("IsRef", 1)
    for index, point in enumerate(points, start=1):
        target = result.bin(index)
        target.setVal(float(point[value_key]))
        target.setErr(-float(point["stat"]), float(point["stat"]), "stat")
        systematic = float(point["systematic"])
        target.setErr(-systematic, systematic, "syst")
    return result


def _write_yoda(objects: Iterable[Any], destination: Path) -> Path:
    yoda = _import_yoda()
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="sidis-ref-", dir=destination.parent) as temporary:
        plain = Path(temporary) / destination.with_suffix("").name
        yoda.write(list(objects), str(plain))
        compressed = Path(temporary) / destination.name
        with plain.open("rb") as source, compressed.open("wb") as target:
            with gzip.GzipFile(filename="", mode="wb", fileobj=target, mtime=0) as stream:
                shutil.copyfileobj(source, stream)
        os.replace(compressed, destination)
    return destination


def write_reference_yoda(measurement: str, snapshot: Mapping[str, Any] | None = None) -> Path:
    snapshot = snapshot or validate_vendored(measurement)
    yoda = _import_yoda()
    objects: list[Any] = []
    if measurement == "COMPASS_2010_I862410":
        for dataset in snapshot["datasets"].values():
            objects.append(_estimate(yoda, snapshot["binning"]["x"], dataset["points"],
                                     dataset["rivet_path"], "a1"))
    else:
        for dataset in snapshot["datasets"].values():
            flat = _estimate(yoda, range(len(dataset["points"])+1), dataset["points"],
                             dataset["flat_rivet_path"], "value")
            flat.setAnnotation("Layout", "published flattened cells")
            objects.append(flat)
            for item in dataset.get("slices", []):
                points = [dataset["points"][int(index)] for index in item["flat_bins"]]
                objects.append(_estimate(yoda, item["edges"], points, item["rivet_path"], "value"))
        for projection in snapshot.get("readable_projections", {}).values():
            result = _estimate(
                yoda, projection["edges"], projection["points"],
                projection["rivet_path"], "value",
            )
            result.setAnnotation(
                "ProjectionPolicy",
                "integrate generator numerator and DIS denominator separately",
            )
            objects.append(result)
    return _write_yoda(objects, ROOT / REFERENCE_YODA_PATHS[measurement])


def write_normalized_snapshot(measurement: str) -> Path:
    if measurement == "HERMES_2013_I1208547":
        snapshot, audit = normalize_hermes()
        _atomic_json(ROOT / AUDIT_PATHS[measurement], audit)
    else:
        snapshot = normalized_from_raw(measurement)
        if measurement == "COMPASS_2010_I862410":
            _atomic_json(ROOT / AUDIT_PATHS[measurement], paper_2010_audit())
    return _atomic_json(ROOT / REFERENCE_PATHS[measurement], snapshot)


def fetch_and_validate(measurement: str, cache_directory: Path | None = None) -> list[Path]:
    snapshot = validate_vendored(measurement)
    outputs: list[Path] = []
    cache = cache_directory or ROOT / "campaigns/phenomenology/_data_cache" / measurement
    cache.mkdir(parents=True, exist_ok=True)
    if measurement in HEPDATA_MEASUREMENTS:
        for entry in _manifest(measurement)["tables"]:
            request = urllib.request.Request(entry["url"], headers={"User-Agent": "HerwigPol-reference-data/1.0"})
            with urllib.request.urlopen(request, timeout=90) as response:
                payload = response.read()
            if _sha256(payload) != entry["sha256"]:
                raise SIDISTrancheDataError(f"Remote source checksum changed: {entry['url']}")
            destination = cache / f"table-{int(entry['number']):03d}.json"
            destination.write_bytes(payload)
            outputs.append(destination)
    write_reference_yoda(measurement, snapshot)
    return outputs


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("measurement", choices=sorted(MEASUREMENTS))
    parser.add_argument("--write-snapshot", action="store_true")
    parser.add_argument("--write-yoda", action="store_true")
    parser.add_argument("--fetch", action="store_true")
    parser.add_argument("--full-covariance", action="store_true")
    args = parser.parse_args()
    if args.write_snapshot:
        write_normalized_snapshot(args.measurement)
    snapshot = validate_vendored(args.measurement, full_covariance=args.full_covariance)
    if args.fetch:
        fetch_and_validate(args.measurement)
    elif args.write_yoda:
        write_reference_yoda(args.measurement, snapshot)
    print(f"validated {args.measurement}")
