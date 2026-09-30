#!/usr/bin/env python3
"""Reproduce the pinned HERMES Born A_parallel references and statistical covariance.

These are the unfolded Born column, not the measured-asymmetry column.  No
depolarization, beam/target polarization, nuclear or normalization correction
is applied to published values.  Asymmetries and their covariance alone cannot
provide a phase-space-integrated asymmetry without unpolarized yield weights.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DATA = Path("data/experimental/HERMES_2007_I726689")
REFERENCE_PATH = DATA / "born-apar-reference.json"
MEASUREMENT = "HERMES_2007_I726689"
PUBLICATION = "https://arxiv.org/abs/hep-ex/0609039"
THESIS = "https://www.desy.de/~w3hermes/05-LIB/m.ehrenfried.thesis.ps.gz"
THESIS_SHA256 = "3dabaf34773fd856d97e3420a95abf23b9258dee5f79d3ecd6d82cad710b28a7"
FIGURE = "https://www.desy.de/~w3hermes/TRANS/lara.g1_long-4.eps.gz"
FIGURE_SHA256 = "3d490b04ede8572ba7f93dc0d41fa1870b715fe78f72572083950574dc154751"

# Appendix C, p. 97 of the collaboration-listed analysis thesis supplies the
# numerical boundaries. Figure 4 of the final paper retains these internal
# boundaries, removes the first thesis x cell, and imposes 0.18 < Q2 < 20.
PHYSICAL_X_EDGES = [
    0.00406, 0.00733, 0.01177, 0.01677, 0.02124, 0.02948,
    0.03619, 0.04442, 0.05681, 0.07265, 0.09291, 0.11882,
    0.15195, 0.19432, 0.24850, 0.31780, 0.40641, 0.51974,
    0.66466, 0.90000,
]
DISPLAY_X_EDGES = [
    0.0041, 0.0073, 0.0118, 0.0168, 0.0212, 0.0295,
    0.0362, 0.0444, 0.0568, 0.0727, 0.0929, 0.119,
    0.152, 0.194, 0.249, 0.318, 0.406, 0.520, 0.665, 0.9,
]
Q2_EDGES = [[0.18, 1.0] for _ in range(4)] + [
    [0.18, 1.0, 20.0] for _ in range(4)
] + [
    [1.0, 1.505, 2.265, 20.0],
    [1.0, 1.620, 2.623, 20.0],
    [1.0, 1.740, 3.026, 20.0],
    [1.0, 1.882, 3.491, 20.0],
    [1.0, 2.061, 4.032, 20.0],
    [1.0, 2.237, 4.614, 20.0],
    [1.0, 2.657, 5.491, 20.0],
    [1.0, 3.305, 6.645, 20.0],
    [1.0, 4.093, 7.967, 20.0],
    [1.0, 5.043, 9.458, 20.0],
    [1.0, 7.655, 12.528, 20.0],
]

SOURCES = {
    7: (10813, "faa0f3fc562ebd8fc96148a57f43b95367ff0b37a1b747c1f740f6b8613412ae"),
    8: (10818, "8b0d66a33e0aea0c00ae3b5cc0f7ade61ead7cd09f7ebb5cbebc09bb23abfed3"),
    18: (10869, "58ac14270c80bcc44153d7c6af4e739fb447c1c29f42746125c1d4f50aa42b7e"),
    19: (10876, "765085652a070ad0ad096918b7281f9f4d50a5a66374811d708ac91caa4e7c4a"),
}
MISSING_X_MEANS = {4: 0.0190, 10: 0.0403, 12: 0.0506}


class BornReferenceError(RuntimeError):
    """A reference source, bin identity or normalized product failed validation."""


def _source(table: int) -> dict[str, Any]:
    identifier, digest = SOURCES[table]
    return {
        "name": f"HEPData Table {table}",
        "url": f"https://www.hepdata.net/record/data/11211/{identifier}/1",
        "doi": f"10.17182/hepdata.11211.v1/t{table}",
        "sha256": digest,
        "raw": str(DATA / f"raw-born-apar-table-{table}.json"),
    }


def _load_table(table: int, root: Path) -> dict[str, Any]:
    source = _source(table)
    payload = (root / source["raw"]).read_bytes()
    if hashlib.sha256(payload).hexdigest() != source["sha256"]:
        raise BornReferenceError(f"Checksum changed for {source['name']}")
    value = json.loads(payload)
    if value.get("name") != f"Table {table}" or value.get("doi") != source["doi"]:
        raise BornReferenceError(f"Table identity changed for {source['name']}")
    expected = 45 if table in (7, 8) else 2025
    if len(value.get("values", [])) != expected:
        raise BornReferenceError(f"Row count changed for {source['name']}")
    return value


def _uncertainty(cell: dict[str, Any], label: str) -> float:
    matches = [item for item in cell["errors"] if item.get("label") == label]
    if len(matches) != 1 or "symerror" not in matches[0]:
        raise BornReferenceError(f"Missing or ambiguous {label} uncertainty")
    return float(matches[0]["symerror"])


def _correlation(table: dict[str, Any]) -> list[list[float]]:
    matrix = [[None for _ in range(45)] for _ in range(45)]
    for row in table["values"]:
        i, j = [int(float(item["value"])) - 1 for item in row["x"]]
        if not (0 <= i < 45 and 0 <= j < 45) or matrix[i][j] is not None:
            raise BornReferenceError("Invalid or repeated correlation coordinates")
        matrix[i][j] = float(row["y"][0]["value"])
    for i in range(45):
        for j in range(45):
            value = matrix[i][j]
            if value is None or not math.isfinite(value) or abs(value) > 1:
                raise BornReferenceError("Incomplete or invalid correlation matrix")
            if value != matrix[j][i] or (i == j and value != 1):
                raise BornReferenceError("Non-symmetric correlation matrix or non-unit diagonal")
    return matrix


def build_snapshot(root: Path = ROOT) -> dict[str, Any]:
    figure = (root / DATA / "raw-born-apar-fig4.eps.gz").read_bytes()
    if hashlib.sha256(figure).hexdigest() != FIGURE_SHA256:
        raise BornReferenceError("Checksum changed for published Figure 4")
    if b"1340 1032 m 73 X" not in gzip.decompress(figure):
        raise BornReferenceError("Published Figure 4 cell boundary missing")
    datasets, covariance = [], {}
    for target, asymmetry_table, correlation_table, normalization in (
        ("P", 7, 18, 0.052), ("D", 8, 19, 0.05),
    ):
        table = _load_table(asymmetry_table, root)
        headers = [item["name"] for item in table["headers"]]
        if headers != ["X", "Y", "Q**2 [GEV**2]",
                       "ASYM(NAME=PARALLEL,C=MEASURED)",
                       "ASYM(NAME=PARALLEL,C=BORN)"]:
            raise BornReferenceError("Asymmetry table column identities changed")
        rows = iter(enumerate(table["values"], start=1))
        all_points = []
        for index, edges in enumerate(Q2_EDGES):
            points = []
            for local_bin in range(1, len(edges)):
                global_bin, row = next(rows)
                x, y, q2 = row["x"]
                if [float(x["low"]), float(x["high"])] != DISPLAY_X_EDGES[index:index+2]:
                    raise BornReferenceError("Published x-cell identity changed")
                born = [cell for cell in row["y"] if cell["group"] == 1]
                if len(born) != 1:
                    raise BornReferenceError("Born asymmetry column missing")
                cell = born[0]
                point = {
                    "bin": local_bin, "global_bin": global_bin,
                    "x_mean": float(x.get("value", MISSING_X_MEANS.get(global_bin, math.nan))),
                    "y_mean": float(y["value"]), "q2_mean": float(q2["value"]),
                    "q2_low": edges[local_bin-1], "q2_high": edges[local_bin],
                    "value": float(cell["value"]), "stat": _uncertainty(cell, "stat"),
                    "systematic_combined": _uncertainty(cell, "sys"),
                    "systematics": {"published_total": _uncertainty(cell, "sys")},
                    "supported_prediction": edges[local_bin-1] >= 1.0,
                    "x_mean_source": "HEPData" if "value" in x else "paper Table XI/XII",
                }
                if not (edges[local_bin-1] < point["q2_mean"] < edges[local_bin]):
                    raise BornReferenceError("Published Q2 mean outside authoritative cell")
                if not all(math.isfinite(point[key]) for key in ("x_mean", "value", "stat")):
                    raise BornReferenceError("Non-finite asymmetry point")
                points.append(point)
                all_points.append(point)
            datasets.append({
                "id": f"Born{target}_X{index+1:02d}",
                "observable": f"A_parallel_{target}", "target": target,
                "selection": f"Born{target}_X{index+1:02d}", "plot_axis": "q2_mean",
                "rivet_path": f"/REF/{MEASUREMENT}/d{asymmetry_table:02d}-x{index+1:02d}-y02",
                "x_low": PHYSICAL_X_EDGES[index], "x_high": PHYSICAL_X_EDGES[index+1],
                "display_x_low": DISPLAY_X_EDGES[index], "display_x_high": DISPLAY_X_EDGES[index+1],
                "bin_edges": edges, "points": points,
                "units": {"x": "dimensionless", "q2_mean": "GeV2", "value": "dimensionless"},
                "provenance": {
                    **_source(asymmetry_table), "publication": PUBLICATION,
                    "publication_table": "Table XI" if target == "P" else "Table XII",
                    "normalization_uncertainty": normalization,
                    "normalization_uncertainty_included_in": "published systematic uncertainty",
                    "correlation_table": _source(correlation_table),
                    "correction_convention": "Published unfolded Born values retained unchanged",
                },
            })
        correlation = _correlation(_load_table(correlation_table, root))
        covariance[target] = {
            "global_bin_order": list(range(1, 46)), "correlation": correlation,
            "matrix": [[all_points[i]["stat"] * all_points[j]["stat"] * correlation[i][j]
                        for j in range(45)] for i in range(45)],
            "scope": "Statistical only; no systematic or joint proton/deuteron covariance supplied",
            "provenance": _source(correlation_table),
        }
    return {
        "schema_version": 2, "measurement": MEASUREMENT,
        "observable": "Born A_parallel", "datasets": datasets,
        "statistical_covariance": covariance,
        "source_files": [{"path": _source(i)["raw"], "sha256": SOURCES[i][1]} for i in SOURCES] + [
            {"path": str(DATA / "raw-born-apar-fig4.eps.gz"), "sha256": FIGURE_SHA256},
        ],
        "provenance": {"publication": PUBLICATION, "sources": [_source(i) for i in SOURCES],
                       "born_column": "group=1; measured group=0 is excluded"},
        "binning": {
            "physical_x_edges": PHYSICAL_X_EDGES, "display_x_edges": DISPLAY_X_EDGES,
            "q2_edges_by_x": Q2_EDGES,
            "internal_boundary_source": THESIS, "thesis_sha256": THESIS_SHA256,
            "internal_boundary_location": "Appendix C, printed page 97",
            "publication_boundary_cross_check": FIGURE, "figure_sha256": FIGURE_SHA256,
            "notes": "Thesis numerical internal edges agree with published Figure 4 within EPS coordinate quantization. Final publication removes the first thesis x cell and changes global Q2 support from 0.1-22 to 0.18-20 GeV2. Published x limits are rounded display labels; physical slices use thesis numerical x edges.",
        },
        "projection_reference_policy": "No derived phase-space-integrated data: published asymmetries and covariance do not provide unpolarized yield/cross-section weights. Both MC projections are prediction-only.",
    }


def validate_vendored(root: Path = ROOT) -> dict[str, Any]:
    expected = build_snapshot(root)
    actual = json.loads((root / REFERENCE_PATH).read_text(encoding="utf-8"))
    if actual != expected:
        raise BornReferenceError("Born reference snapshot differs from pinned source reconstruction")
    return actual


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", action="store_true", help="Regenerate normalized reference from pinned sources")
    args = parser.parse_args()
    if args.build:
        (ROOT / REFERENCE_PATH).write_text(json.dumps(build_snapshot(), indent=2) + "\n", encoding="utf-8")
    else:
        validate_vendored()
    print(f"Validated 90 Born asymmetries, 74 supported predictions and two 45x45 statistical covariances: {REFERENCE_PATH}")


if __name__ == "__main__":
    main()
