#!/usr/bin/env python3
"""Pinned STAR-jet and HERMES-SIDIS reference-data support.

This module is intentionally separate from the older phenomenology snapshots:
the STAR records contain covariance blocks spanning several observables, while
the HERMES source is an immutable APS ZIP archive rather than HEPData.  Public
entry points mirror ``phenomenology_reference_data`` so the campaign runner can
dispatch without special-casing the command-line workflow.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import math
import os
import shutil
import tempfile
import urllib.request
import zipfile
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


SCRIPT_PATH = Path(__file__).resolve()
DISPOL_ROOT = SCRIPT_PATH.parents[1]


class NewReferenceDataError(RuntimeError):
    """A pinned STAR/HERMES source or normalized snapshot failed validation."""


MEASUREMENTS = {
    "STAR_2021_I1850855",
    "STAR_2022_I1949588",
    "HERMES_2019_I1698889",
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
    "STAR_2021_I1850855": "analyses/rivet/pp/STAR_2021_I1850855.yoda.gz",
    "STAR_2022_I1949588": "analyses/rivet/pp/STAR_2022_I1949588.yoda.gz",
    "HERMES_2019_I1698889": "analyses/rivet/dis/HERMES_2019_I1698889.yoda.gz",
}

STAR_SPECIFICATIONS: dict[str, dict[str, Any]] = {
    "STAR_2021_I1850855": {
        "record_doi": "10.17182/hepdata.104836.v1",
        "record_id": 104836,
        "inspire_id": 1850855,
        "sqrt_s_gev": 200.0,
        "table_count": 21,
        "datasets": {
            "inclusive_forward": {"table": 4, "coordinate_group": 0, "value_group": 1},
            "inclusive_central": {"table": 5, "coordinate_group": 0, "value_group": 1},
            "dijet_same_sign": {"table": 6, "coordinate_group": 0, "value_group": 2},
            "dijet_opposite_sign": {"table": 7, "coordinate_group": 0, "value_group": 2},
            "inclusive_combined": {
                "table": 8, "coordinate_group": 0, "value_group": 2,
                "alternate_projection": True,
            },
        },
        "primary": [
            "inclusive_forward", "inclusive_central",
            "dijet_same_sign", "dijet_opposite_sign",
        ],
        "correlation_blocks": {
            9: ("inclusive_forward", "inclusive_forward"),
            10: ("inclusive_central", "inclusive_central"),
            11: ("inclusive_combined", "inclusive_combined"),
            12: ("inclusive_central", "inclusive_forward"),
            13: ("inclusive_forward", "dijet_same_sign"),
            14: ("inclusive_forward", "dijet_opposite_sign"),
            15: ("inclusive_central", "dijet_same_sign"),
            16: ("inclusive_central", "dijet_opposite_sign"),
            17: ("inclusive_combined", "dijet_same_sign"),
            18: ("inclusive_combined", "dijet_opposite_sign"),
            19: ("dijet_same_sign", "dijet_same_sign"),
            20: ("dijet_opposite_sign", "dijet_opposite_sign"),
            21: ("dijet_opposite_sign", "dijet_same_sign"),
        },
        "relative_luminosity_absolute": 7.0e-4,
        "polarization_relative": 0.061,
    },
    "STAR_2022_I1949588": {
        "record_doi": "10.17182/hepdata.114778.v1",
        "record_id": 114778,
        "inspire_id": 1949588,
        "sqrt_s_gev": 510.0,
        "table_count": 20,
        "datasets": {
            "inclusive": {"table": 1, "coordinate_group": 1, "value_group": 2},
            "dijet_A": {"table": 2, "coordinate_group": 0, "value_group": 1},
            "dijet_B": {"table": 3, "coordinate_group": 0, "value_group": 1},
            "dijet_C": {"table": 4, "coordinate_group": 0, "value_group": 1},
            "dijet_D": {"table": 5, "coordinate_group": 0, "value_group": 1},
        },
        "primary": ["inclusive", "dijet_A", "dijet_B", "dijet_C", "dijet_D"],
        "correlation_blocks": {
            6: ("inclusive", "inclusive"),
            7: ("inclusive", "dijet_A"),
            8: ("inclusive", "dijet_B"),
            9: ("inclusive", "dijet_C"),
            10: ("inclusive", "dijet_D"),
            11: ("dijet_A", "dijet_A"),
            12: ("dijet_A", "dijet_B"),
            13: ("dijet_A", "dijet_C"),
            14: ("dijet_A", "dijet_D"),
            15: ("dijet_B", "dijet_B"),
            16: ("dijet_B", "dijet_C"),
            17: ("dijet_B", "dijet_D"),
            18: ("dijet_C", "dijet_C"),
            19: ("dijet_C", "dijet_D"),
            20: ("dijet_D", "dijet_D"),
        },
        "relative_luminosity_absolute": 4.7e-4,
        "polarization_relative": 0.064,
    },
}

HERMES_ARCHIVE = "data/phenomenology/HERMES_2019_I1698889/raw-aps-supplement.zip"
HERMES_URL = (
    "https://journals.aps.org/prd/supplemental/"
    "10.1103/PhysRevD.99.112001/DB12134_supplemental.zip"
)
HERMES_ARCHIVE_SHA256 = (
    "05c38ee9dca04015248677e73a18f8cdb9c5df30a9ea6a6bbde15720a37aae45"
)
HERMES_X_EDGES = [0.023, 0.04, 0.055, 0.075, 0.1, 0.14, 0.2, 0.3, 0.4, 0.6]
HERMES_COARSE_X_EDGES = [0.023, 0.055, 0.1, 0.6]
HERMES_Z_EDGES = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
HERMES_PT_EDGES = [0.0, 0.15, 0.3, 0.4, 0.5, 0.6, 2.0]
HERMES_3D_X_EDGES = [0.023, 0.04, 0.055, 0.075, 0.1, 0.14, 0.2, 0.3, 0.4, 0.6]
HERMES_3D_Z_EDGES = [0.2, 0.35, 0.5, 0.8]
HERMES_3D_PT_EDGES = [0.0, 0.3, 0.5, 2.0]


def _hermes_mean_bin(value: float, edges: Sequence[float], label: str) -> int:
    """Locate a published mean in the discrete HERMES projection grid."""

    for index, (low, high) in enumerate(zip(edges, edges[1:])):
        if low <= value < high:
            return index
    raise NewReferenceDataError(
        f"HERMES {label} mean {value} is outside the published grid"
    )


def _hermes_xzpt_flat_bin(x: float, z: float, pt: float) -> int:
    """Map a sparse APS x-z-pT row to the complete 9x3x3 Rivet grid.

    The APS tables number only populated rows, starting at one, so their
    ``bin#`` column is not a canonical flat-grid coordinate.  The row means
    identify the corresponding discrete cells.  The main publication gives
    the final split as x=0.4, and the last-bin means cluster around x=0.44.
    The x=0.45 value repeated in the APS text-file headers is retained as a
    documented source discrepancy rather than used as the analysis boundary.
    """

    ix = _hermes_mean_bin(x, HERMES_3D_X_EDGES, "3D x")
    iz = _hermes_mean_bin(z, HERMES_3D_Z_EDGES, "3D z")
    ipt = _hermes_mean_bin(pt, HERMES_3D_PT_EDGES, "3D pT")
    return ix*9 + iz*3 + ipt


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise NewReferenceDataError(f"Could not read {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise NewReferenceDataError(f"Expected a JSON object in {path}")
    return value


def _dependent(row: Mapping[str, Any], group: int) -> Mapping[str, Any]:
    for cell in row.get("y", []):
        if int(cell.get("group", -1)) == group:
            return cell
    raise NewReferenceDataError(f"Missing dependent-variable group {group}")


def _error(cell: Mapping[str, Any], *labels: str) -> float:
    normalized = {label.lower() for label in labels}
    for item in cell.get("errors", []):
        if str(item.get("label", "")).lower() in normalized and "symerror" in item:
            token = item["symerror"]
            if isinstance(token, str) and token.endswith("%"):
                return float(token[:-1]) / 100.0
            return float(token)
    raise NewReferenceDataError(f"Missing symmetric uncertainty {labels}")


def _hepdata_table(payload: bytes, expected: Mapping[str, Any]) -> dict[str, Any]:
    if sha256_bytes(payload) != expected["sha256"]:
        raise NewReferenceDataError(
            f"Vendored checksum mismatch for {expected['path']}"
        )
    try:
        table = json.loads(payload.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise NewReferenceDataError(f"Invalid HEPData JSON: {exc}") from exc
    headers = [str(item.get("name")) for item in table.get("headers", [])]
    if str(table.get("name")) != expected["name"]:
        raise NewReferenceDataError(f"HEPData table name changed for {expected['path']}")
    if str(table.get("doi")) != expected["doi"]:
        raise NewReferenceDataError(f"HEPData table DOI changed for {expected['path']}")
    if headers != list(expected["headers"]):
        raise NewReferenceDataError(f"HEPData headers changed for {expected['path']}")
    if len(table.get("values", [])) != int(expected["rows"]):
        raise NewReferenceDataError(f"HEPData row count changed for {expected['path']}")
    return table


def _manifest(measurement: str) -> dict[str, Any]:
    path = DISPOL_ROOT / SOURCE_MANIFEST_PATHS[measurement]
    value = _json(path)
    if value.get("measurement") != measurement:
        raise NewReferenceDataError(f"Source manifest identity mismatch in {path}")
    return value


def _star_tables(measurement: str) -> dict[int, dict[str, Any]]:
    specification = STAR_SPECIFICATIONS[measurement]
    manifest = _manifest(measurement)
    if manifest.get("record_doi") != specification["record_doi"]:
        raise NewReferenceDataError(f"Pinned HEPData version changed for {measurement}")
    if int(manifest.get("record_version", -1)) != 1:
        raise NewReferenceDataError(f"Unexpected HEPData version for {measurement}")
    record_path = DISPOL_ROOT / str(manifest.get("record_path", ""))
    if (
        not record_path.is_file()
        or sha256_bytes(record_path.read_bytes())
        != str(manifest.get("record_sha256", ""))
    ):
        raise NewReferenceDataError(
            f"Vendored HEPData record checksum mismatch for {measurement}"
        )
    entries = manifest.get("tables", [])
    if len(entries) != specification["table_count"]:
        raise NewReferenceDataError(f"Expected all tables for {measurement}")
    expected_numbers = set(range(1, int(specification["table_count"]) + 1))
    numbers = {int(entry["number"]) for entry in entries}
    if numbers != expected_numbers:
        raise NewReferenceDataError(
            f"HEPData table inventory changed for {measurement}"
        )
    tables: dict[int, dict[str, Any]] = {}
    for entry in entries:
        number = int(entry["number"])
        payload = (DISPOL_ROOT / entry["path"]).read_bytes()
        tables[number] = _hepdata_table(payload, entry)
    return tables


def _star_analysis_binning(measurement: str) -> dict[str, list[float]] | None:
    """Load publication-only operational event bins when HEPData omits them.

    The STAR 510-GeV HEPData horizontal ranges are uncertainties on the
    unfolded parton-level coordinates.  The event intervals themselves are
    tabulated only in the publication, so they are stored in a separately
    checksum-pinned, reviewable snapshot.
    """

    if measurement != "STAR_2022_I1949588":
        return None
    manifest = _manifest(measurement)
    entry = manifest.get("analysis_binning")
    if not isinstance(entry, Mapping):
        raise NewReferenceDataError(
            f"Missing publication binning manifest for {measurement}"
        )
    path = DISPOL_ROOT / str(entry.get("path", ""))
    if (
        not path.is_file()
        or sha256_bytes(path.read_bytes()) != str(entry.get("sha256", ""))
    ):
        raise NewReferenceDataError(
            f"Publication binning checksum mismatch for {measurement}"
        )
    snapshot = _json(path)
    if snapshot.get("measurement") != measurement:
        raise NewReferenceDataError(
            f"Publication binning identity mismatch for {measurement}"
        )
    datasets = snapshot.get("datasets")
    expected = set(STAR_SPECIFICATIONS[measurement]["datasets"])
    if not isinstance(datasets, Mapping) or set(datasets) != expected:
        raise NewReferenceDataError(
            f"Publication binning inventory changed for {measurement}"
        )
    output: dict[str, list[float]] = {}
    for identifier, raw_edges in datasets.items():
        edges = [float(value) for value in raw_edges]
        point_count = len(
            _star_tables(measurement)[
                int(STAR_SPECIFICATIONS[measurement]["datasets"][identifier]["table"])
            ]["values"]
        )
        if (
            len(edges) != point_count + 1
            or any(first >= second for first, second in zip(edges, edges[1:]))
        ):
            raise NewReferenceDataError(
                f"Invalid publication binning for {measurement}/{identifier}"
            )
        output[str(identifier)] = edges
    return output


def _voronoi_edges(coordinates: Sequence[float]) -> list[float]:
    if len(coordinates) < 2:
        raise NewReferenceDataError("At least two coordinates are needed for binning")
    edges = [coordinates[0] - 0.5*(coordinates[1]-coordinates[0])]
    edges.extend(0.5*(a+b) for a, b in zip(coordinates[:-1], coordinates[1:]))
    edges.append(coordinates[-1] + 0.5*(coordinates[-1]-coordinates[-2]))
    if edges[0] < 0.0:
        edges[0] = 0.0
    return [float(value) for value in edges]


def _star_point_to_point_systematic(
    measurement: str, reported_systematic: float
) -> float:
    """Remove any separately profiled global term from a reported systematic.

    The STAR 510-GeV tables report a total systematic uncertainty that already
    contains the 4.7e-4 relative-luminosity contribution.  The accompanying
    point-to-point correlation matrices explicitly exclude that contribution,
    and the campaign profiles it as a global nuisance.  STAR 200-GeV HEPData
    supplies the point-to-point, luminosity, and polarization terms separately.
    """

    if measurement != "STAR_2022_I1949588":
        return reported_systematic
    relative_luminosity = float(
        STAR_SPECIFICATIONS[measurement]["relative_luminosity_absolute"]
    )
    point_to_point_variance = (
        reported_systematic**2 - relative_luminosity**2
    )
    if point_to_point_variance < -1.0e-18:
        raise NewReferenceDataError(
            "STAR 510-GeV reported systematic is smaller than its "
            "relative-luminosity component"
        )
    return math.sqrt(max(0.0, point_to_point_variance))


def _star_dataset(
    measurement: str,
    identifier: str,
    table: Mapping[str, Any],
    specification: Mapping[str, Any],
    analysis_binning: Mapping[str, Sequence[float]] | None = None,
) -> dict[str, Any]:
    points: list[dict[str, Any]] = []
    coordinate_group = int(specification["coordinate_group"])
    value_group = int(specification["value_group"])
    for index, row in enumerate(table["values"], start=1):
        coordinate = _dependent(row, coordinate_group)
        value = _dependent(row, value_group)
        stat = _error(value, "stat")
        reported_systematic = _error(value, "syst", "sys")
        systematic = _star_point_to_point_systematic(
            measurement, reported_systematic
        )
        point = {
            "point": index,
            "coordinate": float(coordinate["value"]),
            "coordinate_systematic": _error(coordinate, "syst", "sys"),
            "value": float(value["value"]),
            "stat": stat,
            "systematic": systematic,
            "systematic_combined": systematic,
            "point_to_point_total": math.hypot(stat, systematic),
            "uncertainties": {
                "statistical": stat,
                "point_to_point_systematic": systematic,
            },
        }
        if measurement == "STAR_2022_I1949588":
            point[
                "reported_systematic_including_relative_luminosity"
            ] = reported_systematic
        points.append(point)
    coordinates = [point["coordinate"] for point in points]
    display_edges = _voronoi_edges(coordinates)
    if measurement == "STAR_2021_I1850855":
        # Tables 1--3 define the detector intervals whose unfolded results
        # are quoted at the parton-level representative coordinates in
        # Tables 4--8.  The final inclusive interval (31.6--37.3 GeV) is part
        # of the published STAR interval sequence but has a zero JP1 yield
        # and is omitted from the Figure-1 yield table.
        if identifier.startswith("inclusive"):
            edges = [
                6.0, 7.1, 8.4, 9.9, 11.7, 13.8, 16.3, 19.2,
                22.7, 26.8, 31.6, 37.3,
            ]
            source_tables = [1, 2]
        else:
            edges = [17.0, 19.0, 23.0, 28.0, 34.0, 41.0, 58.0, 82.0]
            source_tables = [3]
        if len(edges) != len(points) + 1:
            raise NewReferenceDataError(
                f"Published STAR interval count does not match {identifier}"
            )
        for index, point in enumerate(points):
            point["low"], point["high"] = edges[index:index+2]
            point["coordinate_low"] = (
                point["coordinate"] - point["coordinate_systematic"]
            )
            point["coordinate_high"] = (
                point["coordinate"] + point["coordinate_systematic"]
            )
        binning_source = (
            "Published analysis intervals from HEPData Tables "
            + ", ".join(str(value) for value in source_tables)
            + "; unfolded points retain the parton-level representative "
              "coordinates from the asymmetry table"
        )
    elif measurement == "STAR_2022_I1949588":
        if analysis_binning is None or identifier not in analysis_binning:
            raise NewReferenceDataError(
                f"Missing publication event bins for {measurement}/{identifier}"
            )
        edges = [float(value) for value in analysis_binning[identifier]]
        coordinate_ranges = []
        for row in table["values"]:
            cell = row["x"][0]
            low, high = float(cell["low"]), float(cell["high"])
            if identifier == "inclusive":
                low *= 0.5*STAR_SPECIFICATIONS[measurement]["sqrt_s_gev"]
                high *= 0.5*STAR_SPECIFICATIONS[measurement]["sqrt_s_gev"]
            coordinate_ranges.append((low, high))
        if len(edges) != len(points) + 1:
            raise NewReferenceDataError(
                f"Publication event-bin count does not match {identifier}"
            )
        for index, (point, coordinate_range) in enumerate(
            zip(points, coordinate_ranges)
        ):
            point["low"], point["high"] = edges[index:index+2]
            point["coordinate_low"], point["coordinate_high"] = coordinate_range
            midpoint = 0.5 * sum(coordinate_range)
            tolerance = 5.0e-3 * max(1.0, abs(float(point["coordinate"])))
            if abs(midpoint - float(point["coordinate"])) > tolerance:
                raise NewReferenceDataError(
                    f"HEPData coordinate range no longer brackets "
                    f"{identifier} point {index + 1}"
                )
        binning_source = (
            "Published detector-level event intervals from Phys. Rev. D 105 "
            "(2022) 092011, Tables II/III and VI/VII; predictions are "
            "displayed at the unfolded parton-level HEPData coordinates"
        )
    else:
        raise NewReferenceDataError(f"Unknown STAR measurement {measurement}")
    observable = "inclusive_jet_pt" if identifier.startswith("inclusive") else "dijet_mass"
    return {
        "id": identifier,
        "observable": observable,
        "rivet_path": f"/{measurement}/{identifier}_ALL",
        "bin_edges": edges,
        "display_bin_edges": display_edges,
        "binning_source": binning_source,
        "coordinate_range_source": (
            "HEPData horizontal jet-energy-scale uncertainty"
        ),
        "points": points,
        "alternate_projection": bool(specification.get("alternate_projection", False)),
    }


def _coordinate_index(dataset: Mapping[str, Any], value: float) -> int:
    coordinates = [float(point["coordinate"]) for point in dataset["points"]]
    nearest = min(range(len(coordinates)), key=lambda index: abs(coordinates[index]-value))
    tolerance = 2.0e-5 * max(1.0, abs(value))
    if abs(coordinates[nearest]-value) > tolerance:
        raise NewReferenceDataError(
            f"Correlation coordinate {value} does not match {dataset['id']}"
        )
    return nearest


def _star_covariance(
    measurement: str,
    datasets: Mapping[str, Mapping[str, Any]],
    tables: Mapping[int, Mapping[str, Any]],
) -> dict[str, Any]:
    specification = STAR_SPECIFICATIONS[measurement]
    primary = list(specification["primary"])
    offsets: dict[str, int] = {}
    labels: list[str] = []
    points: list[Mapping[str, Any]] = []
    for identifier in primary:
        offsets[identifier] = len(points)
        dataset_points = list(datasets[identifier]["points"])
        points.extend(dataset_points)
        labels.extend(f"{identifier}:{index+1}" for index in range(len(dataset_points)))
    dimension = len(points)
    covariance = [[0.0 for _ in range(dimension)] for _ in range(dimension)]
    for index, point in enumerate(points):
        covariance[index][index] = float(point["point_to_point_total"])**2

    blocks: list[dict[str, Any]] = []
    for table_number, (first_id, second_id) in specification["correlation_blocks"].items():
        table = tables[table_number]
        basis = "point_to_point_total"
        if measurement == "STAR_2022_I1949588" and table_number >= 11:
            basis = "systematic"
        blocks.append({
            "table": table_number,
            "first": first_id,
            "second": second_id,
            "uncertainty_basis": basis,
            "rows": len(table["values"]),
        })
        if first_id not in offsets or second_id not in offsets:
            continue
        first_dataset, second_dataset = datasets[first_id], datasets[second_id]
        for row in table["values"]:
            first_coordinate = float(row["x"][0]["value"])
            second_coordinate = float(row["x"][1]["value"])
            rho = float(_dependent(row, 0)["value"])
            first_local = _coordinate_index(first_dataset, first_coordinate)
            second_local = _coordinate_index(second_dataset, second_coordinate)
            first_global = offsets[first_id] + first_local
            second_global = offsets[second_id] + second_local
            sigma_first = float(first_dataset["points"][first_local][basis])
            sigma_second = float(second_dataset["points"][second_local][basis])
            value = rho*sigma_first*sigma_second
            if first_global == second_global:
                # Statistical variance remains on the diagonal when a STAR
                # 510-GeV block is explicitly labelled systematic-only.
                value = float(points[first_global]["point_to_point_total"])**2
            covariance[first_global][second_global] = value
            covariance[second_global][first_global] = value

    correlation = [[0.0 for _ in range(dimension)] for _ in range(dimension)]
    for row in range(dimension):
        for column in range(dimension):
            denominator = math.sqrt(covariance[row][row]*covariance[column][column])
            correlation[row][column] = (
                covariance[row][column]/denominator if denominator else 0.0
            )
    return {
        "ordering": labels,
        "dimension": dimension,
        "covariance": covariance,
        "correlation": correlation,
        "blocks": blocks,
        "excludes_global_nuisances": True,
    }


def normalize_star(measurement: str) -> dict[str, Any]:
    tables = _star_tables(measurement)
    specification = STAR_SPECIFICATIONS[measurement]
    analysis_binning = _star_analysis_binning(measurement)
    datasets = {
        identifier: _star_dataset(
            measurement,
            identifier,
            tables[int(dataset_spec["table"])],
            dataset_spec,
            analysis_binning,
        )
        for identifier, dataset_spec in specification["datasets"].items()
    }
    covariance = _star_covariance(measurement, datasets, tables)
    if covariance["dimension"] != (36 if measurement.endswith("1850855") else 63):
        raise NewReferenceDataError(f"Unexpected primary point count for {measurement}")
    manifest = _manifest(measurement)
    return {
        "schema_version": 5,
        "measurement": measurement,
        "observable": "A_LL",
        "sqrt_s_gev": specification["sqrt_s_gev"],
        "primary_datasets": list(specification["primary"]),
        "datasets": datasets,
        "primary_covariance": covariance,
        "global_uncertainties": {
            "relative_luminosity_absolute":
                specification["relative_luminosity_absolute"],
            "polarization_relative": specification["polarization_relative"],
        },
        "goodness_of_fit": {
            "simultaneous_datasets": list(specification["primary"]),
            "excluded_alternate_projections": [
                identifier for identifier, dataset in datasets.items()
                if dataset["alternate_projection"]
            ],
            "global_nuisances_profiled_separately": True,
        },
        "provenance": {
            "source": "HEPData",
            "record_doi": specification["record_doi"],
            "record_id": specification["record_id"],
            "inspire_id": specification["inspire_id"],
            "version": 1,
            "retrieved": manifest["retrieved"],
            "source_manifest": SOURCE_MANIFEST_PATHS[measurement],
        },
    }


def _zip_members(payload: bytes) -> dict[str, bytes]:
    if payload[:2] != b"PK":
        prefix = payload[:200].decode("utf-8", errors="replace")
        raise NewReferenceDataError(
            "APS source is not a ZIP archive (HTML/Cloudflare responses are "
            f"rejected): {prefix!r}"
        )
    try:
        with zipfile.ZipFile(io.BytesIO(payload)) as archive:
            members = {
                name: archive.read(name)
                for name in archive.namelist()
                if not name.endswith("/") and not name.endswith(".DS_Store")
                and not name.startswith("__MACOSX/")
            }
    except zipfile.BadZipFile as exc:
        raise NewReferenceDataError(f"Invalid APS supplemental ZIP: {exc}") from exc
    if not members:
        raise NewReferenceDataError("APS supplemental ZIP contains no data")
    return members


def _numeric_rows(payload: bytes) -> list[list[float]]:
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise NewReferenceDataError(f"Non-UTF8 HERMES table: {exc}") from exc
    rows: list[list[float]] = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        try:
            rows.append([float(token) for token in stripped.split()])
        except ValueError as exc:
            raise NewReferenceDataError(f"Invalid HERMES numeric row {line!r}") from exc
    return rows


def _hermes_identity(filename: str) -> tuple[str, str]:
    token = Path(filename).name
    if token.startswith("pipos_"):
        species = "piplus"
    elif token.startswith("pineg_"):
        species = "piminus"
    elif token.startswith("kpos_"):
        species = "kplus"
    elif token.startswith("kneg_"):
        species = "kminus"
    elif token.startswith("pi_"):
        species = "pi_charge_difference"
    elif token.startswith("k_"):
        species = "k_charge_difference"
    elif token.startswith("h_"):
        species = "h_charge_difference"
    else:
        raise NewReferenceDataError(f"Unrecognized HERMES species in {token}")
    target = "proton" if "_pro" in token else "deuteron"
    return target, species


def _hermes_projection(member: str) -> str:
    relative = member.split("DC85_supplemental/", 1)[-1]
    directory = relative.split("/", 1)[0]
    return {
        "x": "x",
        "xz": "xz",
        "xPt": "xpt",
        "xzPt": "xzpt",
        "difference": "difference",
        "cosphi": "cosphi",
    }[directory]


def _matrix(rows: Sequence[Sequence[float]], expected: int, member: str) -> list[list[float]]:
    matrix = [[float(value) for value in row] for row in rows]
    if len(matrix) != expected or any(len(row) != expected for row in matrix):
        raise NewReferenceDataError(
            f"HERMES covariance dimension mismatch in {member}: "
            f"expected {expected}x{expected}"
        )
    for row in range(expected):
        if abs(matrix[row][row]-1.0) > 5.0e-4:
            raise NewReferenceDataError(f"Non-unit HERMES correlation diagonal in {member}")
        for column in range(expected):
            if abs(matrix[row][column]-matrix[column][row]) > 5.0e-4:
                raise NewReferenceDataError(f"Asymmetric HERMES correlation in {member}")
    return matrix


def _covariance_from_correlation(
    correlation: Sequence[Sequence[float]], errors: Sequence[float]
) -> list[list[float]]:
    return [
        [
            float(correlation[row][column])*float(errors[row])*float(errors[column])
            for column in range(len(errors))
        ]
        for row in range(len(errors))
    ]


def _hermes_point(projection: str, row: Sequence[float], point: int) -> dict[str, Any]:
    published_row: int | None = None
    if projection == "x":
        if len(row) != 11:
            raise NewReferenceDataError("HERMES x table column count changed")
        x, q2, z, pt, ratio, a1, a1stat, a1sys, apar, aparstat, aparsys = row
        bin_id = point
        w2 = None
    elif projection in {"xz", "xpt"}:
        if len(row) != 13:
            raise NewReferenceDataError(f"HERMES {projection} column count changed")
        (bin_id, x, q2, w2, z, pt, ratio, a1, a1stat, a1sys,
         apar, aparstat, aparsys) = row
        bin_id = int(bin_id)
    elif projection == "xzpt":
        if len(row) != 12:
            raise NewReferenceDataError("HERMES xzPt table column count changed")
        (bin_id, x, q2, z, pt, ratio, a1, a1stat, a1sys,
         apar, aparstat, aparsys) = row
        published_row = int(bin_id)
        bin_id = _hermes_xzpt_flat_bin(x, z, pt)
        w2 = None
    else:
        raise NewReferenceDataError(f"Not an A1/Aparallel projection: {projection}")
    result = {
        "point": point+1,
        "flat_bin": int(bin_id),
        "x": x,
        "q2": q2,
        "w2": w2,
        "z": z,
        "pt": pt,
        "a1_over_aparallel": ratio,
        "a1": a1,
        "a1_stat": a1stat,
        "a1_systematic": a1sys,
        "aparallel": apar,
        "aparallel_stat": aparstat,
        "aparallel_systematic": aparsys,
        "systematic_combined": aparsys,
    }
    if published_row is not None:
        result["published_row"] = published_row
    return result


def _hermes_dataset(
    member: str, payload: bytes, members: Mapping[str, bytes]
) -> tuple[str, dict[str, Any]]:
    projection = _hermes_projection(member)
    target, species = _hermes_identity(member)
    rows = _numeric_rows(payload)
    identifier = f"{target}_{species}_{projection}"
    if projection in {"x", "xz", "xpt", "xzpt"}:
        points = [_hermes_point(projection, row, index) for index, row in enumerate(rows)]
        if projection == "xzpt":
            flat_bins = [int(point["flat_bin"]) for point in points]
            if len(set(flat_bins)) != len(flat_bins):
                raise NewReferenceDataError(
                    f"HERMES sparse 3D rows do not map uniquely in {member}"
                )
        correlation_name = member.replace("/data_tables/", "/correlations/")
        replacements = {
            "_1D_A1Apar.txt": "_1D_corr.txt",
            "_2D-xz_A1Apar.txt": "_2D-xz_corr.txt",
            "_2D-xPt_A1Apar.txt": "_2D-xPt_corr.txt",
            "_3D_A1Apar.txt": "_3D_corr.txt",
            "_3D_A1Apartxt.txt": "_3D_corr.txt",
        }
        for suffix, replacement in replacements.items():
            if correlation_name.endswith(suffix):
                correlation_name = correlation_name[:-len(suffix)] + replacement
                break
        if correlation_name in members:
            correlation = _matrix(
                _numeric_rows(members[correlation_name]), len(points), correlation_name
            )
            status = "published"
            a1_covariance = _covariance_from_correlation(
                correlation, [point["a1_stat"] for point in points]
            )
            aparallel_covariance = _covariance_from_correlation(
                correlation, [point["aparallel_stat"] for point in points]
            )
        else:
            correlation = None
            a1_covariance = None
            aparallel_covariance = None
            status = "not supplied in official APS archive"
        if projection == "x":
            flat_edges = HERMES_X_EDGES
        elif projection == "xz":
            flat_edges = list(range(22))
        elif projection == "xpt":
            flat_edges = list(range(19))
        else:
            flat_edges = list(range(82))
        return identifier, {
            "id": identifier,
            "target": target,
            "species": species,
            "projection": projection,
            "observable": "A_parallel",
            "rivet_path": f"/HERMES_2019_I1698889/{identifier}",
            "flat_bin_edges": flat_edges,
            "points": points,
            "correlation": correlation,
            "a1_statistical_covariance": a1_covariance,
            "aparallel_statistical_covariance": aparallel_covariance,
            "covariance_status": status,
            "covariance_scope": (
                "within this projection only; no cross-projection covariance published"
            ),
        }
    if projection == "difference":
        if any(len(row) != 5 for row in rows):
            raise NewReferenceDataError("HERMES charge-difference column count changed")
        points = [{
            "point": index+1,
            "x": row[0],
            "q2": row[1],
            "a1": row[2],
            "a1_stat": row[3],
            "a1_systematic": row[4],
            "systematic_combined": row[4],
        } for index, row in enumerate(rows)]
        return identifier, {
            "id": identifier,
            "target": target,
            "species": species,
            "projection": projection,
            "observable": "A1_charge_difference",
            "rivet_path": f"/HERMES_2019_I1698889/{identifier}",
            "bin_edges": HERMES_X_EDGES,
            "points": points,
            "covariance_status": "not supplied in official APS archive",
            "model_note": (
                "The supplement publishes charge-difference A1, not A_parallel; "
                "comparison requires the published R1999 depolarization conversion."
            ),
        }
    if projection == "cosphi":
        if any(len(row) != 10 for row in rows):
            raise NewReferenceDataError("HERMES cos(phi) column count changed")
        points = [{
            "point": index+1,
            "binning": int(row[0]),
            "target_code": int(row[1]),
            "hadron_code": int(row[2]),
            "x": row[3],
            "q2": row[4],
            "z": row[5],
            "pt": row[6],
            "aparallel_cosphi_amplitude": row[7],
            "stat": row[8],
            "systematic": row[9],
            "systematic_combined": row[9],
        } for index, row in enumerate(rows)]
        return identifier, {
            "id": identifier,
            "target": target,
            "species": species,
            "projection": projection,
            "observable": "A_parallel_cosphi_amplitude",
            "source_column_label": "2<cos(phi)>",
            "rivet_path": (
                f"/HERMES_2019_I1698889/PUBLISHED_AParallelCosPhi/"
                f"{identifier}"
            ),
            "points": points,
            "diagnostic_only": True,
            "comparison_status": (
                "Published spin-asymmetry fit amplitude; no compatible "
                "generator prediction is currently implemented."
            ),
        }
    raise NewReferenceDataError(f"Unsupported HERMES projection {projection}")


def _validated_hermes_members(payload: bytes) -> tuple[dict[str, bytes], dict[str, Any]]:
    if sha256_bytes(payload) != HERMES_ARCHIVE_SHA256:
        raise NewReferenceDataError(
            "HERMES APS supplement checksum changed: expected "
            f"{HERMES_ARCHIVE_SHA256}, got {sha256_bytes(payload)}"
        )
    members = _zip_members(payload)
    manifest = _manifest("HERMES_2019_I1698889")
    expected = {entry["path"]: entry for entry in manifest.get("members", [])}
    if set(members) != set(expected):
        missing = sorted(set(expected)-set(members))
        added = sorted(set(members)-set(expected))
        raise NewReferenceDataError(
            f"HERMES supplemental inventory changed; missing={missing}, added={added}"
        )
    for name, member_payload in members.items():
        entry = expected[name]
        if sha256_bytes(member_payload) != entry["sha256"]:
            raise NewReferenceDataError(f"HERMES member checksum changed: {name}")
        rows = _numeric_rows(member_payload)
        if len(rows) != int(entry["rows"]):
            raise NewReferenceDataError(f"HERMES row count changed: {name}")
        if any(len(row) != int(entry["columns"]) for row in rows):
            raise NewReferenceDataError(f"HERMES column count changed: {name}")
    return members, manifest


def normalize_hermes() -> dict[str, Any]:
    payload = (DISPOL_ROOT/HERMES_ARCHIVE).read_bytes()
    members, manifest = _validated_hermes_members(payload)
    datasets: dict[str, Any] = {}
    diagnostics: dict[str, Any] = {}
    for name in sorted(members):
        if "/data_tables/" not in name:
            continue
        identifier, dataset = _hermes_dataset(name, members[name], members)
        target = diagnostics if dataset.get("diagnostic_only") else datasets
        if identifier in target:
            raise NewReferenceDataError(f"Duplicate HERMES dataset {identifier}")
        target[identifier] = dataset
    projection_counts: dict[str, int] = {}
    for dataset in datasets.values():
        projection = str(dataset["projection"])
        projection_counts[projection] = (
            projection_counts.get(projection, 0) + len(dataset["points"])
        )
    expected_counts = {"x": 54, "xz": 126, "xpt": 108, "xzpt": 477, "difference": 36}
    if projection_counts != expected_counts:
        raise NewReferenceDataError(
            f"HERMES projection inventory changed: {projection_counts}"
        )
    if sum(len(dataset["points"]) for dataset in diagnostics.values()) != 180:
        raise NewReferenceDataError("HERMES azimuthal-moment inventory changed")
    return {
        "schema_version": 5,
        "measurement": "HERMES_2019_I1698889",
        "observable": "identified-hadron A_parallel",
        "beam_energy_gev": 27.6,
        "binning": {
            "x": HERMES_X_EDGES,
            "coarse_x": HERMES_COARSE_X_EDGES,
            "z": HERMES_Z_EDGES,
            "pt_gev": HERMES_PT_EDGES,
            "three_dimensional_x": HERMES_3D_X_EDGES,
            "three_dimensional_z": HERMES_3D_Z_EDGES,
            "three_dimensional_pt_gev": HERMES_3D_PT_EDGES,
            "azimuthal_x_fine": [0.023, 0.040, 0.055, 0.075, 0.140, 0.600],
            "azimuthal_x_coarse": [0.023, 0.100, 0.600],
            "azimuthal_z_fine": [0.20, 0.32, 0.44, 0.56, 0.68, 0.80],
            "azimuthal_z_coarse": [0.20, 0.40, 0.80],
            "azimuthal_pt_gev": [0.0, 0.30, 0.40, 0.50, 0.60, 2.0],
        },
        "selection": {
            "q2_min_gev2": 1.0,
            "w2_min_gev2": 10.0,
            "y_max": 0.85,
            "lepton_theta_rad": [0.04, 0.22],
            "xf_min": 0.1,
            "z": [0.2, 0.8],
            "xz_extended_z": [0.1, 0.8],
            "proton_pion_momentum_gev": [4.0, 13.8],
            "deuteron_hadron_momentum_gev": [2.0, 15.0],
        },
        "datasets": datasets,
        "diagnostics": diagnostics,
        "projection_point_counts": projection_counts,
        "deuteron_impulse_approximation": {
            "unpolarized": "(p+n)/2",
            "longitudinal": "(p+n)/2",
            "published_nucleon_to_nucleus_polarization_ratio": 0.926,
            "comparison_note": (
                "The published deuteron A_parallel values already include "
                "the 1/f_D nucleon-polarization correction, so no additional "
                "0.926 factor is applied to free-p+n theory."
            ),
        },
        "goodness_of_fit": {
            "combine_projections": False,
            "reason": "No cross-projection covariance is published",
            "missing_covariance_dataset": "deuteron_kminus_xpt",
        },
        "global_uncertainties": {
            "proton_polarization_relative": 0.066,
            "deuteron_polarization_relative": 0.057,
            "treatment": (
                "Recorded as publication metadata. The APS tables provide "
                "only the combined point systematic, so these components "
                "are not profiled separately to avoid double counting."
            ),
        },
        "provenance": {
            "source": "official APS supplemental archive",
            "url": HERMES_URL,
            "archive_sha256": HERMES_ARCHIVE_SHA256,
            "publication_doi": "10.1103/PhysRevD.99.112001",
            "arxiv": "1810.07054",
            "retrieved": manifest["retrieved"],
            "source_manifest": SOURCE_MANIFEST_PATHS["HERMES_2019_I1698889"],
            "normalization_overrides": [
                {
                    "field": "binning.three_dimensional_x[-2]",
                    "normalized_value": 0.4,
                    "supplemental_header_value": 0.45,
                    "reason": (
                        "The publication binning table gives 0.4 and the "
                        "reported last-bin means cluster near x=0.44; the "
                        "0.45 value in the APS text-file headers is treated "
                        "as a header typo."
                    ),
                }
            ],
        },
    }


def normalized_from_raw(measurement: str) -> dict[str, Any]:
    if measurement in STAR_SPECIFICATIONS:
        return normalize_star(measurement)
    if measurement == "HERMES_2019_I1698889":
        return normalize_hermes()
    raise NewReferenceDataError(f"Unknown STAR/HERMES measurement {measurement}")


def assert_normalized_snapshot_matches(
    generated: Any,
    snapshot: Any,
    path: str = "$",
) -> None:
    """Compare normalized records across platforms without accepting drift.

    The STAR covariance contains square roots whose final binary digit can
    differ between otherwise compatible libm implementations.  Structure,
    strings, integers, and booleans remain exact; floating-point values allow
    only a few units in the last place.  Physics-scale changes therefore
    remain hard failures while macOS/Linux roundoff does not stale a pinned
    snapshot.
    """

    if type(generated) is not type(snapshot):
        raise NewReferenceDataError(
            f"Normalized snapshot type mismatch at {path}: "
            f"{type(generated).__name__} != {type(snapshot).__name__}"
        )
    if isinstance(generated, Mapping):
        if set(generated) != set(snapshot):
            difference = sorted(set(generated) ^ set(snapshot))
            raise NewReferenceDataError(
                f"Normalized snapshot key mismatch at {path}: {difference}"
            )
        for key in generated:
            assert_normalized_snapshot_matches(
                generated[key], snapshot[key], f"{path}.{key}"
            )
        return
    if isinstance(generated, list):
        if len(generated) != len(snapshot):
            raise NewReferenceDataError(
                f"Normalized snapshot length mismatch at {path}: "
                f"{len(generated)} != {len(snapshot)}"
            )
        for index, (actual, expected) in enumerate(zip(generated, snapshot)):
            assert_normalized_snapshot_matches(
                actual, expected, f"{path}[{index}]"
            )
        return
    if isinstance(generated, float):
        if not math.isclose(
            generated, snapshot, rel_tol=1.0e-15, abs_tol=1.0e-24
        ):
            raise NewReferenceDataError(
                f"Normalized snapshot float mismatch at {path}: "
                f"{generated!r} != {snapshot!r}"
            )
        return
    if generated != snapshot:
        raise NewReferenceDataError(
            f"Normalized snapshot mismatch at {path}: "
            f"{generated!r} != {snapshot!r}"
        )


def _validate_psd(matrix: Sequence[Sequence[float]], label: str) -> None:
    try:
        import numpy as np
    except ImportError as exc:
        raise NewReferenceDataError("NumPy is required for covariance validation") from exc
    array = np.asarray(matrix, dtype=float)
    if not np.allclose(array, array.T, rtol=0.0, atol=1.0e-12):
        raise NewReferenceDataError(f"{label} is not symmetric")
    eigenvalues = np.linalg.eigvalsh(array)
    tolerance = max(1.0, float(np.max(np.abs(eigenvalues))))*1.0e-10
    if float(eigenvalues[0]) < -tolerance:
        raise NewReferenceDataError(
            f"{label} is not positive semidefinite; min eigenvalue={eigenvalues[0]}"
        )


def validate_vendored(measurement: str) -> dict[str, Any]:
    if measurement not in MEASUREMENTS:
        raise NewReferenceDataError(f"Unknown STAR/HERMES measurement {measurement}")
    generated = normalized_from_raw(measurement)
    snapshot = _json(DISPOL_ROOT/REFERENCE_PATHS[measurement])
    assert_normalized_snapshot_matches(generated, snapshot)
    if measurement in STAR_SPECIFICATIONS:
        _validate_psd(
            snapshot["primary_covariance"]["covariance"],
            f"{measurement} primary covariance",
        )
    else:
        for dataset in snapshot["datasets"].values():
            covariance = dataset.get("aparallel_statistical_covariance")
            if covariance is not None:
                _validate_psd(covariance, dataset["id"])
    return snapshot


def fetch_and_validate(
    measurement: str,
    cache_directory: Path | None = None,
    source_file: Path | None = None,
) -> list[Path]:
    validate_vendored(measurement)
    cache = cache_directory or (
        DISPOL_ROOT/"campaigns"/"phenomenology"/"_data_cache"/measurement
    )
    cache.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    if measurement == "HERMES_2019_I1698889":
        if source_file is not None:
            payload = source_file.expanduser().resolve().read_bytes()
        else:
            request = urllib.request.Request(
                HERMES_URL,
                headers={"User-Agent": "HerwigPol-reference-data/1.0"},
            )
            try:
                with urllib.request.urlopen(request, timeout=60) as response:
                    payload = response.read()
            except OSError as exc:
                raise NewReferenceDataError(
                    f"Could not download APS archive ({exc}); use "
                    "fetch-data --source-file PATH"
                ) from exc
        _validated_hermes_members(payload)
        if payload != (DISPOL_ROOT/HERMES_ARCHIVE).read_bytes():
            raise NewReferenceDataError(
                "Checksum-valid APS payload differs from vendored archive"
            )
        destination = cache/"source-01.zip"
        destination.write_bytes(payload)
        outputs.append(destination)
    else:
        manifest = _manifest(measurement)
        for entry in manifest["tables"]:
            request = urllib.request.Request(
                entry["url"], headers={"User-Agent": "HerwigPol-reference-data/1.0"}
            )
            try:
                with urllib.request.urlopen(request, timeout=60) as response:
                    payload = response.read()
            except OSError as exc:
                raise NewReferenceDataError(
                    f"Could not download {entry['url']}: {exc}"
                ) from exc
            _hepdata_table(payload, entry)
            vendored = (DISPOL_ROOT/entry["path"]).read_bytes()
            if payload != vendored:
                raise NewReferenceDataError(
                    f"Checksum-valid HEPData table differs from {entry['path']}"
                )
            destination = cache/f"table-{int(entry['number']):02d}.json"
            destination.write_bytes(payload)
            outputs.append(destination)
    write_reference_yoda(measurement)
    return outputs


def _import_yoda() -> Any:
    try:
        import yoda  # type: ignore
    except (ImportError, OSError) as exc:
        raise NewReferenceDataError(
            "YODA Python bindings are unavailable; load herwig/pol"
        ) from exc
    return yoda


def _write_yoda(yoda: Any, objects: Sequence[Any], destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="polarized-data-", dir=destination.parent) as temporary:
        plain = Path(temporary)/destination.with_suffix("").name
        yoda.write(list(objects), str(plain))
        if destination.suffix == ".gz":
            compressed = Path(temporary)/destination.name
            with plain.open("rb") as source, compressed.open("wb") as target:
                with gzip.GzipFile(
                    filename="", mode="wb", fileobj=target, mtime=0
                ) as stream:
                    shutil.copyfileobj(source, stream)
            os.replace(compressed, destination)
        else:
            os.replace(plain, destination)


def _scatter(
    yoda: Any,
    path: str,
    points: Iterable[
        tuple[float, float, float] | tuple[float, float, float, float]
    ],
    annotations: Mapping[str, Any],
) -> Any:
    result = yoda.Scatter2D("/REF" + path)
    result.setAnnotation("IsRef", 1)
    for key, value in annotations.items():
        result.setAnnotation(str(key), value)
    for point in points:
        if len(point) == 3:
            x, value, total_error = point
            x_error = 0.0
        else:
            x, value, total_error, x_error = point
        result.addPoint(
            float(x), float(value), float(x_error), float(total_error)
        )
    return result


def write_reference_yoda(
    measurement: str, snapshot: Mapping[str, Any] | None = None
) -> Path:
    snapshot = snapshot or validate_vendored(measurement)
    yoda = _import_yoda()
    objects: list[Any] = []
    if measurement in STAR_SPECIFICATIONS:
        for dataset in snapshot["datasets"].values():
            objects.append(_scatter(
                yoda,
                dataset["rivet_path"],
                (
                    (
                        point["coordinate"],
                        point["value"],
                        math.hypot(point["stat"], point["systematic"]),
                        point["coordinate_systematic"],
                    )
                    for point in dataset["points"]
                ),
                {
                    "Observable": "A_LL",
                    "ExperimentalError": "statistical and point-to-point systematic",
                    "GlobalUncertainties": json.dumps(
                        snapshot["global_uncertainties"], sort_keys=True
                    ),
                    "PublishedCoordinates": 1,
                    "AlternateProjection": int(dataset["alternate_projection"]),
                },
            ))
    else:
        for dataset in snapshot["datasets"].values():
            observable = dataset["observable"]
            if dataset["projection"] == "difference":
                point_rows = (
                    (
                        point["x"],
                        point["a1"],
                        math.hypot(point["a1_stat"], point["a1_systematic"]),
                    )
                    for point in dataset["points"]
                )
            else:
                point_rows = (
                    (
                        point["x"] if dataset["projection"] == "x"
                        else point["flat_bin"]+0.5,
                        point["aparallel"],
                        math.hypot(
                            point["aparallel_stat"],
                            point["aparallel_systematic"],
                        ),
                    )
                    for point in dataset["points"]
                )
            objects.append(_scatter(
                yoda, dataset["rivet_path"], point_rows,
                {
                    "Observable": observable,
                    "Target": dataset["target"],
                    "Species": dataset["species"],
                    "Projection": dataset["projection"],
                    "CovarianceStatus": dataset["covariance_status"],
                },
            ))
        for dataset in snapshot["diagnostics"].values():
            objects.append(_scatter(
                yoda,
                dataset["rivet_path"],
                (
                    (
                        point["point"]-0.5,
                        point["aparallel_cosphi_amplitude"],
                        math.hypot(point["stat"], point["systematic"]),
                    )
                    for point in dataset["points"]
                ),
                {
                    "Observable": "A_parallel_cosphi_amplitude",
                    "SourceColumnLabel": "2<cos(phi)>",
                    "CompatibleGeneratorPrediction": 0,
                    "DiagnosticOnly": 1,
                    "Target": dataset["target"],
                    "Species": dataset["species"],
                },
            ))
    destination = DISPOL_ROOT/REFERENCE_YODA_PATHS[measurement]
    _write_yoda(yoda, objects, destination)
    return destination


def write_normalized_snapshot(measurement: str) -> Path:
    """Regenerate a normalized snapshot from checksum-pinned vendored bytes."""

    if measurement in STAR_SPECIFICATIONS:
        snapshot = normalize_star(measurement)
    elif measurement == "HERMES_2019_I1698889":
        snapshot = normalize_hermes()
    else:
        raise NewReferenceDataError(f"Unknown measurement {measurement}")
    destination = DISPOL_ROOT / REFERENCE_PATHS[measurement]
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload = (
        json.dumps(snapshot, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    with tempfile.NamedTemporaryFile(
        prefix=destination.name + ".",
        dir=destination.parent,
        delete=False,
    ) as stream:
        temporary = Path(stream.name)
        stream.write(payload)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, destination)
    return destination


def build_source_manifest_from_vendored(measurement: str) -> dict[str, Any]:
    """Build the initial reviewable source manifest from already-vendored bytes.

    This is a repository-maintenance helper, not a refresh operation.  Normal
    campaign commands only read and enforce the resulting manifest.
    """
    if measurement in STAR_SPECIFICATIONS:
        specification = STAR_SPECIFICATIONS[measurement]
        directory = DISPOL_ROOT/f"data/phenomenology/{measurement}/raw"
        record_path = directory/"record.json"
        record_payload = record_path.read_bytes()
        record = json.loads(record_payload)
        if record["record"]["hepdata_doi"] != specification["record_doi"]:
            raise NewReferenceDataError("Wrong HEPData record supplied")
        table_metadata = {
            str(item["doi"]): item for item in record["data_tables"]
        }
        tables = []
        for number in range(1, int(specification["table_count"])+1):
            path = directory/f"t{number:02d}.json"
            payload = path.read_bytes()
            table = json.loads(payload)
            metadata = table_metadata[str(table["doi"])]
            tables.append({
                "number": number,
                "path": str(path.relative_to(DISPOL_ROOT)),
                "url": (
                    f"https://www.hepdata.net/record/data/"
                    f"{specification['record_id']}/{metadata['id']}/1"
                ),
                "sha256": sha256_bytes(payload),
                "name": table["name"],
                "doi": table["doi"],
                "headers": [str(item["name"]) for item in table["headers"]],
                "rows": len(table["values"]),
            })
        result = {
            "schema_version": 1,
            "measurement": measurement,
            "source": "HEPData",
            "record_doi": specification["record_doi"],
            "record_version": 1,
            "record_path": str(record_path.relative_to(DISPOL_ROOT)),
            "record_sha256": sha256_bytes(record_payload),
            "retrieved": "2026-07-22",
            "tables": tables,
        }
        if measurement == "STAR_2022_I1949588":
            binning_path = (
                DISPOL_ROOT
                / "data/phenomenology/STAR_2022_I1949588/"
                  "analysis-binning.json"
            )
            result["analysis_binning"] = {
                "path": str(binning_path.relative_to(DISPOL_ROOT)),
                "publication_doi": "10.1103/PhysRevD.105.092011",
                "sha256": sha256_bytes(binning_path.read_bytes()),
            }
        return result
    if measurement == "HERMES_2019_I1698889":
        archive = DISPOL_ROOT/HERMES_ARCHIVE
        payload = archive.read_bytes()
        if sha256_bytes(payload) != HERMES_ARCHIVE_SHA256:
            raise NewReferenceDataError("Wrong APS archive supplied")
        members = _zip_members(payload)
        entries = []
        for name, member_payload in sorted(members.items()):
            rows = _numeric_rows(member_payload)
            entries.append({
                "path": name,
                "sha256": sha256_bytes(member_payload),
                "bytes": len(member_payload),
                "rows": len(rows),
                "columns": len(rows[0]) if rows else 0,
            })
        return {
            "schema_version": 1,
            "measurement": measurement,
            "source": "official APS supplemental archive",
            "publication_doi": "10.1103/PhysRevD.99.112001",
            "url": HERMES_URL,
            "archive_path": HERMES_ARCHIVE,
            "archive_sha256": HERMES_ARCHIVE_SHA256,
            "retrieved": "2026-07-22",
            "members": entries,
        }
    raise NewReferenceDataError(f"Unknown measurement {measurement}")
