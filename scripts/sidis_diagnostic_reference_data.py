#!/usr/bin/env python3
"""Frozen paper-source data for the SIDIS diagnostic additions.

COMPASS 2013 and 2014 publish their complete numerical tables in the arXiv
sources.  This module reconstructs the normalized snapshots directly from
those checksum-pinned sources.  HERMES 2013 is intentionally absent here:
only its exact observable definition is tracked because the official 21,600
moment values and covariance endpoint were unavailable and no values are
invented.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import os
import re
import shutil
import tarfile
import tempfile
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]


class SIDISDiagnosticDataError(RuntimeError):
    """A pinned diagnostic source or normalized snapshot failed validation."""


MEASUREMENTS = {"COMPASS_2013_I1236358", "COMPASS_2014_I1278730"}
REFERENCE_PATHS = {
    item: f"data/phenomenology/{item}/reference.json" for item in MEASUREMENTS
}
REFERENCE_YODA_PATHS = {
    item: f"analyses/rivet/dis/{item}.yoda.gz" for item in MEASUREMENTS
}
SOURCE_MANIFEST_PATHS = {
    item: f"data/phenomenology/{item}/source-manifest.json" for item in MEASUREMENTS
}

Z_EDGES_2013 = [.20, .25, .30, .35, .40, .50, .60, .70, .80]
PT2_EDGES_2013 = [.01 + .02 * index for index in range(36)] + [.7225]

X_EDGES_2014 = [.003, .008, .013, .020, .032, .050, .080, .130]
Z_EDGES_2014 = [.20, .25, .30, .34, .38, .42, .49, .63, .85]
PT_EDGES_2014 = [.10, .20, .27, .33, .39, .46, .55, .64, .77, 1.00]
X3_EDGES_2014 = [.003, .012, .020, .038, .130]
Z3_EDGES_2014 = [.20, .25, .32, .40, .55, .70, .85]
PT3_EDGES_2014 = [.10, .30, .50, .64, 1.00]


def _load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SIDISDiagnosticDataError(f"Could not read {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise SIDISDiagnosticDataError(f"Expected a JSON object in {path}")
    return value


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _source_text(measurement: str) -> tuple[str, Mapping[str, Any]]:
    manifest = _load_json(ROOT / SOURCE_MANIFEST_PATHS[measurement])
    if manifest.get("measurement") != measurement:
        raise SIDISDiagnosticDataError(
            f"Source manifest identity mismatch for {measurement}"
        )
    source = manifest.get("source")
    if not isinstance(source, Mapping):
        raise SIDISDiagnosticDataError(f"Missing source entry for {measurement}")
    path = ROOT / str(source["path"])
    payload = path.read_bytes()
    if len(payload) != int(source["bytes"]) or _sha256(payload) != source["sha256"]:
        raise SIDISDiagnosticDataError(f"Source checksum mismatch for {measurement}")
    member = str(source["member"])
    try:
        with tarfile.open(path, "r:*") as archive:
            stream = archive.extractfile(member)
            if stream is None:
                raise KeyError(member)
            text = stream.read().decode("utf-8", errors="strict")
    except (tarfile.TarError, KeyError, UnicodeDecodeError) as exc:
        raise SIDISDiagnosticDataError(
            f"Could not extract {member} from {path}: {exc}"
        ) from exc
    return text, manifest


def _table(text: str, label: str) -> str:
    anchor = rf"\label{{{label}}}"
    position = text.find(anchor)
    if position < 0:
        raise SIDISDiagnosticDataError(f"Could not find table label {label}")
    begin_tokens = (
        r"\begin{table}", r"\begin{table*}",
        r"\begin{sidewaystable}", r"\begin{sidewaystable*}",
    )
    start = max(text.rfind(token, 0, position) for token in begin_tokens)
    if start < 0:
        raise SIDISDiagnosticDataError(f"Could not find start of table {label}")
    ends = [
        found for token in (
            r"\end{table}", r"\end{table*}",
            r"\end{sidewaystable}", r"\end{sidewaystable*}",
        )
        if (found := text.find(token, position)) >= 0
    ]
    if not ends:
        raise SIDISDiagnosticDataError(f"Could not find end of table {label}")
    return text[start:min(ends)]


def _atomic_json(path: Path, payload: Mapping[str, Any]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    encoded = (
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False,
                   allow_nan=False) + "\n"
    ).encode("utf-8")
    with tempfile.NamedTemporaryFile(
        prefix=path.name + ".", dir=path.parent, delete=False
    ) as stream:
        temporary = Path(stream.name)
        stream.write(encoded)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)
    return path


def _assert_equal(generated: Any, snapshot: Any, path: str = "$") -> None:
    if type(generated) is not type(snapshot):
        raise SIDISDiagnosticDataError(f"Snapshot type mismatch at {path}")
    if isinstance(generated, Mapping):
        if set(generated) != set(snapshot):
            raise SIDISDiagnosticDataError(f"Snapshot key mismatch at {path}")
        for key in generated:
            _assert_equal(generated[key], snapshot[key], f"{path}.{key}")
    elif isinstance(generated, list):
        if len(generated) != len(snapshot):
            raise SIDISDiagnosticDataError(f"Snapshot length mismatch at {path}")
        for index, (left, right) in enumerate(zip(generated, snapshot)):
            _assert_equal(left, right, f"{path}[{index}]")
    elif isinstance(generated, float):
        if not math.isclose(generated, snapshot, rel_tol=1.e-14, abs_tol=1.e-24):
            raise SIDISDiagnosticDataError(
                f"Snapshot float mismatch at {path}: {generated} != {snapshot}"
            )
    elif generated != snapshot:
        raise SIDISDiagnosticDataError(
            f"Snapshot mismatch at {path}: {generated!r} != {snapshot!r}"
        )


def _parenthesized_value(field: str) -> tuple[float, float]:
    match = re.search(
        r"([+-]?\s*(?:\d+\.\d*|\.\d+|\d+))\s*\(\s*(\d+)\s*\)",
        field,
    )
    if match is None:
        raise SIDISDiagnosticDataError(f"Could not parse value(error) field {field!r}")
    token = match.group(1).replace(" ", "")
    decimals = len(token.split(".", 1)[1]) if "." in token else 0
    return float(token), int(match.group(2)) * 10.0 ** (-decimals)


def _parse_2013_dis_cells(text: str) -> list[dict[str, float]]:
    block = _table(text, "Table1")
    rows: list[dict[str, float]] = []
    for match in re.finditer(r"(?m)^\s*(\d+)\s*&(.+?)\\\\\s*$", block):
        fields = [item.strip() for item in match.group(2).split("&")]
        if len(fields) != 6:
            continue
        try:
            values = [float(re.sub(r"[^0-9.+-]", "", item)) for item in fields]
        except ValueError:
            continue
        rows.append({
            "bin": int(match.group(1)),
            "x_low": values[0], "x_high": values[1], "x_mean": values[2],
            "q2_low": values[3], "q2_high": values[4], "q2_mean": values[5],
        })
    if len(rows) != 23 or [item["bin"] for item in rows] != list(range(1, 24)):
        raise SIDISDiagnosticDataError("COMPASS 2013 DIS-cell table changed")
    return rows


def _parse_2013_slopes(text: str, label: str) -> list[list[tuple[float, float]]]:
    block = _table(text, label)
    rows: list[list[tuple[float, float]]] = []
    for match in re.finditer(r"(?m)^\s*(\d+)\s*&(.+?)\\\\\s*$", block):
        fields = [item.strip() for item in match.group(2).split("&")]
        if len(fields) != 8:
            continue
        try:
            parsed = [_parenthesized_value(item) for item in fields]
        except SIDISDiagnosticDataError:
            continue
        if int(match.group(1)) != len(rows) + 1:
            raise SIDISDiagnosticDataError(f"COMPASS 2013 {label} row order changed")
        rows.append(parsed)
    if len(rows) != 23:
        raise SIDISDiagnosticDataError(f"COMPASS 2013 {label} row count changed")
    return rows


def _raw_azimuthal_objects(identifier: str) -> dict[str, str]:
    return {
        "numerator": f"MomentNumerator_{identifier}",
        "denominator": f"DepolarizationDenominator_{identifier}",
        "covariance_positive": f"CovariancePositive_{identifier}",
        "covariance_negative": f"CovarianceNegative_{identifier}",
    }


def normalize_compass_2013() -> dict[str, Any]:
    measurement = "COMPASS_2013_I1236358"
    text, manifest = _source_text(measurement)
    dis_cells = _parse_2013_dis_cells(text)
    slopes = {
        "hplus": _parse_2013_slopes(text, "Table2"),
        "hminus": _parse_2013_slopes(text, "Table3"),
    }
    datasets: dict[str, Any] = {}
    for charge, rows in slopes.items():
        points: list[dict[str, Any]] = []
        slices: list[dict[str, Any]] = []
        for dis_index, (cell, values) in enumerate(zip(dis_cells, rows)):
            flat_bins: list[int] = []
            for z_index, (value, error) in enumerate(values):
                flat_bin = dis_index * 8 + z_index
                flat_bins.append(flat_bin)
                points.append({
                    **cell,
                    "flat_bin": flat_bin,
                    "z_low": Z_EDGES_2013[z_index],
                    "z_high": Z_EDGES_2013[z_index + 1],
                    "z_mean": .5 * (Z_EDGES_2013[z_index] + Z_EDGES_2013[z_index + 1]),
                    "value": value,
                    "stat": error,
                    "systematic": 0.0,
                })
            slices.append({
                "id": f"{charge}_xq2_{dis_index + 1:02d}",
                "edges": list(Z_EDGES_2013),
                "flat_bins": flat_bins,
                "plotted_dimension": "z",
                "rivet_path": (
                    f"/{measurement}/MeanPt2_{charge}_xq2_{dis_index + 1:02d}"
                ),
            })
        datasets[charge] = {
            "id": charge,
            "observable": "low-pT exponential inverse slope",
            "published_target": "D",
            "edges": [float(index) for index in range(185)],
            "rivet_path": f"/{measurement}/MeanPt2_{charge}_cells",
            "raw_objects": {
                "spectrum": f"HadronYield_{charge}_pt2_cells",
            },
            "spectrum_binning": {
                "dis_cells": 23,
                "z_edges": list(Z_EDGES_2013),
                "pt2_edges": list(PT2_EDGES_2013),
                "fit_pt_gev": [.1, .85],
                "flattening": "((dis_cell * 8) + z_bin) * 36 + pt2_bin",
            },
            "points": points,
            "slices": slices,
        }
    return {
        "schema_version": 1,
        "measurement": measurement,
        "observable": "<pT^2> from A exp(-pT^2/<pT^2>)",
        "beam": {"pid": "mu+", "energy_gev": 160.0},
        "selection": {
            "q2_min_gev2": 1.0, "w2_min_gev2": 25.0,
            "y": [.1, .9], "z": [.2, .8], "fit_pt_gev": [.1, .85],
            "hadron": "stable unidentified charged hadron; pion-mass z convention",
        },
        "fit_policy": {
            "model": "A exp(-pT2/slope)",
            "method": "weighted linear least squares in log spectrum",
            "minimum_positive_bins": 3,
            "normalization": (
                "fit target-combined differential hadron yield; the inclusive-DIS "
                "normalization and z width are constant within each fitted cell and "
                "therefore cancel exactly from the slope"
            ),
        },
        "datasets": datasets,
        "target_outputs": {"D": {"P": .5, "N": .5}},
        "systematics": (
            "paper tables quote fit errors only and provide no separate inverse-"
            "slope systematic; a fully correlated multiplicity normalization "
            "component would cancel from an exponential inverse slope"
        ),
        "provenance": {
            "source_manifest": SOURCE_MANIFEST_PATHS[measurement],
            "source_sha256": manifest["source"]["sha256"],
            "paper_tables": [1, 2, 3],
            "hepdata_record": 61432,
        },
    }


def _range(field: str) -> tuple[float, float]:
    numbers = re.findall(r"(?:\d+\.\d*|\.\d+|\d+)", field)
    if len(numbers) < 2:
        raise SIDISDiagnosticDataError(f"Could not parse range {field!r}")
    return float(numbers[0]), float(numbers[1])


def _value_error(field: str) -> tuple[float, float]:
    match = re.search(
        r"([+-]?\s*(?:\d+\.\d*|\.\d+|\d+))\s*\\pm\s*"
        r"([+-]?\s*(?:\d+\.\d*|\.\d+|\d+))",
        field,
    )
    if match is None:
        raise SIDISDiagnosticDataError(f"Could not parse value +/- error {field!r}")
    return float(match.group(1).replace(" ", "")), abs(
        float(match.group(2).replace(" ", ""))
    )


def _parse_2014_1d(text: str, label: str) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for line in _table(text, label).splitlines():
        if r"\pm" not in line:
            continue
        fields = line.split("&")
        if len(fields) < 3:
            continue
        try:
            low, high = _range(fields[0])
            positive, positive_error = _value_error(fields[1])
            negative, negative_error = _value_error(fields[2])
        except SIDISDiagnosticDataError:
            continue
        rows.append({
            "low": low, "high": high,
            "hplus": positive, "hplus_stat": positive_error,
            "hminus": negative, "hminus_stat": negative_error,
        })
    if len(rows) != 24:
        raise SIDISDiagnosticDataError(
            f"COMPASS 2014 {label} expected 24 1D rows, found {len(rows)}"
        )
    return rows


def _parse_2014_3d(text: str, label: str) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    current_z: tuple[float, float] | None = None
    for line in _table(text, label).splitlines():
        if r"\pm" not in line:
            continue
        fields = line.split("&")
        if len(fields) < 4:
            continue
        try:
            if re.search(r"\d", fields[0]):
                current_z = _range(fields[0])
            if current_z is None:
                raise SIDISDiagnosticDataError("3D z range omitted before first row")
            pt_low, pt_high = _range(fields[1])
            positive, positive_error = _value_error(fields[2])
            negative, negative_error = _value_error(fields[3])
        except SIDISDiagnosticDataError:
            continue
        rows.append({
            "z_low": current_z[0], "z_high": current_z[1],
            "pt_low": pt_low, "pt_high": pt_high,
            "hplus": positive, "hplus_stat": positive_error,
            "hminus": negative, "hminus_stat": negative_error,
        })
    if len(rows) != 24:
        raise SIDISDiagnosticDataError(
            f"COMPASS 2014 {label} expected 24 3D rows, found {len(rows)}"
        )
    return rows


def _point(value: float, stat: float, **coordinates: float | int) -> dict[str, Any]:
    return {
        **coordinates,
        "value": value,
        "stat": stat,
        "systematic": 2.0 * stat,
    }


def normalize_compass_2014() -> dict[str, Any]:
    measurement = "COMPASS_2014_I1278730"
    text, manifest = _source_text(measurement)
    datasets: dict[str, Any] = {}
    for harmonic, label in ((1, "tab:values_1d_c"), (2, "tab:values_1d_c2")):
        rows = _parse_2014_1d(text, label)
        offsets = {"x": 0, "z": 7, "pt": 15}
        edges_by_projection = {
            "x": X_EDGES_2014, "z": Z_EDGES_2014, "pt": PT_EDGES_2014,
        }
        sizes = {"x": 7, "z": 8, "pt": 9}
        for projection in ("x", "z", "pt"):
            edges = edges_by_projection[projection]
            selected = rows[offsets[projection]:offsets[projection] + sizes[projection]]
            for charge in ("hplus", "hminus"):
                identifier = f"{charge}_cos{harmonic}_{projection}"
                points = [
                    _point(
                        row[charge], row[f"{charge}_stat"],
                        flat_bin=index,
                        **{
                            f"{projection}_low": row["low"],
                            f"{projection}_high": row["high"],
                            f"{projection}_mean": .5 * (row["low"] + row["high"]),
                        },
                    )
                    for index, row in enumerate(selected)
                ]
                datasets[identifier] = {
                    "id": identifier,
                    "charge": charge,
                    "harmonic": harmonic,
                    "projection": projection,
                    "published_target": "D",
                    "edges": list(edges),
                    "rivet_path": f"/{measurement}/AUU_{identifier}",
                    "raw_objects": _raw_azimuthal_objects(identifier),
                    "points": points,
                    "slices": [],
                }

    for harmonic, prefix in ((1, "tab:values_3d_c"), (2, "tab:values_3d_c2")):
        rows_by_x = [
            _parse_2014_3d(text, f"{prefix}_{index}") for index in range(1, 5)
        ]
        for charge in ("hplus", "hminus"):
            identifier = f"{charge}_cos{harmonic}_x_z_pt"
            points: list[dict[str, Any]] = []
            slices: list[dict[str, Any]] = []
            for x_index, rows in enumerate(rows_by_x):
                for row_index, row in enumerate(rows):
                    z_index, pt_index = divmod(row_index, 4)
                    expected_z = (Z3_EDGES_2014[z_index], Z3_EDGES_2014[z_index + 1])
                    expected_pt = (PT3_EDGES_2014[pt_index], PT3_EDGES_2014[pt_index + 1])
                    if not all(math.isclose(left, right) for left, right in zip(
                        (row["z_low"], row["z_high"], row["pt_low"], row["pt_high"]),
                        (*expected_z, *expected_pt),
                    )):
                        raise SIDISDiagnosticDataError(
                            f"COMPASS 2014 {identifier} 3D bin order changed"
                        )
                    flat_bin = (x_index * 6 + z_index) * 4 + pt_index
                    points.append(_point(
                        row[charge], row[f"{charge}_stat"], flat_bin=flat_bin,
                        x_low=X3_EDGES_2014[x_index],
                        x_high=X3_EDGES_2014[x_index + 1],
                        x_mean=.5 * (X3_EDGES_2014[x_index] + X3_EDGES_2014[x_index + 1]),
                        z_low=row["z_low"], z_high=row["z_high"],
                        z_mean=.5 * (row["z_low"] + row["z_high"]),
                        pt_low=row["pt_low"], pt_high=row["pt_high"],
                        pt_mean=.5 * (row["pt_low"] + row["pt_high"]),
                    ))
                for z_index in range(6):
                    slices.append({
                        "id": f"{identifier}_x{x_index + 1:02d}_z{z_index + 1:02d}",
                        "edges": list(PT3_EDGES_2014),
                        "flat_bins": [
                            (x_index * 6 + z_index) * 4 + pt_index
                            for pt_index in range(4)
                        ],
                        "plotted_dimension": "pt",
                        "rivet_path": (
                            f"/{measurement}/AUU_{identifier}_x{x_index + 1:02d}"
                            f"_z{z_index + 1:02d}"
                        ),
                    })
            points.sort(key=lambda item: int(item["flat_bin"]))
            datasets[identifier] = {
                "id": identifier,
                "charge": charge,
                "harmonic": harmonic,
                "projection": "x_z_pt",
                "published_target": "D",
                "edges": [float(index) for index in range(97)],
                "rivet_path": f"/{measurement}/AUU_{identifier}_cells",
                "raw_objects": _raw_azimuthal_objects(identifier),
                "points": points,
                "slices": slices,
            }

    if len(datasets) != 16 or sum(len(item["points"]) for item in datasets.values()) != 480:
        raise SIDISDiagnosticDataError("COMPASS 2014 normalized inventory changed")
    return {
        "schema_version": 1,
        "measurement": measurement,
        "observable": "A_UU^cos(phi) and A_UU^cos(2phi)",
        "beam": {"pid": "mu+", "energy_gev": 160.0},
        "selection": {
            "q2_min_gev2": 1.0, "w2_min_gev2": 25.0,
            "x": [.003, .13], "y": [.2, .9],
            "virtual_photon_angle_max_rad": .060,
            "z": [.2, .85], "pt_gev": [.1, 1.0],
            "hadron": "stable unidentified charged hadron; pion-mass z convention",
        },
        "estimator": {
            "numerator": "sum 2 cos(n phi_h)",
            "denominator": "sum epsilon_n(y)",
            "epsilon_1": "2(2-y)sqrt(1-y)/(1+(1-y)^2)",
            "epsilon_2": "2(1-y)/(1+(1-y)^2)",
            "same_event_covariance": True,
        },
        "datasets": datasets,
        "target_outputs": {"D": {"P": .5, "N": .5}},
        "systematics": (
            "point-to-point systematic uncertainty is twice the tabulated "
            "statistical uncertainty, as stated for both integrated and 3D results"
        ),
        "excluded_observable": (
            "A_LU^sin(phi_h) is not evaluated in the 00 samples because it "
            "requires the polarized muon beam"
        ),
        "provenance": {
            "source_manifest": SOURCE_MANIFEST_PATHS[measurement],
            "source_sha256": manifest["source"]["sha256"],
            "paper_tables": list(range(2, 13)),
            "hepdata_record": 64754,
        },
    }


def normalized_from_raw(measurement: str) -> dict[str, Any]:
    if measurement == "COMPASS_2013_I1236358":
        return normalize_compass_2013()
    if measurement == "COMPASS_2014_I1278730":
        return normalize_compass_2014()
    raise SIDISDiagnosticDataError(f"Unknown SIDIS diagnostic {measurement}")


def validate_vendored(
    measurement: str, *, full_covariance: bool = False
) -> dict[str, Any]:
    if full_covariance:
        raise SIDISDiagnosticDataError(
            "No full-covariance mode is defined for the SIDIS diagnostics"
        )
    generated = normalized_from_raw(measurement)
    snapshot = _load_json(ROOT / REFERENCE_PATHS[measurement])
    _assert_equal(generated, snapshot)
    return snapshot


def _import_yoda() -> Any:
    try:
        import yoda  # type: ignore
    except (ImportError, OSError) as exc:
        raise SIDISDiagnosticDataError(
            "YODA Python bindings are unavailable; load the HerwigPol environment"
        ) from exc
    return yoda


def _estimate(
    yoda: Any, edges: Sequence[float], points: Sequence[Mapping[str, Any]], path: str
) -> Any:
    result = yoda.BinnedEstimate1D(list(map(float, edges)), "/REF" + path)
    result.setAnnotation("IsRef", 1)
    for index, point in enumerate(points, start=1):
        target = result.bin(index)
        target.setVal(float(point["value"]))
        target.setErr(-float(point["stat"]), float(point["stat"]), "stat")
        systematic = float(point["systematic"])
        target.setErr(-systematic, systematic, "syst")
    return result


def write_reference_yoda(
    measurement: str, snapshot: Mapping[str, Any] | None = None
) -> Path:
    snapshot = dict(snapshot or validate_vendored(measurement))
    yoda = _import_yoda()
    objects: list[Any] = []
    for dataset in snapshot["datasets"].values():
        objects.append(_estimate(
            yoda, dataset["edges"], dataset["points"], dataset["rivet_path"]
        ))
        for slice_spec in dataset.get("slices", []):
            points = [
                dataset["points"][int(index)] for index in slice_spec["flat_bins"]
            ]
            objects.append(_estimate(
                yoda, slice_spec["edges"], points, slice_spec["rivet_path"]
            ))
    destination = ROOT / REFERENCE_YODA_PATHS[measurement]
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix="sidis-diagnostic-ref-", dir=destination.parent
    ) as directory:
        plain = Path(directory) / destination.with_suffix("").name
        yoda.write(objects, str(plain))
        compressed = Path(directory) / destination.name
        with plain.open("rb") as source, compressed.open("wb") as target:
            with gzip.GzipFile(
                filename="", mode="wb", fileobj=target, mtime=0
            ) as stream:
                shutil.copyfileobj(source, stream)
        os.replace(compressed, destination)
    return destination


def fetch_and_validate(
    measurement: str,
    cache_directory: Path | None = None,
    source_file: Path | None = None,
) -> list[Path]:
    del cache_directory
    if source_file is not None:
        raise SIDISDiagnosticDataError(
            "The diagnostic paper sources are already checksum-pinned"
        )
    snapshot = validate_vendored(measurement)
    yoda_path = write_reference_yoda(measurement, snapshot)
    manifest = _load_json(ROOT / SOURCE_MANIFEST_PATHS[measurement])
    return [
        ROOT / str(manifest["source"]["path"]),
        ROOT / REFERENCE_PATHS[measurement],
        yoda_path,
    ]


def refresh(measurement: str) -> Path:
    snapshot = normalized_from_raw(measurement)
    return _atomic_json(ROOT / REFERENCE_PATHS[measurement], snapshot)
