#!/usr/bin/env python3
"""Construct HERMES Born A_parallel projections conditional on fitted UU weights.

This is a derived cell estimator, not an additional HERMES measurement.  It
retains the unfolded Born values and their full statistical covariance.  The
fit supplies unpolarized cross sections; it supplies no spin information.
Splitting a published cell assumes its A_parallel is constant within that cell.
"""
from __future__ import annotations

import argparse
import csv
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Callable, Mapping


ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "data/experimental/HERMES_2007_I726689/born-apar-reference.json"
MEASUREMENT = "HERMES_2007_I726689"
Q2_EDGES = (1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 20.0)
Integrator = Callable[[str, float, float, float, float], float]


class ProjectionError(RuntimeError):
    """Inputs or derived projection do not satisfy the declared cell contract."""


def _finite(value: Any, label: str) -> float:
    if isinstance(value, bool):
        raise ProjectionError(f"Invalid {label}")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ProjectionError(f"Invalid {label}") from exc
    if not math.isfinite(number):
        raise ProjectionError(f"Non-finite {label}")
    return number


def _source_cells(reference: Mapping[str, Any], target: str) -> list[dict[str, Any]]:
    cells = []
    for dataset in reference["datasets"]:
        if dataset["target"] != target:
            continue
        for point in dataset["points"]:
            cells.append({
                **point, "selection": dataset["selection"],
                "x_low": dataset["x_low"], "x_high": dataset["x_high"],
            })
    cells.sort(key=lambda cell: cell["global_bin"])
    if [cell["global_bin"] for cell in cells] != list(range(1, 46)):
        raise ProjectionError(f"Target {target} must have exactly the 45 published cells")
    for cell in cells:
        for key in ("x_low", "x_high", "q2_low", "q2_high", "value", "stat",
                    "systematic_combined", "x_mean", "q2_mean"):
            cell[key] = _finite(cell[key], f"{target} cell {cell['global_bin']} {key}")
        if not (cell["x_low"] < cell["x_high"] and cell["q2_low"] < cell["q2_high"]):
            raise ProjectionError("Cell boundaries must be ordered")
        if cell["stat"] < 0 or cell["systematic_combined"] < 0:
            raise ProjectionError("Cell uncertainties must be nonnegative")
        if bool(cell["supported_prediction"]) != (cell["q2_low"] >= 1.0):
            raise ProjectionError("A cell support mask disagrees with the Q2=1 boundary")
    return cells


def _covariance(reference: Mapping[str, Any], target: str) -> list[list[float]]:
    block = reference["statistical_covariance"][target]
    if block["global_bin_order"] != list(range(1, 46)):
        raise ProjectionError("Statistical covariance global-bin order changed")
    matrix = block["matrix"]
    if len(matrix) != 45 or any(len(row) != 45 for row in matrix):
        raise ProjectionError("Expected a 45 by 45 statistical covariance")
    matrix = [[_finite(value, "statistical covariance") for value in row] for row in matrix]
    for i in range(45):
        if matrix[i][i] < 0:
            raise ProjectionError("Negative diagonal statistical variance")
        for j in range(i):
            if not math.isclose(matrix[i][j], matrix[j][i], rel_tol=1.e-12, abs_tol=1.e-16):
                raise ProjectionError("Statistical covariance is not symmetric")
    return matrix


def _mc_cells(summary: Mapping[str, Any], cells: list[dict[str, Any]]) -> list[dict[str, float] | None]:
    lookup = {}
    for row in summary["bins"]:
        if not str(row.get("selection", "")).startswith("Born"):
            continue
        identity = (row["selection"], row["bin"])
        if identity in lookup:
            raise ProjectionError(f"Repeated theory cell {identity}")
        lookup[identity] = row
    result = []
    for cell in cells:
        identity = (cell["selection"], cell["bin"])
        row = lookup.get(identity)
        if row is None:
            raise ProjectionError(f"Theory summary is missing {identity}")
        if bool(row.get("supported")) != bool(cell["supported_prediction"]):
            raise ProjectionError(f"Theory/reference support mismatch for {identity}")
        if not all(math.isclose(_finite(row[key], key), cell[source], abs_tol=1.e-12, rel_tol=0)
                   for key, source in (("bin_low", "q2_low"), ("bin_high", "q2_high"))):
            raise ProjectionError(f"Theory/reference Q2 boundaries differ for {identity}")
        if not cell["supported_prediction"]:
            result.append(None)
            continue
        value = _finite(row["a_parallel"], f"{identity} A_parallel")
        error = _finite(row["a_parallel_stat"], f"{identity} A_parallel error")
        if error < 0:
            raise ProjectionError("Negative theory statistical error")
        result.append({"value": value, "stat": error})
    return result


def projected_covariance(weights: list[list[float]], covariance: list[list[float]]) -> list[list[float]]:
    """Propagate all source-cell correlations, including reused split cells."""
    count = len(covariance)
    if any(len(row) != count for row in covariance) or any(len(row) != count for row in weights):
        raise ProjectionError("Projection/covariance dimensions disagree")
    product = [[sum(row[j] * covariance[j][k] for j in range(count))
                for k in range(count)] for row in weights]
    return [[sum(left[k] * right[k] for k in range(count)) for right in weights]
            for left in product]


def _projection_rows(cells: list[dict[str, Any]], target: str) -> list[dict[str, Any]]:
    slices = sorted(set((cell["x_low"], cell["x_high"]) for cell in cells
                        if cell["supported_prediction"]))
    if len(slices) != 15:
        raise ProjectionError("Expected 15 supported physical x slices")
    rows = []
    for minimum in (1.0, 4.0):
        suffix = f"Q2GT{int(minimum)}"
        for index, (low, high) in enumerate(slices, start=1):
            rows.append({
                "id": f"{target}_X_{suffix}_{index:02d}", "target": target,
                "projection": "x", "selection": suffix, "bin": index,
                "bin_low": low, "bin_high": high, "x_low": low, "x_high": high,
                "q2_low": minimum, "q2_high": 20.0,
            })
        edges = [value for value in Q2_EDGES if value >= minimum]
        for index, (low, high) in enumerate(zip(edges, edges[1:]), start=1):
            rows.append({
                "id": f"{target}_Q2_{suffix}_{index:02d}", "target": target,
                "projection": "q2", "selection": suffix, "bin": index,
                "bin_low": low, "bin_high": high,
                "x_low": slices[0][0], "x_high": slices[-1][1],
                "q2_low": low, "q2_high": high,
            })
    return rows


def build_projection_snapshot(
    reference: Mapping[str, Any], summary: Mapping[str, Any], integrator: Integrator,
    *, fit_provenance: Mapping[str, Any] | None = None,
    controls: Mapping[str, Integrator] | None = None,
) -> dict[str, Any]:
    """Apply independent fitted-UU weights to data and the same MC cell ratios.

    integrator integrates the cross section, not its density, inside the
    rectangular cell intersection and the shared fiducial selection.  Source
    asymmetries therefore require no bin-width factors.  The production
    sigma_uu_pb density columns are deliberately unused by this estimator.
    """
    if reference.get("measurement") != MEASUREMENT or summary.get("measurement") != MEASUREMENT:
        raise ProjectionError("Reference and theory must identify HERMES_2007_I726689")

    @lru_cache(maxsize=None)
    def integral(target: str, xlow: float, xhigh: float, qlow: float, qhigh: float) -> float:
        value = _finite(integrator(target, xlow, xhigh, qlow, qhigh), "fitted UU integral")
        if value < 0:
            raise ProjectionError("Fitted UU integrals must be nonnegative")
        return value

    targets = {}
    for target in ("P", "D"):
        cells = _source_cells(reference, target)
        covariance = _covariance(reference, target)
        mc = _mc_cells(summary, cells)
        rows = _projection_rows(cells, target)
        matrix = []
        for row in rows:
            yields = [0.0] * 45
            split = []
            intersections = []
            for index, cell in enumerate(cells):
                if not cell["supported_prediction"]:
                    continue
                xlow, xhigh = max(row["x_low"], cell["x_low"]), min(row["x_high"], cell["x_high"])
                qlow, qhigh = max(row["q2_low"], cell["q2_low"]), min(row["q2_high"], cell["q2_high"])
                if xlow >= xhigh or qlow >= qhigh:
                    continue
                yields[index] = integral(target, xlow, xhigh, qlow, qhigh)
                if yields[index] == 0:
                    continue
                intersections.append({"global_bin": cell["global_bin"],
                                      "x_low": xlow, "x_high": xhigh,
                                      "q2_low": qlow, "q2_high": qhigh,
                                      "fitted_uu_pb": yields[index]})
                if (xlow, xhigh, qlow, qhigh) != (cell["x_low"], cell["x_high"], cell["q2_low"], cell["q2_high"]):
                    split.append(cell["global_bin"])
            total = sum(yields)
            weights = [value / total for value in yields] if total else [0.0] * 45
            matrix.append(weights)
            supported = total > 0
            center = (sum(weight * cell["x_mean"] for weight, cell in zip(weights, cells))
                      if row["projection"] == "x" and supported
                      else math.sqrt(row["bin_low"] * row["bin_high"]))
            if not row["bin_low"] <= center <= row["bin_high"]:
                raise ProjectionError("Derived marker lies outside its projection bin")
            row.update(
                supported=supported, marker=center,
                marker_convention=("fitted-weight average of published full-cell mean x"
                                   if row["projection"] == "x" else "geometric target-bin center"),
                fitted_uu_pb=total,
                source_global_bins=[item["global_bin"] for item in intersections],
                split_global_bins=split, cell_constant_assumption=bool(split),
                intersections=intersections,
                a_parallel=(sum(weight * cell["value"] for weight, cell in zip(weights, cells))
                            if supported else None),
                systematic_upper_bound=(sum(abs(weight) * cell["systematic_combined"]
                                            for weight, cell in zip(weights, cells)) if supported else None),
                systematic_diagonal_diagnostic=(math.sqrt(sum(weight * weight * cell["systematic_combined"] ** 2
                                                              for weight, cell in zip(weights, cells)))
                                                if supported else None),
                mc_a_parallel=(sum(weight * prediction["value"] for weight, prediction in zip(weights, mc)
                                   if prediction is not None) if supported else None),
                mc_stat_independent_cell_approx=(math.sqrt(sum(weight * weight * prediction["stat"] ** 2
                                                                for weight, prediction in zip(weights, mc)
                                                                if prediction is not None)) if supported else None),
            )
        projected = projected_covariance(matrix, covariance)
        for index, row in enumerate(rows):
            variance = projected[index][index]
            if variance < -1.e-14:
                raise ProjectionError("Projected statistical variance is negative")
            row["stat"] = math.sqrt(max(0., variance)) if row["supported"] else None
            row["stat_plus_systematic_upper_bound"] = (
                math.hypot(row["stat"], row["systematic_upper_bound"]) if row["supported"] else None
            )
        # Null covariance entries preserve masks rather than suggesting measured zero errors.
        projected = [[value if rows[i]["supported"] and rows[j]["supported"] else None
                      for j, value in enumerate(values)] for i, values in enumerate(projected)]
        target_datasets = [dataset for dataset in reference["datasets"] if dataset["target"] == target]
        targets[target] = {
            "source_global_bin_order": list(range(1, 46)),
            "source_cells": cells,
            "projection_order": [row["id"] for row in rows],
            "bins": rows, "weights": matrix,
            "statistical_covariance": projected,
            "normalization_uncertainty": target_datasets[0]["provenance"]["normalization_uncertainty"],
            "normalization_policy": "Already included in published cell systematic errors; neither added nor subtracted",
            "covariance_scope": "Full supplied source statistical covariance; conditional on fixed fitted UU weights",
        }
    snapshot = {
        "schema_version": 1, "measurement": MEASUREMENT,
        "observable": "Derived Born A_parallel cell projections with independent fitted UU weights",
        "campaign_tag": summary.get("tag"), "targets": targets,
        "fit": dict(fit_provenance or {}),
        "assumptions": [
            "Unfolded Born values are used unchanged; no D, polarization, dilution, tensor or nuclear correction is reapplied to data",
            "The fit supplies only unpolarized cross-section weights; no Herwig spin prediction enters experimental central values",
            "Whole-cell x projections at Q2>=1 require no within-cell asymmetry interpolation",
            "Split cells in Q2 projections and Q2>=4 controls assume A_parallel constant within each published cell",
            "The same projection matrix is applied to the data and production Born-cell MC asymmetries; these MC values differ from event-level projections",
            "Published 45 by 45 statistical covariance is propagated; cross-target statistical covariance is unavailable",
            "Systematic upper bounds follow arbitrary allowed correlations with published diagonal errors; the diagonal result is a diagnostic assumption",
            "Published normalization uncertainty is already included in systematic errors and is not added or subtracted again",
            "Fit parameter covariance and within-cell spin-shape uncertainty are not included in conditional experimental bars",
            "MC statistical errors use an independent-cell approximation; full MC cell covariance is unavailable",
            "Unsupported Q2<1 source cells retain zero projection weights and are never extrapolated",
        ],
        "uncertainties": {
            "experimental_statistical": "W C_stat W^T, including all published off-diagonal correlations",
            "experimental_systematic_upper_bound": "sum_i abs(W_ri) s_i; not an inferred covariance or confidence level",
            "experimental_systematic_diagonal_diagnostic": "sqrt(sum_i W_ri^2 s_i^2); requires independent source systematics",
            "stat_plus_systematic_upper_bound": "Quadrature combination of statistical error and systematic upper bound; excludes reconstruction model errors",
            "fit_weight_uncertainty": "Not included in conditional errors; requires separate fit covariance/model assessment",
            "mc_statistical": "Independent source-cell approximation; no full MC covariance claim",
        },
    }
    if controls:
        snapshot["controls"] = {}
        for name, control_integrator in controls.items():
            control = build_projection_snapshot(reference, summary, control_integrator)
            snapshot["controls"][name] = {
                "scope": "Reconstruction sensitivity; not a one-sigma error or an independent experimental measurement",
                "targets": control["targets"],
            }
            for target in ("P", "D"):
                central_rows = snapshot["targets"][target]["bins"]
                varied_rows = control["targets"][target]["bins"]
                for central, varied in zip(central_rows, varied_rows):
                    valid = central["supported"] and varied["supported"]
                    variation = {
                        "supported": varied["supported"],
                        "data_shift": varied["a_parallel"] - central["a_parallel"] if valid else None,
                        "mc_shift": varied["mc_a_parallel"] - central["mc_a_parallel"] if valid else None,
                        "data_minus_mc_shift": ((varied["a_parallel"] - varied["mc_a_parallel"])
                                                - (central["a_parallel"] - central["mc_a_parallel"])) if valid else None,
                        "fitted_uu_pb": varied["fitted_uu_pb"],
                    }
                    central.setdefault("weight_controls", {})[name] = variation
    return snapshot


def write_outputs(snapshot: Mapping[str, Any], output: Path) -> None:
    if output.exists() and not output.is_dir():
        raise ProjectionError(f"Output path is not a directory: {output}")
    if output.exists() and any(output.iterdir()):
        raise ProjectionError(f"Output directory is not empty: {output}; choose a new reconstruction directory")
    output.mkdir(parents=True, exist_ok=True)
    (output / "integrated.json").write_text(json.dumps(snapshot, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    columns = ["id", "target", "projection", "selection", "bin", "bin_low", "bin_high", "marker", "supported",
               "a_parallel", "stat", "systematic_upper_bound", "systematic_diagonal_diagnostic",
               "stat_plus_systematic_upper_bound", "mc_a_parallel", "mc_stat_independent_cell_approx",
               "fitted_uu_pb", "cell_constant_assumption", "source_global_bins", "split_global_bins",
               "fit_domain_extrapolated_uu_fraction"]
    control_names = list(snapshot.get("controls", {}))
    control_fields = ("data_shift", "mc_shift", "data_minus_mc_shift")
    columns.extend(f"{name}_{field}" for name in control_names for field in control_fields)
    with (output / "integrated.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for target in ("P", "D"):
            for row in snapshot["targets"][target]["bins"]:
                values = {key: row.get(key) for key in columns}
                for key in ("source_global_bins", "split_global_bins"):
                    values[key] = " ".join(str(value) for value in values[key])
                for name in control_names:
                    for field in control_fields:
                        values[f"{name}_{field}"] = row["weight_controls"][name][field]
                writer.writerow(values)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, default=REFERENCE)
    parser.add_argument("--summary", type=Path, required=True, help="Existing production postprocess/summary.json")
    parser.add_argument("--output", type=Path, required=True, help="New, empty reconstruction directory")
    parser.add_argument("--model", choices=("GD11", "ALLM97_HYBRID"), default="GD11")
    parser.add_argument("--quadrature-order", type=int, default=32)
    parser.add_argument("--no-plots", action="store_true", help="Write numerical reconstruction only")
    parser.add_argument("--no-controls", action="store_true", help="Skip fit/R/acceptance/convergence sensitivity scans")
    args = parser.parse_args()
    from hermes_unpolarized_fit import DEFAULT_CUTS, MODEL_METADATA, SNAPSHOT_PATH, integrate_cross_section
    if args.reference.resolve() == REFERENCE.resolve():
        from hermes_born_reference_data import validate_vendored
        reference = validate_vendored(ROOT)
    else:
        reference = json.loads(args.reference.read_text(encoding="utf-8"))
    summary = json.loads(args.summary.read_text(encoding="utf-8"))
    def make_integrator(**changes: Any) -> Integrator:
        options = {"model": args.model, "quadrature_order": args.quadrature_order}
        options.update(changes)
        def integrator(target: str, xlow: float, xhigh: float, qlow: float, qhigh: float) -> float:
            return integrate_cross_section(target, xlow, xhigh, qlow, qhigh, **options)
        return integrator
    control_integrators = None if args.no_controls else {
        "ALLM97_HYBRID": make_integrator(model="ALLM97_HYBRID"),
        "R1990": make_integrator(r_model="R1990"),
        "W2GT4": make_integrator(w2_min=4.0),
        "quadrature_double": make_integrator(quadrature_order=2 * args.quadrature_order),
    }
    fit = {"model": args.model, "quadrature_order": args.quadrature_order,
           "cuts": DEFAULT_CUTS, "metadata": MODEL_METADATA,
           "controls": {
               "ALLM97_HYBRID": "Alternative unpolarized F2 shape; hybrid deuteron uses the GD11 D/P ratio",
               "R1990": "Alternative R parameterization in the fitted unpolarized cross section",
               "W2GT4": "Stricter W2 domain; changes reconstructed acceptance and is not a same-fiducial uncertainty",
               "quadrature_double": "Numerical convergence check at twice the central quadrature order",
           } if control_integrators else {},
    }
    snapshot = build_projection_snapshot(reference, summary, make_integrator(), fit_provenance=fit,
                                         controls=control_integrators)
    if control_integrators:
        for target in ("P", "D"):
            for row in snapshot["targets"][target]["bins"]:
                central, strict = row["fitted_uu_pb"], row["weight_controls"]["W2GT4"]["fitted_uu_pb"]
                fraction = 1.0 - strict / central if central else None
                if fraction is not None and not -1.e-12 <= fraction <= 1.0 + 1.e-12:
                    raise ProjectionError("Strict W2-domain integral is not a subset of the central acceptance")
                row["fraction_of_fitted_uu_below_W2_4"] = max(0., min(1., fraction)) if fraction is not None else None
                row["fit_domain_extrapolated_uu_fraction"] = row["fraction_of_fitted_uu_below_W2_4"]
    snapshot["input_files"] = {
        "reference": {"path": str(args.reference.resolve()), "sha256": hashlib.sha256(args.reference.read_bytes()).hexdigest()},
        "summary": {"path": str(args.summary.resolve()), "sha256": hashlib.sha256(args.summary.read_bytes()).hexdigest()},
        "unpolarized_fit": {"path": str(SNAPSHOT_PATH.resolve()), "sha256": hashlib.sha256(SNAPSHOT_PATH.read_bytes()).hexdigest()},
    }
    write_outputs(snapshot, args.output)
    if not args.no_plots:
        from hermes_apar_integrated_plots import render
        render(snapshot, args.output)
    print(f"Wrote 84 derived P/D projection bins and full statistical covariance to {args.output}")


if __name__ == "__main__":
    try:
        main()
    except (ProjectionError, ValueError) as exc:
        raise SystemExit(f"error: {exc}") from exc
