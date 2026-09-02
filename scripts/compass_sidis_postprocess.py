#!/usr/bin/env python3
"""Normalized-bin estimators for the COMPASS SIDIS measurements.

This module contains no Rivet or Herwig execution logic.  It operates on the
normalized ``BinSeries`` objects produced after shard combination and after
the POSNLO and NEGNLO samples have been added.  Keeping the arithmetic here
makes the target, helicity, covariance, and masking contracts directly unit
testable.
"""

from __future__ import annotations

import math
from typing import Any, Mapping, Sequence

import run_experimental_campaign as experimental


HELICITIES = ("PP", "PM", "MP", "MM")
LONGITUDINAL_SIGN = {"PP": 1.0, "PM": -1.0, "MP": -1.0, "MM": 1.0}
DEUTERON_POLARIZATION_FACTOR = 0.925


def _check_series(
    collections: Sequence[Mapping[str, experimental.BinSeries]],
    labels: set[str],
) -> experimental.BinSeries:
    if not collections or any(set(collection) != labels for collection in collections):
        raise experimental.CampaignError(
            "COMPASS SIDIS component labels do not match the estimator"
        )
    first = next(iter(collections[0].values()))
    for collection in collections:
        for series in collection.values():
            if not experimental._same_edges(first.edges, series.edges):
                raise experimental.CampaignError(
                    "COMPASS SIDIS components have inconsistent bin edges"
                )
    return first


def a1_target_combination(
    ordinary: Mapping[str, experimental.BinSeries],
    inverse_depolarization: Mapping[str, experimental.BinSeries],
    covariance_proxy: Mapping[str, experimental.BinSeries],
    target_weights: Mapping[str, float],
    longitudinal_target_scale: float,
) -> dict[str, Any]:
    """Form ``A1`` for an explicit target-component combination.

    The ordinary yield is the denominator and the longitudinal difference is
    formed from event-wise inverse-depolarization yields.  Each physical
    helicity receives the averaging factor one quarter.  Target weights and
    the longitudinal scale are descriptor data, making proton and isoscalar
    analyses use exactly the same tested arithmetic.
    """
    if not target_weights or not math.isfinite(float(longitudinal_target_scale)):
        raise experimental.CampaignError("SIDIS A1 target weights/scale are invalid")
    if any(not math.isfinite(float(value)) for value in target_weights.values()):
        raise experimental.CampaignError("SIDIS A1 target weights must be finite")
    labels = {
        f"{target}:{helicity}"
        for target in target_weights
        for helicity in HELICITIES
    }
    first = _check_series(
        (ordinary, inverse_depolarization, covariance_proxy), labels
    )
    denominator_coefficients = {
        f"{target}:{helicity}": .25 * float(weight)
        for target, weight in target_weights.items()
        for helicity in HELICITIES
    }
    numerator_coefficients = {
        f"{target}:{helicity}": .25 * float(weight)
        * float(longitudinal_target_scale) * LONGITUDINAL_SIGN[helicity]
        for target, weight in target_weights.items()
        for helicity in HELICITIES
    }
    values: list[float | None] = []
    errors: list[float | None] = []
    numerator_values: list[float] = []
    numerator_variances: list[float] = []
    denominator_values: list[float] = []
    denominator_variances: list[float] = []
    covariances: list[float] = []
    for index in range(len(first.values)):
        numerator = sum(
            coefficient * inverse_depolarization[label].values[index]
            for label, coefficient in numerator_coefficients.items()
        )
        numerator_variance = sum(
            coefficient**2 * inverse_depolarization[label].variances[index]
            for label, coefficient in numerator_coefficients.items()
        )
        denominator = sum(
            coefficient * ordinary[label].values[index]
            for label, coefficient in denominator_coefficients.items()
        )
        denominator_variance = sum(
            coefficient**2 * ordinary[label].variances[index]
            for label, coefficient in denominator_coefficients.items()
        )
        covariance = sum(
            numerator_coefficients[label]
            * denominator_coefficients[label]
            * covariance_proxy[label].variances[index]
            for label in labels
        )
        if denominator <= 0.0:
            value, error = None, None
        else:
            value, error = experimental.ratio_with_covariance(
                numerator,
                numerator_variance,
                denominator,
                denominator_variance,
                covariance,
            )
        values.append(value)
        errors.append(error)
        numerator_values.append(numerator)
        numerator_variances.append(numerator_variance)
        denominator_values.append(denominator)
        denominator_variances.append(denominator_variance)
        covariances.append(covariance)

    return {
        "edges": list(first.edges),
        "values": values,
        "errors": errors,
        "sigma_ll_over_d": experimental.BinSeries(
            list(first.edges), numerator_values, numerator_variances
        ),
        "sigma_uu": experimental.BinSeries(
            list(first.edges), denominator_values, denominator_variances
        ),
        "numerator_denominator_covariance": covariances,
    }


def a1_deuteron(
    ordinary: Mapping[str, experimental.BinSeries],
    inverse_depolarization: Mapping[str, experimental.BinSeries],
    covariance_proxy: Mapping[str, experimental.BinSeries],
) -> dict[str, Any]:
    """Backward-compatible schema-5 COMPASS 2009 estimator."""
    return a1_target_combination(
        ordinary, inverse_depolarization, covariance_proxy,
        {"P": .5, "N": .5}, DEUTERON_POLARIZATION_FACTOR,
    )


def multiplicity_target_combination(
    numerators: Mapping[str, experimental.BinSeries],
    denominators: Mapping[str, experimental.BinSeries],
    covariance_proxies: Mapping[str, experimental.BinSeries],
    target_weights: Mapping[str, float],
    density_widths: Sequence[float],
) -> dict[str, Any]:
    """Form a density after an explicit target-component sum.

    Non-finite or non-positive inclusive-DIS denominators are represented by
    masked bins.  The same-event numerator--denominator covariance is retained
    before dividing the ratio and its Monte Carlo error by the published cell
    width in ``z``.
    """

    labels = set(target_weights)
    first = _check_series((numerators, denominators, covariance_proxies), labels)
    if not labels or len(density_widths) != len(first.values) or any(
        not math.isfinite(float(width)) or float(width) <= 0.0
        for width in density_widths
    ):
        raise experimental.CampaignError(
            "SIDIS density widths do not match the flattened cells"
        )

    values: list[float | None] = []
    errors: list[float | None] = []
    numerator_values: list[float] = []
    denominator_values: list[float] = []
    covariances: list[float] = []
    for index, width in enumerate(density_widths):
        numerator = sum(float(weight) * numerators[target].values[index]
                        for target, weight in target_weights.items())
        numerator_variance = sum(float(weight)**2 * numerators[target].variances[index]
                                 for target, weight in target_weights.items())
        denominator = sum(float(weight) * denominators[target].values[index]
                          for target, weight in target_weights.items())
        denominator_variance = sum(float(weight)**2 * denominators[target].variances[index]
                                   for target, weight in target_weights.items())
        covariance = sum(float(weight)**2 * covariance_proxies[target].variances[index]
                         for target, weight in target_weights.items())
        if denominator <= 0.0:
            value, error = None, None
        else:
            value, error = experimental.ratio_with_covariance(
                numerator,
                numerator_variance,
                denominator,
                denominator_variance,
                covariance,
            )
            if value is not None and error is not None:
                value /= float(width)
                error /= float(width)
        values.append(value)
        errors.append(error)
        numerator_values.append(numerator)
        denominator_values.append(denominator)
        covariances.append(covariance)

    return {
        "edges": list(first.edges),
        "values": values,
        "errors": errors,
        "isoscalar_numerator": numerator_values,
        "isoscalar_denominator": denominator_values,
        "numerator_denominator_covariance": covariances,
    }


def multiplicity_isoscalar(
    numerators: Mapping[str, experimental.BinSeries],
    denominators: Mapping[str, experimental.BinSeries],
    covariance_proxies: Mapping[str, experimental.BinSeries],
    z_widths: Sequence[float],
) -> dict[str, Any]:
    """Backward-compatible schema-5 COMPASS 2017 estimator."""
    return multiplicity_target_combination(
        numerators, denominators, covariance_proxies,
        {"P": .5, "N": .5}, z_widths,
    )


def azimuthal_target_combination(
    numerators: Mapping[str, experimental.BinSeries],
    denominators: Mapping[str, experimental.BinSeries],
    covariance_positive: Mapping[str, experimental.BinSeries],
    covariance_negative: Mapping[str, experimental.BinSeries],
    target_weights: Mapping[str, float],
) -> dict[str, Any]:
    """Form a cosine amplitude after target combination.

    Rivet stores a signed same-event numerator--denominator covariance as the
    difference of the ``sumW2`` arrays of two non-negative proxy histograms.
    Target components are combined before the ratio, exactly as for the
    multiplicity estimators.
    """

    labels = set(target_weights)
    first = _check_series(
        (numerators, denominators, covariance_positive, covariance_negative),
        labels,
    )
    values: list[float | None] = []
    errors: list[float | None] = []
    covariances: list[float] = []
    for index in range(len(first.values)):
        numerator = sum(
            float(weight) * numerators[target].values[index]
            for target, weight in target_weights.items()
        )
        numerator_variance = sum(
            float(weight) ** 2 * numerators[target].variances[index]
            for target, weight in target_weights.items()
        )
        denominator = sum(
            float(weight) * denominators[target].values[index]
            for target, weight in target_weights.items()
        )
        denominator_variance = sum(
            float(weight) ** 2 * denominators[target].variances[index]
            for target, weight in target_weights.items()
        )
        covariance = sum(
            float(weight) ** 2
            * (
                covariance_positive[target].variances[index]
                - covariance_negative[target].variances[index]
            )
            for target, weight in target_weights.items()
        )
        if denominator <= 0.0:
            value, error = None, None
        else:
            value, error = experimental.ratio_with_covariance(
                numerator,
                numerator_variance,
                denominator,
                denominator_variance,
                covariance,
            )
        values.append(value)
        errors.append(error)
        covariances.append(covariance)
    return {
        "edges": list(first.edges),
        "values": values,
        "errors": errors,
        "numerator_denominator_covariance": covariances,
    }


def pt2_slope_target_combination(
    spectra: Mapping[str, experimental.BinSeries],
    target_weights: Mapping[str, float],
    pt2_edges: Sequence[float],
    slope_cells: int,
    minimum_positive_bins: int = 3,
) -> dict[str, Any]:
    """Fit ``A exp(-pT2/slope)`` in every target-combined spectrum cell.

    The inclusive-DIS denominator and z-bin width are constant throughout one
    fitted spectrum, so neither changes its exponential inverse slope.  The
    fit therefore uses the normalized, signed-NLO hadron yield directly after
    the target sum.  A weighted straight-line fit to ``log(dY/dpT2)`` gives the
    slope and its Monte Carlo error.
    """

    labels = set(target_weights)
    first = _check_series((spectra,), labels)
    bins_per_cell = len(pt2_edges) - 1
    if (
        bins_per_cell < minimum_positive_bins
        or slope_cells <= 0
        or len(first.values) != slope_cells * bins_per_cell
        or any(
            not math.isfinite(float(high))
            or float(high) <= float(low)
            for low, high in zip(pt2_edges, pt2_edges[1:])
        )
    ):
        raise experimental.CampaignError(
            "SIDIS pT2 spectrum layout does not match the slope fit"
        )

    combined_values = [
        sum(
            float(weight) * spectra[target].values[index]
            for target, weight in target_weights.items()
        )
        for index in range(len(first.values))
    ]
    combined_variances = [
        sum(
            float(weight) ** 2 * spectra[target].variances[index]
            for target, weight in target_weights.items()
        )
        for index in range(len(first.values))
    ]
    values: list[float | None] = []
    errors: list[float | None] = []
    retained_bins: list[list[int]] = []
    for cell in range(slope_cells):
        x_values: list[float] = []
        log_values: list[float] = []
        weights: list[float] = []
        retained: list[int] = []
        for pt2_bin, (low, high) in enumerate(zip(pt2_edges, pt2_edges[1:])):
            index = cell * bins_per_cell + pt2_bin
            width = float(high) - float(low)
            value = combined_values[index] / width
            variance = combined_variances[index] / (width * width)
            if value <= 0.0 or variance <= 0.0 or not all(
                math.isfinite(item) for item in (value, variance)
            ):
                continue
            x_values.append(.5 * (float(low) + float(high)))
            log_values.append(math.log(value))
            weights.append(value * value / variance)
            retained.append(pt2_bin)
        retained_bins.append(retained)
        if len(retained) < minimum_positive_bins:
            values.append(None)
            errors.append(None)
            continue
        total_weight = sum(weights)
        weighted_x = sum(w * x for w, x in zip(weights, x_values))
        weighted_y = sum(w * y for w, y in zip(weights, log_values))
        weighted_xx = sum(w * x * x for w, x in zip(weights, x_values))
        weighted_xy = sum(
            w * x * y for w, x, y in zip(weights, x_values, log_values)
        )
        determinant = total_weight * weighted_xx - weighted_x * weighted_x
        if determinant <= 0.0:
            values.append(None)
            errors.append(None)
            continue
        exponent = (
            total_weight * weighted_xy - weighted_x * weighted_y
        ) / determinant
        if exponent >= 0.0 or not math.isfinite(exponent):
            values.append(None)
            errors.append(None)
            continue
        exponent_variance = total_weight / determinant
        slope = -1.0 / exponent
        slope_error = math.sqrt(exponent_variance) / (exponent * exponent)
        values.append(slope)
        errors.append(slope_error)
    return {
        "edges": [float(index) for index in range(slope_cells + 1)],
        "values": values,
        "errors": errors,
        "retained_pt2_bins": retained_bins,
    }


def charge_ratio_target_combination(
    negative_yields: Mapping[str, experimental.BinSeries],
    positive_yields: Mapping[str, experimental.BinSeries],
    covariance_proxies: Mapping[str, experimental.BinSeries],
    target_weights: Mapping[str, float],
) -> dict[str, Any]:
    """Form a negative/positive yield ratio after target combination.

    The two charge yields and their same-event covariance are combined across
    proton/neutron components before division.  Unit density widths make this
    an ordinary dimensionless ratio while retaining the tested normalized-bin
    signed-NLO arithmetic used by the multiplicity estimator.
    """

    first = _check_series(
        (negative_yields, positive_yields, covariance_proxies),
        set(target_weights),
    )
    result = multiplicity_target_combination(
        negative_yields, positive_yields, covariance_proxies,
        target_weights, [1.0] * len(first.values),
    )
    result["combined_negative_yield"] = result.pop("isoscalar_numerator")
    result["combined_positive_yield"] = result.pop("isoscalar_denominator")
    return result


def _covariance_result(
    residual: Sequence[float],
    covariance: Any,
    retained: Sequence[str],
    description: str,
) -> dict[str, Any]:
    try:
        import numpy as np
    except ImportError as exc:
        raise experimental.CampaignError(
            "NumPy is required for COMPASS covariance calculations"
        ) from exc
    matrix = np.asarray(covariance, dtype=float)
    vector = np.asarray(residual, dtype=float)
    inverse = np.linalg.pinv(matrix, hermitian=True, rcond=1.0e-12)
    chi2 = float(vector @ inverse @ vector)
    try:
        decorrelated = np.linalg.solve(np.linalg.cholesky(matrix), vector)
    except np.linalg.LinAlgError:
        decorrelated = np.full(len(vector), np.nan)
    return {
        "points": len(vector),
        "retained_points": list(retained),
        "chi2_correlated": chi2,
        "covariance": description,
        "decorrelated_pulls": [
            float(value) if math.isfinite(float(value)) else None
            for value in decorrelated
        ],
    }


def _diagonal_plus_rank_one_result(
    residual: Sequence[float],
    diagonal_variance: Sequence[float],
    nuisance: Sequence[float],
    retained: Sequence[str],
    description: str,
) -> dict[str, Any]:
    """Evaluate ``diag(d) + u u^T`` without constructing a dense matrix.

    The published 2017 systematic model has exactly this form.  The
    Sherman--Morrison identity gives the correlated chi-square in linear
    time.  The reported pulls use the symmetric inverse square root of the
    whitened rank-one covariance, so their squared norm is the same
    correlated chi-square.
    """

    try:
        import numpy as np
    except ImportError as exc:
        raise experimental.CampaignError(
            "NumPy is required for COMPASS covariance calculations"
        ) from exc
    vector = np.asarray(residual, dtype=float)
    diagonal = np.asarray(diagonal_variance, dtype=float)
    correlated = np.asarray(nuisance, dtype=float)
    if (
        vector.shape != diagonal.shape
        or vector.shape != correlated.shape
        or np.any(~np.isfinite(diagonal))
        or np.any(diagonal <= 0.0)
    ):
        raise experimental.CampaignError(
            "COMPASS diagonal covariance must be finite, positive, and "
            "match the retained prediction bins"
        )

    whitened_residual = vector / np.sqrt(diagonal)
    whitened_nuisance = correlated / np.sqrt(diagonal)
    nuisance_norm2 = float(whitened_nuisance @ whitened_nuisance)
    projection = float(whitened_nuisance @ whitened_residual)
    chi2 = float(
        whitened_residual @ whitened_residual
        - projection**2 / (1.0 + nuisance_norm2)
    )
    if nuisance_norm2 > 0.0:
        inverse_sqrt_shift = (
            1.0 / math.sqrt(1.0 + nuisance_norm2) - 1.0
        ) / nuisance_norm2
        decorrelated = (
            whitened_residual
            + inverse_sqrt_shift * projection * whitened_nuisance
        )
    else:
        decorrelated = whitened_residual
    return {
        "points": len(vector),
        "retained_points": list(retained),
        "chi2_correlated": max(chi2, 0.0),
        "covariance": description,
        "decorrelated_pulls": [float(value) for value in decorrelated],
    }


def a1_goodness_of_fit(
    predictions: Mapping[str, Mapping[str, Any]],
    snapshot: Mapping[str, Any],
) -> dict[str, Any]:
    """Evaluate the 40 COMPASS 2009 points with the published covariance."""

    try:
        import numpy as np
    except ImportError as exc:
        raise experimental.CampaignError(
            "NumPy is required for COMPASS covariance calculations"
        ) from exc
    order = list(snapshot["point_order"])
    data_by_label: dict[str, Mapping[str, Any]] = {}
    prediction_by_label: dict[str, tuple[float | None, float | None]] = {}
    for species, dataset in snapshot["datasets"].items():
        prediction = predictions[species]
        for index, point in enumerate(dataset["points"]):
            label = f"{species}:x{index + 1:02d}"
            data_by_label[label] = point
            prediction_by_label[label] = (
                prediction["values"][index], prediction["errors"][index]
            )

    retained_indices: list[int] = []
    retained_labels: list[str] = []
    residual: list[float] = []
    mc_variance: list[float] = []
    diagonal_systematic: list[float] = []
    correlated_systematic: list[float] = []
    for global_index, label in enumerate(order):
        theory, error = prediction_by_label[label]
        if theory is None or error is None:
            continue
        point = data_by_label[label]
        retained_indices.append(global_index)
        retained_labels.append(label)
        residual.append(float(theory) - float(point["a1"]))
        mc_variance.append(float(error) ** 2)
        diagonal_systematic.append(float(point["systematic_uncorrelated"]) ** 2)
        correlated_systematic.append(float(point.get(
            "systematic_correlated_8pct",
            point.get("systematic_correlated_6pct", 0.0),
        )))
    if not retained_indices:
        return {"points": 0, "status": "no finite theory bins"}

    published = np.asarray(snapshot["statistical_covariance"], dtype=float)
    covariance = published[np.ix_(retained_indices, retained_indices)]
    covariance += np.diag(np.asarray(mc_variance) + np.asarray(diagonal_systematic))
    nuisance = np.asarray(correlated_systematic)
    covariance += np.outer(nuisance, nuisance)
    return _covariance_result(
        residual,
        covariance,
        retained_labels,
        "published per-x 4x4 statistical blocks plus diagonal residual "
        "systematics and MC statistics plus one shared multiplicative nuisance",
    )


def multiplicity_goodness_of_fit(
    predictions: Mapping[str, Mapping[str, Any]],
    snapshot: Mapping[str, Any],
) -> dict[str, Any]:
    """Evaluate one COMPASS 2017 record with its published error split."""
    retained: list[str] = []
    residual: list[float] = []
    diagonal_variance: list[float] = []
    nuisance: list[float] = []
    for species in sorted(snapshot["datasets"]):
        dataset = snapshot["datasets"][species]
        prediction = predictions[species]
        for index, point in enumerate(dataset["points"]):
            theory = prediction["values"][index]
            error = prediction["errors"][index]
            if theory is None or error is None:
                continue
            retained.append(f"{species}:cell{index + 1:03d}")
            residual.append(float(theory) - float(point["value"]))
            diagonal_variance.append(
                float(point["stat"]) ** 2
                + float(point["systematic_uncorrelated_60pct"]) ** 2
                + float(error) ** 2
            )
            nuisance.append(float(point["systematic_correlated_80pct"]))
    if not retained:
        return {"points": 0, "status": "no finite theory bins"}
    return _diagonal_plus_rank_one_result(
        residual,
        diagonal_variance,
        nuisance,
        retained,
        "diagonal published statistical errors, diagonal 0.6 systematic "
        "component, MC statistics, and one record-wide 0.8 systematic nuisance",
    )


def diagonal_multiplicity_goodness_of_fit(
    predictions: Mapping[str, Mapping[str, Any]],
    snapshot: Mapping[str, Any],
) -> dict[str, Any]:
    """Evaluate multiplicities for records with diagonal released errors."""
    results: dict[str, Any] = {}
    for dataset_id, dataset in snapshot["datasets"].items():
        prediction = predictions[dataset_id]
        residual: list[float] = []
        variance: list[float] = []
        retained: list[str] = []
        for index, point in enumerate(dataset["points"]):
            theory, error = prediction["values"][index], prediction["errors"][index]
            if theory is None or error is None:
                continue
            residual.append(float(theory) - float(point["value"]))
            variance.append(float(point["stat"])**2 + float(point["systematic"])**2
                            + float(error)**2)
            retained.append(f"{dataset_id}:cell{index + 1:04d}")
        if residual:
            if any(not math.isfinite(value) or value <= 0.0 for value in variance):
                raise experimental.CampaignError(
                    f"Non-positive diagonal variance in {dataset_id}"
                )
            pulls = [delta / math.sqrt(var) for delta, var in zip(residual, variance)]
            results[dataset_id] = {
                "points": len(retained),
                "retained_points": retained,
                "chi2_correlated": sum(value * value for value in pulls),
                "covariance": (
                    "diagonal published statistical/systematic errors plus "
                    "MC statistics"
                ),
                "decorrelated_pulls": pulls,
            }
        else:
            results[dataset_id] = {"points": 0, "status": "no finite theory bins"}
    return {
        "independent_datasets": results,
        "policy": (
            "datasets are reported separately because no cross-dataset "
            "covariance was released"
        ),
    }


def charge_ratio_goodness_of_fit(
    predictions: Mapping[str, Mapping[str, Any]],
    snapshot: Mapping[str, Any],
) -> dict[str, Any]:
    """Evaluate each overlapping COMPASS 2020 paper table independently."""

    results: dict[str, Any] = {}
    for dataset_id, dataset in snapshot["datasets"].items():
        prediction = predictions[dataset_id]
        retained: list[str] = []
        residual: list[float] = []
        diagonal_variance: list[float] = []
        nuisance: list[float] = []
        for index, point in enumerate(dataset["points"]):
            theory = prediction["values"][index]
            error = prediction["errors"][index]
            if theory is None or error is None:
                continue
            retained.append(f"{dataset_id}:cell{index + 1:02d}")
            residual.append(float(theory) - float(point["value"]))
            diagonal_variance.append(
                float(point["stat"])**2
                + float(point["systematic_uncorrelated_half"])**2
                + float(error)**2
            )
            nuisance.append(float(point["systematic_correlated_sqrt75"]))
        if retained:
            results[dataset_id] = _diagonal_plus_rank_one_result(
                residual, diagonal_variance, nuisance, retained,
                "diagonal published statistical errors, diagonal half-"
                "systematic component, MC statistics, and one sqrt(0.75) "
                "systematic nuisance within this paper table",
            )
        else:
            results[dataset_id] = {
                "points": 0, "status": "no finite theory bins"
            }
    return {
        "independent_datasets": results,
        "policy": (
            "the overlapping x-z and z-momentum tables are never combined "
            "into one chi-square"
        ),
    }


def hermes_multiplicity_goodness_of_fit(
    predictions: Mapping[str, Mapping[str, Any]],
    snapshot: Mapping[str, Any],
) -> dict[str, Any]:
    """Evaluate each HERMES target/species/binning covariance independently."""
    try:
        import numpy as np
    except ImportError as exc:
        raise experimental.CampaignError(
            "NumPy is required for HERMES covariance calculations"
        ) from exc
    results: dict[str, Any] = {}
    for dataset_id, dataset in snapshot["datasets"].items():
        prediction = predictions[dataset_id]
        indices: list[int] = []
        residual: list[float] = []
        retained: list[str] = []
        diagonal: list[float] = []
        for index, point in enumerate(dataset["points"]):
            theory, error = prediction["values"][index], prediction["errors"][index]
            if theory is None or error is None:
                continue
            indices.append(index)
            residual.append(float(theory) - float(point["value"]))
            retained.append(f"{dataset_id}:cell{index + 1:03d}")
            diagonal.append(float(point["systematic"])**2 + float(error)**2)
        if not indices:
            results[dataset_id] = {"points": 0, "status": "no finite theory bins"}
            continue
        published = np.asarray(dataset["statistical_covariance"], dtype=float)
        covariance = published[np.ix_(indices, indices)] + np.diag(diagonal)
        results[dataset_id] = _covariance_result(
            residual, covariance, retained,
            "released dense statistical covariance plus point-to-point systematic and MC variances",
        )
    return {"independent_datasets": results,
            "policy": "the five overlapping binnings are never combined into one chi-square"}
