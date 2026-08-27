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


def a1_deuteron(
    ordinary: Mapping[str, experimental.BinSeries],
    inverse_depolarization: Mapping[str, experimental.BinSeries],
    covariance_proxy: Mapping[str, experimental.BinSeries],
) -> dict[str, Any]:
    """Form the COMPASS 2009 deuteron ``A1`` estimator.

    The ordinary yield is the unpolarized denominator.  The longitudinal
    helicity difference is formed from the event-wise inverse-depolarization
    yield.  Proton and neutron each carry the isoscalar factor one half, every
    helicity carries the averaging factor one quarter, and the longitudinal
    numerator alone carries the deuteron factor 0.925.
    """

    labels = {
        f"{target}:{helicity}"
        for target in ("P", "N")
        for helicity in HELICITIES
    }
    first = _check_series(
        (ordinary, inverse_depolarization, covariance_proxy), labels
    )
    denominator_coefficients = {label: 0.125 for label in labels}
    numerator_coefficients = {
        label: 0.125
        * DEUTERON_POLARIZATION_FACTOR
        * LONGITUDINAL_SIGN[label.split(":", 1)[1]]
        for label in labels
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


def multiplicity_isoscalar(
    numerators: Mapping[str, experimental.BinSeries],
    denominators: Mapping[str, experimental.BinSeries],
    covariance_proxies: Mapping[str, experimental.BinSeries],
    z_widths: Sequence[float],
) -> dict[str, Any]:
    """Form ``dM/dz`` after the proton/neutron isoscalar sum.

    Non-finite or non-positive inclusive-DIS denominators are represented by
    masked bins.  The same-event numerator--denominator covariance is retained
    before dividing the ratio and its Monte Carlo error by the published cell
    width in ``z``.
    """

    labels = {"P", "N"}
    first = _check_series((numerators, denominators, covariance_proxies), labels)
    if len(z_widths) != len(first.values) or any(
        not math.isfinite(float(width)) or float(width) <= 0.0
        for width in z_widths
    ):
        raise experimental.CampaignError(
            "COMPASS multiplicity z widths do not match the flattened cells"
        )

    values: list[float | None] = []
    errors: list[float | None] = []
    numerator_values: list[float] = []
    denominator_values: list[float] = []
    covariances: list[float] = []
    for index, width in enumerate(z_widths):
        numerator = 0.5 * (
            numerators["P"].values[index] + numerators["N"].values[index]
        )
        numerator_variance = 0.25 * (
            numerators["P"].variances[index]
            + numerators["N"].variances[index]
        )
        denominator = 0.5 * (
            denominators["P"].values[index] + denominators["N"].values[index]
        )
        denominator_variance = 0.25 * (
            denominators["P"].variances[index]
            + denominators["N"].variances[index]
        )
        covariance = 0.25 * (
            covariance_proxies["P"].variances[index]
            + covariance_proxies["N"].variances[index]
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
        correlated_systematic.append(float(point["systematic_correlated_8pct"]))
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
        "systematics and MC statistics plus one shared 8% multiplicative nuisance",
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
