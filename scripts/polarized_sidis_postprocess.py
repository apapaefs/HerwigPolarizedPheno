#!/usr/bin/env python3
"""Normalized-yield arithmetic for polarized SIDIS campaign samples.

The functions in this module are deliberately independent of Rivet and
Herwig execution.  They operate on the normalized ``BinSeries`` objects read
from completed campaign shards, so POSNLO/NEGNLO and target-component
combinations are performed before any asymmetry is constructed.
"""

from __future__ import annotations

import math
from typing import Any, Mapping, Sequence

import run_experimental_campaign as experimental


HELICITIES = ("PP", "PM", "MP", "MM")
LONGITUDINAL_SIGN = {"PP": 1.0, "PM": -1.0, "MP": -1.0, "MM": 1.0}
SINGLE_SPIN_A_SIGN = {"PP": 1.0, "PM": 1.0, "MP": -1.0, "MM": -1.0}
SINGLE_SPIN_B_SIGN = {"PP": 1.0, "PM": -1.0, "MP": 1.0, "MM": -1.0}


def _same_edges(
    first: experimental.BinSeries, second: experimental.BinSeries
) -> None:
    if not experimental._same_edges(first.edges, second.edges):
        raise experimental.CampaignError(
            "Cannot combine SIDIS series with different bin edges"
        )


def add_nlo_components(
    contributions: Mapping[str, experimental.BinSeries],
) -> experimental.BinSeries:
    """Add normalized POSNLO and NEGNLO cross sections."""

    if set(contributions) != {"POSNLO", "NEGNLO"}:
        raise experimental.CampaignError(
            "SIDIS NLO combination requires POSNLO and NEGNLO"
        )
    return experimental.add_independent_series(
        [contributions["POSNLO"], contributions["NEGNLO"]]
    )


def difference_series(
    positive: experimental.BinSeries,
    negative: experimental.BinSeries,
    covariance_proxy: experimental.BinSeries,
) -> experimental.BinSeries:
    """Form a charge-difference yield with its same-event covariance.

    The covariance proxy is filled with ``sqrt(n_+ n_-)`` at event level.
    Consequently its post-normalization variance is the covariance of the two
    multiplicity-weighted yields in each bin.
    """

    _same_edges(positive, negative)
    _same_edges(positive, covariance_proxy)
    variances: list[float] = []
    for positive_variance, negative_variance, covariance in zip(
        positive.variances, negative.variances, covariance_proxy.variances
    ):
        variance = positive_variance + negative_variance - 2.0 * covariance
        tolerance = 1.0e-12 * max(
            1.0, positive_variance, negative_variance, abs(covariance)
        )
        if variance < 0.0 and abs(variance) <= tolerance:
            variance = 0.0
        if variance < 0.0:
            raise experimental.CampaignError(
                "Charge-difference covariance gives a negative variance"
            )
        variances.append(variance)
    return experimental.BinSeries(
        list(positive.edges),
        [
            positive_value - negative_value
            for positive_value, negative_value in zip(
                positive.values, negative.values
            )
        ],
        variances,
    )


def target_coefficients(
    target: str,
) -> tuple[dict[str, float], dict[str, float]]:
    """Return denominator and longitudinal-numerator coefficients.

    Labels are ``P:PP`` etc.  The factor of one quarter defining each helicity
    average is retained explicitly.  For deuterium, each nucleon component
    receives the additional impulse-approximation factor one half.  HERMES
    already divides its published deuteron asymmetries by the nucleon-to-nucleus
    polarization ratio ``f_D=0.926``, so a free-proton-plus-free-neutron theory
    prediction must not apply that factor again.
    """

    if target == "proton":
        denominator = {f"P:{helicity}": 0.25 for helicity in HELICITIES}
    elif target == "deuteron":
        denominator = {
            f"{component}:{helicity}": 0.125
            for component in ("P", "N")
            for helicity in HELICITIES
        }
    else:
        raise experimental.CampaignError(f"Unknown SIDIS target {target!r}")
    longitudinal = {
        label: coefficient
        * LONGITUDINAL_SIGN[label.split(":", 1)[1]]
        for label, coefficient in denominator.items()
    }
    return denominator, longitudinal


def _linear_value(
    samples: Mapping[str, experimental.BinSeries],
    coefficients: Mapping[str, float],
    index: int,
) -> tuple[float, float]:
    return (
        sum(
            coefficients[label] * samples[label].values[index]
            for label in coefficients
        ),
        sum(
            coefficients[label] ** 2 * samples[label].variances[index]
            for label in coefficients
        ),
    )


def asymmetry(
    samples: Mapping[str, experimental.BinSeries],
    target: str,
) -> dict[str, Any]:
    """Construct A_parallel and parity/single-spin closure diagnostics."""

    denominator_coefficients, numerator_coefficients = target_coefficients(target)
    if set(samples) != set(denominator_coefficients):
        raise experimental.CampaignError(
            f"SIDIS {target} sample labels do not match the target model"
        )
    first = next(iter(samples.values()))
    if any(
        not experimental._same_edges(first.edges, item.edges)
        for item in samples.values()
    ):
        raise experimental.CampaignError(
            "SIDIS target components have inconsistent binning"
        )

    values: list[float | None] = []
    errors: list[float | None] = []
    sigma_uu_values: list[float] = []
    sigma_uu_variances: list[float] = []
    sigma_ll_values: list[float] = []
    sigma_ll_variances: list[float] = []
    for index in range(len(first.values)):
        numerator, numerator_variance = _linear_value(
            samples, numerator_coefficients, index
        )
        denominator, denominator_variance = _linear_value(
            samples, denominator_coefficients, index
        )
        covariance = sum(
            numerator_coefficients[label]
            * denominator_coefficients[label]
            * samples[label].variances[index]
            for label in samples
        )
        value, error = experimental.ratio_with_covariance(
            numerator,
            numerator_variance,
            denominator,
            denominator_variance,
            covariance,
        )
        values.append(value)
        errors.append(error)
        sigma_uu_values.append(denominator)
        sigma_uu_variances.append(denominator_variance)
        sigma_ll_values.append(numerator)
        sigma_ll_variances.append(numerator_variance)

    closures: dict[str, dict[str, Any]] = {}
    component_weights = {"P": 1.0} if target == "proton" else {"P": 0.5, "N": 0.5}
    target_helicity: dict[str, experimental.BinSeries] = {}
    for helicity in HELICITIES:
        selected = {
            component: samples[f"{component}:{helicity}"]
            for component in component_weights
        }
        target_helicity[helicity] = experimental.linear_combine_series(
            selected, component_weights
        )
    for first_helicity, second_helicity, label in (
        ("PP", "MM", "Parity_PP_MM"),
        ("PM", "MP", "Parity_PM_MP"),
    ):
        closure_values: list[float | None] = []
        closure_errors: list[float | None] = []
        for index in range(len(first.values)):
            value, error = experimental.parity_residual(
                target_helicity[first_helicity].values[index],
                target_helicity[first_helicity].variances[index],
                target_helicity[second_helicity].values[index],
                target_helicity[second_helicity].variances[index],
            )
            closure_values.append(value)
            closure_errors.append(error)
        closures[label] = {
            "edges": list(first.edges),
            "values": closure_values,
            "errors": closure_errors,
        }
    for label, signs in (
        ("SingleSpinA", SINGLE_SPIN_A_SIGN),
        ("SingleSpinB", SINGLE_SPIN_B_SIGN),
    ):
        numerator = {
            helicity: 0.25 * signs[helicity] for helicity in HELICITIES
        }
        denominator = {helicity: 0.25 for helicity in HELICITIES}
        closure_values: list[float | None] = []
        closure_errors: list[float | None] = []
        for index in range(len(first.values)):
            numerator_value, numerator_variance = _linear_value(
                target_helicity, numerator, index
            )
            denominator_value, denominator_variance = _linear_value(
                target_helicity, denominator, index
            )
            covariance = sum(
                numerator[helicity]
                * denominator[helicity]
                * target_helicity[helicity].variances[index]
                for helicity in HELICITIES
            )
            value, error = experimental.ratio_with_covariance(
                numerator_value,
                numerator_variance,
                denominator_value,
                denominator_variance,
                covariance,
            )
            closure_values.append(value)
            closure_errors.append(error)
        closures[label] = {
            "edges": list(first.edges),
            "values": closure_values,
            "errors": closure_errors,
        }

    return {
        "edges": list(first.edges),
        "values": values,
        "errors": errors,
        "sigma_uu": experimental.BinSeries(
            list(first.edges), sigma_uu_values, sigma_uu_variances
        ),
        "sigma_ll": experimental.BinSeries(
            list(first.edges), sigma_ll_values, sigma_ll_variances
        ),
        "closures": closures,
    }


def azimuthal_moment(
    numerators: Mapping[str, experimental.BinSeries],
    denominators: Mapping[str, experimental.BinSeries],
    covariances: Mapping[str, Sequence[float]],
    target: str,
) -> dict[str, Any]:
    """Form the unpolarized ``2<cos(phi)>`` diagnostic."""

    coefficients, _ = target_coefficients(target)
    if set(numerators) != set(coefficients) or set(denominators) != set(coefficients):
        raise experimental.CampaignError(
            "Azimuthal-moment components do not match target coefficients"
        )
    numerator = experimental.linear_combine_series(numerators, coefficients)
    denominator = experimental.linear_combine_series(denominators, coefficients)
    covariance = [
        sum(
            coefficients[label] ** 2 * float(covariances[label][index])
            for label in coefficients
        )
        for index in range(len(numerator.values))
    ]
    values: list[float | None] = []
    errors: list[float | None] = []
    for index in range(len(numerator.values)):
        value, error = experimental.ratio_with_covariance(
            numerator.values[index],
            numerator.variances[index],
            denominator.values[index],
            denominator.variances[index],
            covariance[index],
        )
        values.append(value)
        errors.append(error)
    return {"edges": list(numerator.edges), "values": values, "errors": errors}


def projection_goodness_of_fit(
    prediction: Mapping[str, Any], dataset: Mapping[str, Any]
) -> dict[str, Any]:
    """Evaluate one projection with its own published covariance.

    Cross-projection covariance is not published by HERMES, so this function
    intentionally returns one result per dataset and has no aggregation mode.
    Point-to-point systematic uncertainties and MC statistical errors are
    added on the diagonal of the published statistical covariance.
    """

    try:
        import numpy as np
    except ImportError as exc:
        raise experimental.CampaignError(
            "NumPy is required for covariance-aware SIDIS goodness of fit"
        ) from exc

    selected_theory: list[float] = []
    selected_data: list[float] = []
    selected_mc_variance: list[float] = []
    selected_systematics: list[float] = []
    retained: list[int] = []
    for point_index, point in enumerate(dataset["points"]):
        flat_bin = int(point.get("flat_bin", point_index))
        if flat_bin >= len(prediction["values"]):
            raise experimental.CampaignError(
                f"Published flat bin {flat_bin} is outside {dataset['id']}"
            )
        theory = prediction["values"][flat_bin]
        theory_error = prediction["errors"][flat_bin]
        if theory is None or theory_error is None:
            continue
        selected_theory.append(float(theory))
        selected_data.append(float(point["aparallel"]))
        selected_mc_variance.append(float(theory_error) ** 2)
        selected_systematics.append(float(point["aparallel_systematic"]))
        retained.append(point_index)
    if not retained:
        return {"points": 0, "status": "no finite theory bins"}

    published = dataset.get("aparallel_statistical_covariance")
    if published:
        statistical = np.asarray(published, dtype=float)[np.ix_(retained, retained)]
        covariance_source = "published statistical covariance"
    else:
        statistical = np.diag(
            [
                float(dataset["points"][index]["aparallel_stat"]) ** 2
                for index in retained
            ]
        )
        covariance_source = "diagonal statistical errors; covariance not supplied"
    covariance = statistical + np.diag(
        np.asarray(selected_systematics, dtype=float) ** 2
        + np.asarray(selected_mc_variance, dtype=float)
    )
    inverse = np.linalg.pinv(covariance, hermitian=True, rcond=1.0e-12)
    residual = np.asarray(selected_theory) - np.asarray(selected_data)
    chi2 = float(residual @ inverse @ residual)
    try:
        decorrelated = np.linalg.solve(np.linalg.cholesky(covariance), residual)
    except np.linalg.LinAlgError:
        decorrelated = np.full(len(residual), np.nan)
    return {
        "points": len(retained),
        "retained_point_indices": [index + 1 for index in retained],
        "chi2_correlated": chi2,
        "chi2_profiled": chi2,
        "global_nuisances": "none published for this measurement",
        "covariance": (
            covariance_source
            + " plus diagonal point-to-point systematics and MC statistics"
        ),
        "decorrelated_pulls": [
            None if not math.isfinite(float(value)) else float(value)
            for value in decorrelated
        ],
    }
