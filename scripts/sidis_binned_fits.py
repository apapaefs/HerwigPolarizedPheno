"""Binned COMPASS estimators after normalized signed-NLO and target sums.

Event second moments include correlations between bins and, for azimuths,
the epsilon sum. No logarithm or sign-dependent removal of spectrum bins is
used. These fits implement published functions/ranges at generator level;
they do not emulate detector response or claim the experiment's unpublished
minimizer settings.
"""
from __future__ import annotations

import math
from typing import Any, Mapping, Sequence

import numpy as np
import run_experimental_campaign as experimental


def _moments(samples, target_weights, dimensions, cells, covariance_proxies):
    labels = set(target_weights)
    if not labels or set(samples) != labels:
        raise experimental.CampaignError("Binned-fit target labels do not match")
    size = cells * dimensions
    first = next(iter(samples.values()))
    values = np.zeros(size)
    variances = np.zeros(size)
    for label, series in samples.items():
        weight = float(target_weights[label])
        if (not math.isfinite(weight) or len(series.values) != size
                or len(series.variances) != size
                or not experimental._same_edges(series.edges, first.edges)):
            raise experimental.CampaignError("Invalid binned-fit input layout")
        values += weight * np.asarray(series.values)
        variances += weight**2 * np.asarray(series.variances)
    covariance = np.zeros((cells, dimensions, dimensions))
    diagonal = np.arange(dimensions)
    covariance[:, diagonal, diagonal] = variances.reshape(cells, dimensions)
    if covariance_proxies is not None:
        pairs = dimensions * (dimensions - 1) // 2
        if set(covariance_proxies) != labels:
            raise experimental.CampaignError("Missing binned-fit covariance target")
        first_covariance = next(iter(covariance_proxies.values()))
        off_diagonal = np.zeros(cells * pairs)
        for label, series in covariance_proxies.items():
            if (len(series.variances) != cells * pairs
                    or not experimental._same_edges(series.edges, first_covariance.edges)):
                raise experimental.CampaignError("Invalid binned-fit covariance layout")
            # The proxy's variance, not its central value, stores the product.
            off_diagonal += float(target_weights[label])**2 * np.asarray(series.variances)
        upper = np.triu_indices(dimensions, 1)
        covariance[:, upper[0], upper[1]] = off_diagonal.reshape(cells, pairs)
        covariance[:, upper[1], upper[0]] = off_diagonal.reshape(cells, pairs)
    return values.reshape(cells, dimensions), covariance


def _precision(covariance):
    if not np.all(np.isfinite(covariance)) or np.any(np.diag(covariance) <= 0.):
        raise ValueError("missing_variance_in_fit_range")
    scale = np.sqrt(np.diag(covariance))
    correlation = covariance / np.outer(scale, scale)
    if np.linalg.cond(correlation) > 1.e12:
        raise ValueError("singular_covariance")
    try:
        np.linalg.cholesky(correlation)
        return np.linalg.solve(correlation, np.eye(len(scale))) / np.outer(scale, scale)
    except np.linalg.LinAlgError as exc:
        raise ValueError("nonpositive_covariance") from exc


def azimuthal_design(phi_edges: Sequence[float]):
    """Average each fitted harmonic over the actual histogram bin."""
    edges = np.asarray(phi_edges, dtype=float)
    widths = np.diff(edges)
    return np.column_stack((
        np.ones(len(widths)),
        (np.sin(edges[1:]) - np.sin(edges[:-1])) / widths,
        (np.sin(2.*edges[1:]) - np.sin(2.*edges[:-1])) / (2.*widths),
        (np.cos(edges[:-1]) - np.cos(edges[1:])) / widths,
    ))


def _azimuthal_fit(inputs, covariance, harmonic, design, retained):
    if not np.all(np.isfinite(inputs)):
        raise ValueError("nonfinite_inputs")
    counts, epsilon_sum = inputs[:-1], float(inputs[-1])
    total = float(np.sum(counts))
    if total <= 0. or epsilon_sum <= 0.:
        raise ValueError("nonpositive_total_or_epsilon")
    precision = _precision(covariance[np.ix_(retained, retained)])
    x = design[retained]
    try:
        response = np.linalg.solve(x.T @ precision @ x, x.T @ precision)
    except np.linalg.LinAlgError as exc:
        raise ValueError("singular_harmonic_fit") from exc
    beta = response @ counts[retained]
    if beta[0] <= 0.:
        raise ValueError("nonpositive_fit_normalization")
    residual = counts[retained] - x @ beta
    coefficient = float(beta[harmonic] / beta[0])
    value = coefficient * total / epsilon_sum
    # A = beta_h * N / (beta_0 * E). N and E share events with fitted counts.
    gradient = np.zeros(len(inputs))
    gradient[:-1] = coefficient / epsilon_sum
    gradient[retained] += total / epsilon_sum * (
        response[harmonic] / beta[0]
        - beta[harmonic] * response[0] / beta[0]**2)
    gradient[-1] = -coefficient * total / epsilon_sum**2
    variance = float(gradient @ covariance @ gradient)
    if not math.isfinite(variance) or variance < -1.e-10 * max(1., value*value):
        raise ValueError("invalid_amplitude_variance")
    diagnostics = {
        "status": "ok", "chi2": max(0., float(residual @ precision @ residual)),
        "ndof": len(retained) - 4, "retained_phi_bins": retained.tolist(),
        "epsilon_mean": epsilon_sum / total,
        "normalization": float(beta[0]),
        "cos1_coefficient": float(beta[1]/beta[0]),
        "cos2_coefficient": float(beta[2]/beta[0]),
        "sin_coefficient_nuisance": float(beta[3]/beta[0]),
    }
    return value, math.sqrt(max(0., variance)), diagnostics, gradient


def azimuthal_target_fit(samples, covariance_proxies, target_weights,
                         harmonic: int, cells: int, policy: Mapping[str, Any]):
    if harmonic not in (1, 2):
        raise experimental.CampaignError("Unsupported COMPASS harmonic")
    edges = policy["phi_edges"]
    if (len(edges) != 17 or not np.allclose(edges, np.linspace(0., 2.*np.pi, 17))
            or policy["excluded_phi_bins"] != [0, 15]):
        raise experimental.CampaignError("COMPASS published phi fit layout changed")
    means, covariance = _moments(samples, target_weights, 17, cells, covariance_proxies)
    design = azimuthal_design(edges)
    retained = np.arange(1, 15)
    values, errors, diagnostics = [], [], []
    for inputs, cov in zip(means, covariance):
        try:
            value, error, detail, _ = _azimuthal_fit(inputs, cov, harmonic, design, retained)
        except ValueError as exc:
            value, error, detail = None, None, {"status": str(exc)}
        values.append(value)
        errors.append(error)
        diagnostics.append(detail)
    return {"edges": [float(i) for i in range(cells+1)], "values": values,
            "errors": errors, "fit_diagnostics": diagnostics,
            "covariance_model": "event second moments including phi/epsilon covariance"}


def exponential_integrals(edges, slope):
    low, high = np.asarray(edges[:-1]), np.asarray(edges[1:])
    # expm1 avoids cancellation for slopes large compared with a bin width.
    return slope * np.exp(-low/slope) * (-np.expm1(-(high-low)/slope))


def _slope_fit(yields, covariance, edges, minimum_region_significance):
    if not np.all(np.isfinite(yields)):
        raise ValueError("nonfinite_inputs")
    precision = _precision(covariance)  # Every bin of the full fit range is used.
    # Test coverage in four contiguous regions, without discarding negative bins.
    # This rejects a positive low-pT fragment plus a noise-only/negative tail.
    region_significances = []
    for region in np.array_split(np.arange(len(yields)), 4):
        integral = float(np.sum(yields[region]))
        variance = float(np.sum(covariance[np.ix_(region, region)]))
        significance = integral / math.sqrt(variance) if variance > 0. else -math.inf
        region_significances.append(significance)
    if min(region_significances) < minimum_region_significance:
        raise ValueError("insufficient_signal_across_fit_range")

    weighted_y = precision @ yields
    def objective(log_slope):
        shape = exponential_integrals(edges, math.exp(log_slope))
        norm = float(shape @ precision @ shape)
        amplitude = float(shape @ weighted_y) / norm
        if amplitude <= 0.:
            return math.inf
        residual = yields-amplitude*shape
        return float(residual @ precision @ residual)

    span = float(edges[-1]-edges[0])
    grid = np.linspace(math.log(span*1.e-4), math.log(span*1.e3), 71)
    scores = [objective(point) for point in grid]
    best = int(np.argmin(scores))
    if best in (0, len(grid)-1) or not math.isfinite(scores[best]):
        raise ValueError("slope_not_identifiable")
    left, right = grid[best-1], grid[best+1]
    ratio = (math.sqrt(5.)-1.)/2.
    c, d = right-ratio*(right-left), left+ratio*(right-left)
    fc, fd = objective(c), objective(d)
    for _ in range(90):
        if fc < fd:
            right, d, fd = d, c, fc
            c = right-ratio*(right-left)
            fc = objective(c)
        else:
            left, c, fc = c, d, fd
            d = left+ratio*(right-left)
            fd = objective(d)
        if right-left < 1.e-10:
            break
    slope = math.exp((left+right)/2.)
    shape = exponential_integrals(edges, slope)
    amplitude = float(shape @ weighted_y)/float(shape @ precision @ shape)
    low, high = np.asarray(edges[:-1]), np.asarray(edges[1:])
    derivative = (np.exp(-low/slope)*(1.+low/slope)
                  - np.exp(-high/slope)*(1.+high/slope))
    jacobian = np.column_stack((shape, amplitude*derivative))
    try:
        parameter_covariance = np.linalg.inv(jacobian.T @ precision @ jacobian)
    except np.linalg.LinAlgError as exc:
        raise ValueError("slope_not_identifiable") from exc
    variance = float(parameter_covariance[1, 1])
    if not math.isfinite(variance) or variance <= 0.:
        raise ValueError("invalid_slope_variance")
    residual = yields-amplitude*shape
    chi2 = max(0., float(residual @ precision @ residual))
    ndof = len(yields)-2
    return slope, math.sqrt(variance), {
        "status": "ok", "chi2": chi2, "ndof": ndof,
        "poor_exponential_shape": chi2/ndof > 5.,
        "region_significances": region_significances,
        "amplitude": amplitude, "negative_bins_retained": int(np.sum(yields < 0.)),
        "fit_range_pt2": [float(edges[0]), float(edges[-1])],
    }


def pt2_slope_target_fit(spectra, target_weights, pt2_edges, cells,
                         covariance_proxies=None, minimum_region_significance=2.):
    edges = np.asarray(pt2_edges, dtype=float)
    bins = len(edges)-1
    if (bins < 8 or cells <= 0 or not np.all(np.isfinite(edges))
            or np.any(np.diff(edges) <= 0.) or edges[0] < 0.
            or not math.isfinite(minimum_region_significance)
            or minimum_region_significance < 0.):
        raise experimental.CampaignError("Invalid exponential fit definition")
    means, covariance = _moments(spectra, target_weights, bins, cells, covariance_proxies)
    values, errors, diagnostics, retained = [], [], [], []
    for yields, cov in zip(means, covariance):
        try:
            value, error, detail = _slope_fit(yields, cov, edges, minimum_region_significance)
            used = list(range(bins))
        except ValueError as exc:
            value, error, detail, used = None, None, {"status": str(exc)}, []
        values.append(value)
        errors.append(error)
        diagnostics.append(detail)
        retained.append(used)
    return {"edges": [float(i) for i in range(cells+1)], "values": values,
            "errors": errors, "retained_pt2_bins": retained,
            "fit_diagnostics": diagnostics,
            "covariance_model": ("event second moments across pT2 bins"
                                 if covariance_proxies is not None else "supplied diagonal errors")}
