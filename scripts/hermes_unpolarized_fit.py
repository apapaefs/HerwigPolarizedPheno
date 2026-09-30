#!/usr/bin/env python3
"""Independent unpolarized Born weights for derived HERMES projections.

These weights use the GD11 unpolarized world-data fits. They do not use
Herwig, a spin-PDF fit, or a depolarization factor. The deuteron fit is per
nucleon. Experimental errors using these weights are conditional on the
central fit: the public, rounded fit covariance is retained but indefinite.
"""
from __future__ import annotations

from functools import lru_cache
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np

SNAPSHOT_PATH = Path(__file__).resolve().parents[1] / "data/experimental/HERMES_2007_I726689/unpolarized-fit.json"
_SNAPSHOT_BYTES = SNAPSHOT_PATH.read_bytes()
SNAPSHOT = json.loads(_SNAPSHOT_BYTES)
DEFAULT_CUTS = dict(SNAPSHOT["cuts"])
MODEL_METADATA = {
    "central": "GD11-P/D F2 with R1998; independent unpolarized world-data fits",
    "control": "ALLM97 proton F2; deuteron hybrid ALLM97-P * GD11-D / GD11-P",
    "alternate_R": "R1990 changes only the unpolarized cross-section conversion",
    "sources": SNAPSHOT["sources"],
    "snapshot_sha256": hashlib.sha256(_SNAPSHOT_BYTES).hexdigest(),
    "fit_validity": SNAPSHOT["fit_validity"],
    "error_policy": SNAPSHOT["error_policy"],
    "covariance_audit": {t: v["covariance_audit"] for t, v in SNAPSHOT["models"].items()},
}
_CONSTANTS = SNAPSHOT["constants"]
_E = DEFAULT_CUTS["beam_energy_GeV"]
_M = DEFAULT_CUTS["nucleon_mass_GeV"]


def _target(target: str) -> str:
    result = {"P": "P", "PROTON": "P", "D": "D", "DEUTERON": "D"}.get(str(target).upper())
    if result is None:
        raise ValueError("The unpolarized fit target must be P or D")
    return result


def _scalar_if_scalar(value):
    return float(value) if np.ndim(value) == 0 else value


def _allm_f2(x, q2, parameters):
    x, q2 = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(q2, dtype=float))
    if np.any(~np.isfinite(x)) or np.any(~np.isfinite(q2)) or np.any((x <= 0) | (x >= 1)) or np.any(q2 <= 0):
        raise ValueError("F2 requires finite 0 < x < 1 and Q2 > 0")
    p = parameters
    t = np.log(np.log((q2 + p[3]) / p[4]) / math.log(p[3] / p[4]))
    xp = (q2 + p[1]) / (q2 / x + p[1])
    xr = (q2 + p[2]) / (q2 / x + p[2])
    ap = p[5] + (p[5] - p[6]) * (1 / (1 + t**p[7]) - 1)
    cp = p[11] + (p[11] - p[12]) * (1 / (1 + t**p[13]) - 1)
    bp = p[8] + p[9] * t**p[10]
    ar = p[14] + p[15] * t**p[16]
    br = p[17] + p[18] * t**p[19]
    cr = p[20] + p[21] * t**p[22]
    return q2 / (q2 + p[0]) * (cp * xp**ap * (1 - x)**bp + cr * xr**ar * (1 - x)**br)


def f2(target: str, x, q2, model: str = "GD11"):
    """F2 per nucleon, using the pinned central fit or a named shape control."""
    target = _target(target)
    if model == "GD11":
        result = _allm_f2(x, q2, SNAPSHOT["models"][target]["parameters"])
    elif model == "ALLM97_HYBRID":
        result = _allm_f2(x, q2, SNAPSHOT["ALLM97_P_parameters"])
        if target == "D":
            result = result * _allm_f2(x, q2, SNAPSHOT["models"]["D"]["parameters"]) / _allm_f2(x, q2, SNAPSHOT["models"]["P"]["parameters"])
    else:
        raise ValueError(f"Unknown unpolarized fit model: {model}")
    return _scalar_if_scalar(result)


def longitudinal_ratio(x, q2, r_model: str = "R1998"):
    """Empirical R=σL/σT used only to convert unpolarized F2 to σUU."""
    x, q2 = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(q2, dtype=float))
    if np.any(~np.isfinite(x)) or np.any(~np.isfinite(q2)) or np.any((x <= 0) | (x >= 1)) or np.any(q2 <= .04):
        raise ValueError("R requires finite 0 < x < 1 and Q2 > 0.04")
    theta = 1 + 12 * q2 / (q2 + 1) * .125**2 / (.125**2 + x*x)
    logq = np.log(q2 / .04)
    if r_model == "R1990":
        a, b, c = SNAPSHOT["R1990_coefficients"]
        result = a / logq * theta + b / q2 + c / (q2*q2 + .09)
    elif r_model == "R1998":
        a, b, c = (SNAPSHOT["R1998_coefficients"][key] for key in ("a", "b", "c"))
        ra = a[0] / logq * theta + a[1] / (q2**4 + a[2]**4)**.25 * (1 + a[3]*x + a[4]*x*x) * x**a[5]
        rb = b[0] / logq * theta + (b[1]/q2 + b[2]/(q2*q2 + .09)) * (1 + b[3]*x + b[4]*x*x) * x**b[5]
        threshold = c[3]*x + c[4]*x*x + c[5]*x**3
        rc = c[0] / logq * theta + c[1] / np.sqrt((q2-threshold)**2 + c[2]**2)
        result = (ra + rb + rc) / 3
    else:
        raise ValueError(f"Unknown R model: {r_model}")
    return _scalar_if_scalar(result)


def differential_cross_section(target: str, x, q2, model: str = "GD11", r_model: str = "R1998"):
    """Born d²σUU/(dx dQ²), in pb/GeV², before fiducial selection.

    HERMES2011 Eq2.6 retains the finite nucleon-mass kinematics and neglects
    the positron mass. There is no longitudinal-spin conversion here.
    """
    x, q2 = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(q2, dtype=float))
    y = q2 / (2 * _M * _E * x)
    r = longitudinal_ratio(x, q2, r_model)
    bracket = 1 - y - q2 / (4 * _E**2) + (y*y + q2 / _E**2) / (2 * (1 + r))
    result = 4 * math.pi * _CONSTANTS["alpha_em"]**2 / q2**2 * f2(target, x, q2, model) / x * bracket * _CONSTANTS["GeV_minus2_to_pb"]
    return _scalar_if_scalar(result)


def _boundary_functions(q2low, q2high, w2_min):
    # Each Q2 boundary is (n0+n1*x)/(d0+d1*x). All intersections are analytic.
    def theta_boundary(theta):
        s = math.sin(theta/2)**2
        return (0., 4*_E*_E*s, 2*_E/_M*s, 1.)
    lower = [(q2low, 0., 1., 0.), (0., 2*_M*_E*DEFAULT_CUTS["y_min"], 1., 0.),
             (0., w2_min-_M*_M, 1., -1.), theta_boundary(DEFAULT_CUTS["theta_min_rad"])]
    upper = [(q2high, 0., 1., 0.), (0., 2*_M*_E*DEFAULT_CUTS["y_max"], 1., 0.),
             theta_boundary(DEFAULT_CUTS["theta_max_rad"])]
    return lower, upper


def _evaluate_boundary(boundary, x):
    n0, n1, d0, d1 = boundary
    return (n0+n1*x)/(d0+d1*x)


def _breakpoints(xlow, xhigh, boundaries):
    roots = [xlow, xhigh]
    for left, right in itertools.combinations(boundaries, 2):
        n0, n1, d0, d1 = left
        m0, m1, e0, e1 = right
        c0 = n0*e0-m0*d0
        c1 = n0*e1+n1*e0-m0*d1-m1*d0
        c2 = n1*e1-m1*d1
        if c2 == 0:
            candidates = [] if c1 == 0 else [-c0/c1]
        else:
            disc = c1*c1-4*c2*c0
            candidates = [] if disc < 0 else [(-c1-math.sqrt(disc))/(2*c2), (-c1+math.sqrt(disc))/(2*c2)]
        roots.extend(v for v in candidates if xlow < v < xhigh)
    return sorted(set(roots))


def q2_limits(x, q2low: float = 1., q2high: float = 20., w2_min: float = 3.24):
    """Analytic actual HERMES fiducial Q2 limits at Bjorken x."""
    x = np.asarray(x, dtype=float)
    lower, upper = _boundary_functions(q2low, q2high, w2_min)
    low = np.maximum.reduce([_evaluate_boundary(b, x) for b in lower])
    high = np.minimum.reduce([_evaluate_boundary(b, x) for b in upper])
    return _scalar_if_scalar(low), _scalar_if_scalar(high)


@lru_cache(maxsize=32)
def _gauss(order):
    if not isinstance(order, int) or order < 4:
        raise ValueError("quadrature_order must be an integer >=4")
    return np.polynomial.legendre.leggauss(order)


def integrate_cross_section(target: str, xlow: float, xhigh: float, q2low: float,
                            q2high: float, quadrature_order: int = 32,
                            model: str = "GD11", r_model: str = "R1998",
                            w2_min: float = 3.24) -> float:
    """Integrate independent σUU over a cell intersection, returning pb.

    Applies the actual v4 Q2/y/W2/theta/x cuts. Gauss integration splits x
    at every analytic boundary crossing and integrates in log(Q2), so empty
    or narrow clipped cells are handled without rectangular cut sampling.
    w2_min=4 is a fit-domain weight sensitivity control, not a new MC cut.
    """
    target = _target(target)
    if model not in {"GD11", "ALLM97_HYBRID"} or r_model not in {"R1998", "R1990"}:
        raise ValueError("Unknown unpolarized model")
    if not all(math.isfinite(v) for v in (xlow, xhigh, q2low, q2high, w2_min)):
        raise ValueError("Cell boundaries must be finite")
    if xhigh < xlow or q2high < q2low or w2_min < DEFAULT_CUTS["w2_min_GeV2"]:
        raise ValueError("Invalid cell boundaries or a W2 cut looser than the campaign")
    nodes, weights = _gauss(quadrature_order)
    xlow, xhigh = max(xlow, DEFAULT_CUTS["x_min"]), min(xhigh, DEFAULT_CUTS["x_max"])
    q2low, q2high = max(q2low, DEFAULT_CUTS["q2_min_GeV2"]), min(q2high, DEFAULT_CUTS["q2_max_GeV2"])
    if xhigh <= xlow or q2high <= q2low:
        return 0.
    lower, upper = _boundary_functions(q2low, q2high, w2_min)
    points = _breakpoints(xlow, xhigh, lower+upper)
    result = 0.
    for xa, xb in zip(points[:-1], points[1:]):
        x = (xa+xb)/2 + (xb-xa)/2 * nodes
        low, high = q2_limits(x, q2low, q2high, w2_min)
        # The boundary splitting guarantees a segment is wholly open or closed.
        if np.all(high <= low):
            continue
        if np.any(high <= low):
            raise RuntimeError("Analytic fiducial boundary segmentation failed")
        logs = np.log(low), np.log(high)
        q2 = np.exp((logs[0][:, None]+logs[1][:, None])/2 + (logs[1]-logs[0])[:, None]/2 * nodes)
        inner = np.sum(differential_cross_section(target, x[:, None], q2, model, r_model) * q2 * weights, axis=1) * (logs[1]-logs[0])/2
        result += float(np.dot(weights, inner)) * (xb-xa)/2
    return result
