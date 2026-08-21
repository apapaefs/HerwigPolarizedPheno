#!/usr/bin/env python3
"""Shared primary-bin policy for STAR 510-GeV comparisons.

The published inclusive display coordinates are not the analysis bin edges.
Selections therefore use one-based analysis-bin indices and independently
validate the declared lower analysis edge.
"""

from __future__ import annotations

import hashlib
import json
import math
from typing import Any, Mapping, Sequence


class ComparisonPolicyError(ValueError):
    """Raised when a configured comparison policy is inconsistent."""


def policy_payload(measurement: Mapping[str, Any]) -> dict[str, Any]:
    value = measurement.get("comparison_policy", {})
    if not isinstance(value, Mapping):
        raise ComparisonPolicyError("comparison_policy must be an object")
    return json.loads(json.dumps(value, sort_keys=True))


def policy_sha256(measurement: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        policy_payload(measurement),
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def primary_bin_mask(
    measurement: Mapping[str, Any],
    snapshot: Mapping[str, Any],
    observable: str,
    size: int,
) -> list[bool]:
    """Return the primary-comparison mask for *observable*.

    Unconfigured observables retain every bin. Configured masks are validated
    against the analysis edges in the vendored reference snapshot.
    """

    policy = policy_payload(measurement)
    specifications = policy.get("primary_bin_masks", {})
    if not isinstance(specifications, Mapping):
        raise ComparisonPolicyError(
            "comparison_policy.primary_bin_masks must be an object"
        )
    specification = specifications.get(observable)
    if specification is None:
        return [True] * size
    if not isinstance(specification, Mapping):
        raise ComparisonPolicyError(
            f"primary-bin policy for {observable} must be an object"
        )
    try:
        first_bin = int(specification["first_bin"])
        expected_low = float(specification["minimum_analysis_low_gev"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ComparisonPolicyError(
            f"primary-bin policy for {observable} needs first_bin and "
            "minimum_analysis_low_gev"
        ) from exc
    if first_bin < 1 or first_bin > size:
        raise ComparisonPolicyError(
            f"primary first bin {first_bin} for {observable} is outside 1..{size}"
        )
    datasets = snapshot.get("datasets", {})
    dataset = datasets.get(observable) if isinstance(datasets, Mapping) else None
    if not isinstance(dataset, Mapping):
        raise ComparisonPolicyError(
            f"reference snapshot has no dataset {observable}"
        )
    edges = dataset.get("bin_edges")
    if not isinstance(edges, Sequence) or isinstance(edges, (str, bytes)):
        raise ComparisonPolicyError(
            f"reference snapshot has no analysis edges for {observable}"
        )
    if len(edges) != size + 1:
        raise ComparisonPolicyError(
            f"{observable} has {size} prediction bins but {len(edges)-1} "
            "reference analysis bins"
        )
    actual_low = float(edges[first_bin - 1])
    if not math.isclose(actual_low, expected_low, rel_tol=0.0, abs_tol=1.0e-12):
        raise ComparisonPolicyError(
            f"{observable} primary bin {first_bin} starts at {actual_low:g} GeV, "
            f"not the configured {expected_low:g} GeV"
        )
    return [index >= first_bin for index in range(1, size + 1)]


def masked_copy(
    prediction: Mapping[str, Any], mask: Sequence[bool]
) -> dict[str, Any]:
    """Copy a prediction and replace excluded values/errors with nulls."""

    if len(prediction.get("values", [])) != len(mask):
        raise ComparisonPolicyError("prediction and primary-bin mask differ in size")
    result = json.loads(json.dumps(prediction))
    for field in ("values", "errors"):
        values = result.get(field)
        if isinstance(values, list):
            result[field] = [
                value if keep else None
                for value, keep in zip(values, mask)
            ]
    return result


def comparison_roles(mask: Sequence[bool]) -> list[str]:
    return ["primary" if keep else "diagnostic_only" for keep in mask]
