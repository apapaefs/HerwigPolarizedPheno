"""Inspect Rivet's sample normalization, never the sign of an observable bin.

This reader intentionally supports the scalar metadata written by Rivet/YODA
2 and 3 only. It does not implement histogram reading or change any YODA file.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence


def read_sign_evidence(path: Path) -> dict[str, Any]:
    payload = path.read_bytes()
    text = (gzip.decompress(payload) if path.suffix == ".gz" else payload).decode("utf-8")
    wanted = {"/_XSEC", "/RAW/_XSEC", "/RAW/_EVTCOUNT"}
    blocks: dict[str, tuple[str, list[list[float]]]] = {}
    active = None
    data = False
    for line in text.splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "BEGIN":
            active = fields[2] if len(fields) == 3 and fields[2] in wanted else None
            if active is not None:
                if active in blocks:
                    raise ValueError(f"Duplicate normalization object {active} in {path}")
                blocks[active] = (fields[1], [])
            data = False
        elif fields[0] == "END":
            active = None
        elif active is not None:
            if fields[0] == "---":
                data = True
            elif data and not line.lstrip().startswith("#") and not fields[0].endswith(":"):
                try:
                    blocks[active][1].append([float(value) for value in fields])
                except ValueError as exc:
                    raise ValueError(f"Invalid normalization data in {path}: {active}") from exc

    def scalar(name: str, types: set[str]) -> list[float]:
        kind, rows = blocks.get(name, ("", []))
        if kind not in types or len(rows) != 1 or len(rows[0]) != 3:
            raise ValueError(f"Missing or unsupported {name} metadata in {path}")
        if not all(math.isfinite(value) for value in rows[0]):
            raise ValueError(f"Nonfinite {name} metadata in {path}")
        return rows[0]

    xsec_types = {"YODA_ESTIMATE0D_V3", "YODA_SCATTER1D_V2", "YODA_SCATTER1D"}
    xsec_path = "/_XSEC" if "/_XSEC" in blocks else "/RAW/_XSEC"
    xsec = scalar(xsec_path, xsec_types)[0]
    if "/_XSEC" in blocks and "/RAW/_XSEC" in blocks:
        raw_xsec = scalar("/RAW/_XSEC", xsec_types)[0]
        if not math.isclose(xsec, raw_xsec, rel_tol=2e-6, abs_tol=0.0):
            raise ValueError(f"Raw and normalized cross-section metadata disagree in {path}")
    sumw, sumw2, entries = scalar(
        "/RAW/_EVTCOUNT", {"YODA_COUNTER_V3", "YODA_COUNTER_V2", "YODA_COUNTER"}
    )
    if xsec == 0.0 or entries <= 0.0 or sumw == 0.0 or sumw2 <= 0.0:
        raise ValueError(f"Empty or zero-normalization sample cannot establish an NLO sign: {path}")
    # These campaigns are unweighted. Requiring unit weights avoids pretending
    # that sumW alone proves the absence of mixed positive/negative weights.
    unit_weights = (math.isclose(abs(sumw), entries, rel_tol=2e-6)
                    and math.isclose(sumw2, entries, rel_tol=2e-6))
    return {
        "sha256": hashlib.sha256(payload).hexdigest(), "cross_section": xsec,
        "sumw": sumw, "sumw2": sumw2, "entries": entries,
        "uniform_unit_weights": unit_weights,
    }


def audit_order_signs(
    jobs: Sequence[Mapping[str, Any]], directory: Path, *, mode: str = "auto",
    configured: Mapping[str, float] | None = None,
) -> dict[str, Any]:
    """Resolve one consistent stored-order convention across every NLO shard.

Rivet finalizes cross-section histograms with crossSection()/sumW(). For a
uniform signed event stream the resulting *normalized* bin sign therefore
follows crossSection(), not necessarily the sign of the raw event weights.
Explicit conventions permit non-unit weights but still require compatible
normalization metadata. Mixed stored conventions are rejected, not averaged.
"""
    if mode not in {"auto", "magnitudes", "signed"}:
        raise ValueError(f"Unknown NLO sign convention: {mode}")
    explicit = ({"POSNLO": 1.0, "NEGNLO": -1.0 if mode == "magnitudes" else 1.0}
                if mode != "auto" else None)
    if configured is not None:
        configured = {key: float(value) for key, value in configured.items()}
        if (set(configured) != {"POSNLO", "NEGNLO"}
                or any(not math.isfinite(value) or value == 0 for value in configured.values())):
            raise ValueError("NLO order_coefficients must explicitly retain both finite, nonzero contributions")
        if explicit is not None and configured != explicit:
            raise ValueError("--nlo-sign-convention conflicts with combination.order_coefficients")
        explicit = configured
    coefficients: dict[str, float] = {}
    counts = {"POSNLO": 0, "NEGNLO": 0}
    ranges: dict[str, list[float]] = {}
    digest = hashlib.sha256()
    weighted_files = 0
    root = directory.resolve()
    nlo_jobs = [job for job in jobs if str(job["order"]) in counts]
    for job in sorted(nlo_jobs, key=lambda item: str(item["output_yoda"])):
        order = str(job["order"])
        path = (root / str(job["output_yoda"])).resolve()
        if not path.is_relative_to(root):
            raise ValueError(f"NLO input path escapes campaign directory: {path}")
        evidence = read_sign_evidence(path)
        if not evidence["uniform_unit_weights"]:
            weighted_files += 1
            if explicit is None:
                raise ValueError(
                    f"Cannot infer a pure-sign, unit-weight stream in {path}; "
                    "verify its generation convention, then specify "
                    "--nlo-sign-convention magnitudes or signed explicitly"
                )
        xsec = evidence["cross_section"]
        if order == "POSNLO" and xsec < 0:
            raise ValueError(f"POSNLO has a negative normalized cross section: {path}")
        inferred = -1.0 if order == "NEGNLO" and xsec > 0 else 1.0
        coefficient = explicit[order] if explicit is not None else inferred
        if math.copysign(1., coefficient) != inferred:
            raise ValueError(f"Requested {order} coefficient {coefficient:g} contradicts stored normalization in {path}")
        if order in coefficients and coefficients[order] != coefficient:
            raise ValueError(f"Mixed stored {order} conventions across campaign shards; first disagreement: {path}")
        coefficients[order] = coefficient
        counts[order] += 1
        bounds = ranges.setdefault(order, [xsec, xsec])
        bounds[:] = [min(bounds[0], xsec), max(bounds[1], xsec)]
        digest.update(json.dumps({"path": str(job["output_yoda"]), "order": order,
                                  "evidence": evidence, "coefficient": coefficient},
                                 sort_keys=True).encode("utf-8"))
    if any(counts.values()) and not all(counts.values()):
        raise ValueError("The NLO sign audit requires both POSNLO and NEGNLO inputs")
    return {
        "requested_mode": mode, "resolution": "configured" if configured is not None else mode,
        "order_coefficients": coefficients, "files_checked": sum(counts.values()),
        "files_by_order": counts, "non_unit_weight_files": weighted_files,
        "cross_section_ranges": ranges, "inputs_sha256": digest.hexdigest(),
        "evidence": "Rivet cross-section normalization and RAW event counters; never observable-bin signs",
    }
