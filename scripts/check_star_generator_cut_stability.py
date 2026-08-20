#!/usr/bin/env python3
"""Gate STAR 510-GeV production with independent generator-cut samples."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence


MEASUREMENT = "STAR_2022_I1949588"
PRIMARY_OBSERVABLES = ("inclusive", "dijet_A", "dijet_B", "dijet_C", "dijet_D")
NOMINAL_VARIATION = "nominal-jets-p000-u000-mu1-mpioff"


class StabilityError(RuntimeError):
    """Raised when a cut-scan input is incomplete or incompatible."""


def _load_json(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise StabilityError(f"Could not read {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise StabilityError(f"Expected a JSON object in {path}")
    return payload


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _finite_number(value: Any) -> float | None:
    if value is None:
        return None
    numeric = float(value)
    return numeric if math.isfinite(numeric) else None


def load_campaign(campaign_dir: Path, expected_cut: float) -> dict[str, Any]:
    manifest_path = campaign_dir / "manifest.json"
    summary_path = campaign_dir / "postprocess" / "summary.json"
    manifest = _load_json(manifest_path)
    summary = _load_json(summary_path)
    if manifest.get("measurement") != MEASUREMENT:
        raise StabilityError(f"{campaign_dir} is not a {MEASUREMENT} campaign")
    if manifest.get("status") != "complete":
        raise StabilityError(f"{campaign_dir} is not complete")
    jobs = manifest.get("jobs", [])
    if not jobs or any(job.get("status") != "success" for job in jobs):
        raise StabilityError(f"{campaign_dir} does not have an all-success shard matrix")
    configuration = manifest.get("configuration", {})
    cut = float(configuration.get("jet_kt_min_gev", -1.0))
    if not math.isclose(cut, expected_cut):
        raise StabilityError(
            f"{campaign_dir} records a {cut:g} GeV cut, expected {expected_cut:g}"
        )
    if int(configuration.get("lo_events", 0)) != 50_000_000:
        raise StabilityError(f"{campaign_dir} is not the 50M-event/helicity scan")
    if int(configuration.get("shards", 0)) != 100:
        raise StabilityError(f"{campaign_dir} does not contain 100 shards/helicity")
    if float(summary.get("jet_kt_min_gev", -1.0)) != expected_cut:
        raise StabilityError(f"{summary_path} has a mismatched generator cut")
    variations = summary.get("variations", {})
    prediction = variations.get(NOMINAL_VARIATION)
    if not isinstance(prediction, Mapping):
        raise StabilityError(f"{summary_path} has no central nominal prediction")
    source = manifest.get("runtime", {}).get("provenance", {}).get("source_control", {})
    return {
        "directory": str(campaign_dir.resolve()),
        "cut_gev": cut,
        "manifest": manifest,
        "summary": summary,
        "prediction": prediction,
        "source_commit": source.get("commit"),
        "manifest_sha256": _sha256(manifest_path),
        "summary_sha256": _sha256(summary_path),
        "initial_seeds": {int(job["initial_seed"]) for job in jobs},
    }


def _series(prediction: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    value = prediction.get(name)
    if not isinstance(value, Mapping):
        raise StabilityError(f"Central prediction is missing {name}")
    for key in ("edges", "values", "errors"):
        if not isinstance(value.get(key), list):
            raise StabilityError(f"{name} has no {key} array")
    if len(value["edges"]) != len(value["values"]) + 1:
        raise StabilityError(f"{name} has inconsistent bin edges")
    if len(value["values"]) != len(value["errors"]):
        raise StabilityError(f"{name} has inconsistent values/errors")
    return value


def compare_campaigns(
    reference: Mapping[str, Any],
    alternate: Mapping[str, Any],
    *,
    gate: bool,
    gate_bins_per_observable: int = 2,
    sigma_relative_floor: float = 0.02,
    asymmetry_absolute_floor: float = 5.0e-4,
    statistical_sigma: float = 3.0,
) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    failures: list[str] = []
    for observable in PRIMARY_OBSERVABLES:
        for quantity, series_name in (
            ("sigma_uu", f"SigmaUU_{observable}"),
            ("a_ll", observable),
        ):
            first = _series(reference["prediction"], series_name)
            second = _series(alternate["prediction"], series_name)
            if first["edges"] != second["edges"]:
                raise StabilityError(f"{series_name} binning differs across cut samples")
            finite_rank = 0
            for index, (a_value_raw, b_value_raw, a_error_raw, b_error_raw) in enumerate(
                zip(first["values"], second["values"], first["errors"], second["errors"])
            ):
                a_value = _finite_number(a_value_raw)
                b_value = _finite_number(b_value_raw)
                a_error = _finite_number(a_error_raw)
                b_error = _finite_number(b_error_raw)
                finite = None not in (a_value, b_value, a_error, b_error)
                is_gate_bin = False
                if finite:
                    is_gate_bin = finite_rank < gate_bins_per_observable
                    finite_rank += 1
                    difference = abs(float(b_value) - float(a_value))
                    combined_error = math.hypot(float(a_error), float(b_error))
                    if quantity == "sigma_uu":
                        scale = max(
                            0.5 * (abs(float(a_value)) + abs(float(b_value))),
                            1.0e-300,
                        )
                        metric = difference / scale
                        tolerance = max(
                            sigma_relative_floor,
                            statistical_sigma * combined_error / scale,
                        )
                    else:
                        metric = difference
                        tolerance = max(
                            asymmetry_absolute_floor,
                            statistical_sigma * combined_error,
                        )
                    passed = metric <= tolerance
                else:
                    difference = combined_error = metric = tolerance = None
                    passed = False
                gated = gate and is_gate_bin
                label = f"{observable}:{quantity}:bin{index + 1}"
                if gated and not passed:
                    failures.append(label)
                rows.append(
                    {
                        "label": label,
                        "observable": observable,
                        "quantity": quantity,
                        "bin": index + 1,
                        "low": first["edges"][index],
                        "high": first["edges"][index + 1],
                        "reference": a_value,
                        "alternate": b_value,
                        "combined_mc_error": combined_error,
                        "difference": difference,
                        "metric": metric,
                        "tolerance": tolerance,
                        "finite": finite,
                        "gated": gated,
                        "passed": passed,
                    }
                )
    gated_rows = [row for row in rows if row["gated"]]
    if gate and len(gated_rows) != 2 * gate_bins_per_observable * len(PRIMARY_OBSERVABLES):
        raise StabilityError("The gate does not contain two finite bins for every STAR observable")
    return {
        "reference_cut_gev": reference["cut_gev"],
        "alternate_cut_gev": alternate["cut_gev"],
        "gate": gate,
        "gate_bins_per_observable": gate_bins_per_observable,
        "passed": not failures,
        "failures": failures,
        "rows": rows,
    }


def build_report(campaigns: Mapping[float, Mapping[str, Any]]) -> dict[str, Any]:
    commits = {campaign.get("source_commit") for campaign in campaigns.values()}
    if None in commits or len(commits) != 1:
        raise StabilityError("Cut-scan campaigns do not share one pinned Git commit")
    cuts = sorted(campaigns)
    for first, second in zip(cuts, cuts[1:]):
        if campaigns[first]["initial_seeds"] & campaigns[second]["initial_seeds"]:
            raise StabilityError("Cut-scan campaigns reuse initial random seeds")
    gate = compare_campaigns(campaigns[3.0], campaigns[4.0], gate=True)
    stress = compare_campaigns(campaigns[4.0], campaigns[5.0], gate=False)
    return {
        "schema_version": 1,
        "measurement": MEASUREMENT,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "source_commit": next(iter(commits)),
        "criteria": {
            "gated_bins": "first two finite bins of every primary observable",
            "sigma_uu": "relative difference <= max(2%, 3 combined MC standard errors)",
            "a_ll": "absolute difference <= max(5e-4, 3 combined MC standard errors)",
            "five_gev": "reported as a non-gating stress test",
        },
        "inputs": {
            f"{cut:g}gev": {
                "directory": campaign["directory"],
                "manifest_sha256": campaign["manifest_sha256"],
                "summary_sha256": campaign["summary_sha256"],
                "source_commit": campaign["source_commit"],
            }
            for cut, campaign in sorted(campaigns.items())
        },
        "gate_passed": bool(gate["passed"]),
        "three_vs_four_gev": gate,
        "four_vs_five_gev_stress": stress,
    }


def markdown_report(report: Mapping[str, Any]) -> str:
    gate = report["three_vs_four_gev"]
    lines = [
        "# STAR 510 GeV generator-cut stability",
        "",
        f"Gate: **{'PASS' if report['gate_passed'] else 'FAIL'}**",
        "",
        "The 3 and 4 GeV samples gate the nominal 4 GeV production. The 5 GeV "
        "sample is reported only as a stress test.",
        "",
        "| Observable | Quantity | Bin | Metric | Tolerance | Result |",
        "|---|---|---:|---:|---:|---|",
    ]
    for row in gate["rows"]:
        if not row["gated"]:
            continue
        metric = "n/a" if row["metric"] is None else f"{row['metric']:.6g}"
        tolerance = "n/a" if row["tolerance"] is None else f"{row['tolerance']:.6g}"
        lines.append(
            f"| {row['observable']} | {row['quantity']} | {row['bin']} | "
            f"{metric} | {tolerance} | {'pass' if row['passed'] else 'FAIL'} |"
        )
    lines.extend(
        [
            "",
            f"Pinned source commit: `{report['source_commit']}`.",
            "",
        ]
    )
    return "\n".join(lines)


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign-3", type=Path, required=True)
    parser.add_argument("--campaign-4", type=Path, required=True)
    parser.add_argument("--campaign-5", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-markdown", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = make_parser().parse_args(argv)
    try:
        campaigns = {
            cut: load_campaign(path, cut)
            for cut, path in (
                (3.0, args.campaign_3),
                (4.0, args.campaign_4),
                (5.0, args.campaign_5),
            )
        }
        report = build_report(campaigns)
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_markdown.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        args.output_markdown.write_text(markdown_report(report), encoding="utf-8")
    except StabilityError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(
        f"STAR generator-cut gate: {'PASS' if report['gate_passed'] else 'FAIL'}; "
        f"report: {args.output_json}"
    )
    return 0 if report["gate_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
