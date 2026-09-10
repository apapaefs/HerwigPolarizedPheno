#!/usr/bin/env python3
"""Audit a completed harder-jet pilot, optionally against unbiased spectra.

This checks compensated-sampling normalization, not full/LHE physics closure.
Reference and pilot must be independent samples with the same physics settings.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path

import run_experimental_campaign as experimental


def normalization_comparison(pilot, reference):
    rows = []
    for variation, outputs in pilot["variations"].items():
        other = reference["variations"][variation]
        for jet in (1, 2):
            for helicity in ("UU", "PP", "PM", "MP", "MM"):
                name = f"Sigma{helicity}_jet{jet}_pt"
                a, b = outputs[name], other[name]
                if a["edges"] != b["edges"]:
                    raise ValueError(f"Reference edges differ for {name}")
                for lower, upper in ((20,30), (30,45)):
                    i = a["edges"].index(lower)
                    if a["edges"][i+1] != upper:
                        raise ValueError(f"Reference does not contain [{lower},{upper})")
                    x, y = a["values"][i], b["values"][i]
                    ex, ey = a["errors"][i], b["errors"][i]
                    if x is None or y is None or y <= 0:
                        raise ValueError(f"Missing positive normalization input {variation}/{name}")
                    ratio = x/y
                    error = math.hypot(ex/y, x*ey/y**2)
                    rows.append({"variation": variation, "observable": name,
                                 "pt_window_gev": [lower,upper], "ratio": ratio,
                                 "independent_error": error,
                                 "relative_95_bound": abs(ratio-1)+1.96*error,
                                 "pull": (x-y)/math.hypot(ex,ey)})
    return rows


def audit(directory: Path, reference_path: Path | None = None, require_html=True):
    manifest = json.loads((directory/"manifest.json").read_text())
    summary = json.loads((directory/"postprocess/summary.json").read_text())
    failures = []
    if manifest["measurement"] != "MC_POLJETSHAPES_HARD":
        failures.append("Not a MC_POLJETSHAPES_HARD campaign")
    for job in manifest["jobs"]:
        if job["status"] != "success":
            failures.append("Incomplete shard: "+job["id"])
            continue
        content = (directory/job["output_yoda"]).read_text()
        for name in ("MC_POLDIJETS", "MC_POLJETSHAPES_HARD"):
            if "/"+name+":" not in content and "/"+name+"/" not in content:
                failures.append("Missing raw namespace "+name+": "+job["id"])
    closures = []
    primary_yields = {}
    for variation, outputs in summary["variations"].items():
        primary_yields[variation] = {}
        for name, item in outputs.items():
            if name.startswith("Shape") and all(v is not None for v in item["values"]):
                integral = sum(v*(hi-lo) for v,lo,hi in zip(item["values"], item["edges"], item["edges"][1:]))
                if integral:
                    closures.append(abs(integral-1))
            if name.startswith(("R32_", "R43_", "ThirdJetVeto_")):
                if any(v is not None and not -1.e-8 <= v <= 1+1.e-8 for v in item["values"]):
                    failures.append("Invalid rate "+variation+"/"+name)
        for window in ("pt20_30", "pt30_45"):
            for jet in (1,2):
                name = f"SigmaUU_{window}_dpsi12_j{jet}_kt1"
                item = outputs[name]
                total = sum(v*(hi-lo) for v,lo,hi in zip(item["values"], item["edges"], item["edges"][1:]))
                primary_yields[variation][name] = total
                if total <= 0:
                    failures.append("Empty baseline primary/secondary observable "+variation+"/"+name)
            for helicity in ("UU", "PP", "PM", "MP", "MM"):
                r32 = outputs[f"R32_{helicity}_{window}"]["values"]
                veto = outputs[f"ThirdJetVeto_{helicity}_{window}"]["values"]
                if any(a is not None and b is not None and abs(a+b-1) > 1e-8 for a,b in zip(r32,veto)):
                    failures.append("R32+veto identity failed")
                if any(a is not None and b is not None and b > a+1e-7 for a,b in zip(r32,r32[1:])):
                    failures.append("Non-monotonic R32")
            for observable in ("ge3_threshold", "ge4_threshold", "pt31_cumulative_tail", "pt41_cumulative_tail"):
                item = outputs[f"SigmaUU_{window}_{observable}"]
                integrals = [v*(hi-lo) for v,lo,hi in zip(item["values"], item["edges"], item["edges"][1:])]
                if any(b > a+max(1e-8, abs(a)*2e-6) for a,b in zip(integrals,integrals[1:])):
                    failures.append("Non-monotonic cumulative scan "+variation+"/"+observable)
            # Event-window selection must close to the companion jet-1 bin.
            lo, hi = (20,30) if window == "pt20_30" else (30,45)
            spectrum = outputs["SigmaUU_jet1_pt"]
            index = spectrum["edges"].index(lo)
            expected = spectrum["values"][index]*(hi-lo)
            actual = outputs[f"SigmaUU_{window}_dijet_rate"]["values"][0]
            if not math.isclose(actual, expected, rel_tol=2e-6, abs_tol=1e-6):
                failures.append("Leading-window rate does not close to companion jet-1 bin")
    if not closures or max(closures) > 1e-8:
        failures.append("Normalized-shape closure failed")
    if require_html:
        for relative in ("plots/html/index.html", "plots/html/focus/index.html"):
            if not (directory/relative).is_file():
                failures.append("Missing "+relative)
    result = {"tag": manifest["tag"], "successful_shards": sum(j["status"] == "success" for j in manifest["jobs"]),
              "requested_events": sum(j["events"] for j in manifest["jobs"]),
              "normalized_shapes_checked": len(closures), "maximum_shape_closure_error": max(closures, default=None),
              "primary_selected_cross_sections_pb": primary_yields,
              "structural_failures": failures, "structural_checks_pass": not failures,
              "sampling": summary.get("sampling", {}),
              "statistics_projection": summary["statistics_projection"]}
    if reference_path:
        reference = json.loads(reference_path.read_text())
        if summary["shower_spin_policy"] != reference["shower_spin_policy"]:
            raise ValueError("Reference has different shower switches")
        for family, instances in summary["analysis_instances"].items():
            if (instances["MC_POLDIJETS"] !=
                    reference["analysis_instances"][family]["MC_POLDIJETS"]):
                raise ValueError("Reference has different companion-jet cuts")
        if summary["jet_kt_min_gev"] != reference["jet_kt_min_gev"]:
            raise ValueError("Reference has a different generation floor")
        rows = normalization_comparison(summary, reference)
        uu = [r for r in rows if r["observable"].startswith("SigmaUU")]
        result["sampling_normalization_check"] = {
            "reference": str(reference_path), "reference_tag": reference["tag"],
            "rows": rows, "maximum_UU_relative_95_bound": max(r["relative_95_bound"] for r in uu),
            "maximum_absolute_UU_pull": max(abs(r["pull"]) for r in uu),
            "interpretation": "Pointwise independent-sample sampling check only; not a simultaneous confidence region or proof of angular LHE equivalence. Correlations between bins/jets and UU/helicities are not fitted as independent points."}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign-dir", type=Path, required=True)
    parser.add_argument("--reference-summary", type=Path)
    parser.add_argument("--without-html", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = audit(args.campaign_dir, args.reference_summary, not args.without_html)
    if args.output:
        experimental.atomic_write_json(args.output, result)
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0 if result["structural_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
