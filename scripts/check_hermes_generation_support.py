#!/usr/bin/env python3
"""Validate HERMES v5 support using completed, actual Herwig/Rivet shards.

Only normalized shard moments are combined: event-count-weighted shard
averages, POSNLO minus NEGNLO, then the four-helicity UU average. This does
not merge raw histograms or use GD11 as a generator prediction.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass, asdict
import gzip
import hashlib
import json
import math
from pathlib import Path
import sys

import run_experimental_campaign as campaign

ANALYSIS = "HERMES_2007_I726689"
HELICITIES = ("PP", "PM", "MP", "MM")
ORDERS = {"POSNLO": 1., "NEGNLO": -1.}
SELECTIONS = {"rectangle": "", "ring_control": "RingControl_"}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_yoda_objects(path):
    """Read only the relevant v3 blocks, retaining the full-precision scale."""
    opener = gzip.open if path.suffix == ".gz" else open
    objects, block = {}, None
    with opener(path, "rt") as stream:
        for line in stream:
            line = line.strip()
            if line.startswith("BEGIN "):
                _, kind, name = line.split(maxsplit=2)
                keep = name.split("/")[-1] in {
                    prefix + suffix for prefix in SELECTIONS.values()
                    for suffix in ("BornSigmaQ2_X19", "Accepted_W2Fine_Q2GT1")
                }
                block = {"kind": kind, "path": name, "rows": []} if keep else None
            elif line.startswith("END "):
                if block is not None:
                    objects[block["path"]] = block
                block = None
            elif block is not None:
                if line.startswith("Edges(A1):"):
                    block["edges"] = json.loads(line.split(":", 1)[1])
                elif line.startswith("ScaledBy:"):
                    block["scale"] = float(line.split(":", 1)[1])
                elif line and line != "---" and line[0] in "-+.0123456789n" and ":" not in line:
                    block["rows"].append([float(value) if value != "---" else math.nan
                                          for value in line.split()])
    return objects


@dataclass
class Moments:
    sum_w: float = 0.
    sum_wq: float = 0.
    variance_w: float = 0.
    covariance_w_wq: float = 0.
    variance_wq: float = 0.
    raw_entries: float = 0.

    def add(self, other, coefficient=1.):
        self.sum_w += coefficient*other.sum_w
        self.sum_wq += coefficient*other.sum_wq
        self.variance_w += coefficient**2*other.variance_w
        self.covariance_w_wq += coefficient**2*other.covariance_w_wq
        self.variance_wq += coefficient**2*other.variance_wq
        self.raw_entries += other.raw_entries
        return self

    def report(self):
        mean = self.sum_wq/self.sum_w if self.sum_w else None
        variance = ((self.variance_wq - 2*mean*self.covariance_w_wq
                     + mean**2*self.variance_w)/self.sum_w**2) if mean is not None else None
        # Production YODA stores these raw moments with six significant
        # figures. Their cancellation can give a tiny negative variance for
        # one-entry bins; bound this by the precision of the source terms.
        rounding_bound = (3e-6*(abs(self.variance_wq)
                          + 2*abs(mean*self.covariance_w_wq)
                          + mean**2*abs(self.variance_w))/self.sum_w**2
                          if mean is not None else 0.)
        if variance is not None and variance < -rounding_bound:
            raise ValueError(f"Invalid propagated mean variance: {variance}")
        return dict(asdict(self), mean_q2=mean,
                    mean_variance_roundoff_clamped=variance is not None and variance < 0.,
                    mean_q2_mc_error=math.sqrt(max(0., variance)) if variance is not None else None,
                    effective_entries=self.sum_w**2/self.variance_w if self.variance_w else 0.)


def normalized_bin(objects, name, index):
    raw = objects.get(f"/RAW/{ANALYSIS}/{name}")
    estimate = objects.get(f"/{ANALYSIS}/{name}")
    if raw is None or estimate is None:
        raise ValueError(f"Missing raw or normalized {name}")
    if raw["kind"] != "YODA_HISTO1D_V3" or estimate["kind"] != "YODA_ESTIMATE1D_V3":
        raise ValueError("Support validation currently requires the production YODA v3 representation")
    if raw["edges"] != estimate["edges"]:
        raise ValueError(f"Raw/normalized binning mismatch: {name}")
    # YODA v3 includes underflow and overflow rows.
    row = raw["rows"][index+1]
    sw, sw2, swq, swq2, entries = row
    scale = estimate.get("scale")
    if scale is None or not math.isfinite(scale):
        raise ValueError(f"Missing finite ScaledBy annotation: {name}")
    width = raw["edges"][index+1] - raw["edges"][index]
    density = estimate["rows"][index+1][0]
    if not math.isclose(density, scale*sw/width, rel_tol=3e-6, abs_tol=1e-12):
        raise ValueError(f"Finalized normalization does not match raw moments: {name}")
    # YODA stores sum(w q) and sum(w q^2), not sum(w^2 q) or sum(w^2 q^2).
    # In a constant-weight stream these follow exactly. Verify this property
    # from the recorded moments before using it, including negative weights.
    if entries and not math.isclose(sw*sw, entries*sw2, rel_tol=3e-6, abs_tol=1e-12):
        raise ValueError(f"Nonconstant event weights prevent exact mean-error recovery: {name}, bin {index}")
    weight = sw/entries if entries else 0.
    return Moments(scale*sw, scale*swq, scale**2*sw2,
                   scale**2*weight*swq, scale**2*weight*swq2, entries)


def inspect_campaign(directory, minimum_effective_entries=100., maximum_ring_mean=6.8):
    directory = directory.resolve()
    manifest_path = directory/"manifest.json"
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("configuration", {}).get("measurement") != ANALYSIS:
        raise ValueError("Expected a HERMES_2007_I726689 campaign")
    groups = campaign._logical_job_groups(manifest)
    expected_shards = int(manifest["configuration"]["shards"])
    report = {
        "schema_version": 1, "campaign": str(directory),
        "validation_scope": "restricted_high_x_probe" if manifest.get("support_probe") else "full_generation_window",
        "support_probe": manifest.get("support_probe"),
        "manifest_sha256": sha256(manifest_path), "configuration": manifest["configuration"],
        "runtime": manifest.get("runtime", {}),
        "normalization": "finalized ScaledBy times raw moments; event-count-weighted shard average; (POSNLO-NEGNLO); sum helicities/4",
        "uncertainty": "independent Poisson moments at fixed shard cross sections; constant event weights checked per bin",
        "reference": {"v4_proton_ring_cell43_mean_q2": 7.07321435,
                      "gd11_full_window_ring_mean_q2_illustration": 6.30570,
                      "note": "The GD11 value is an illustration, not an exact Herwig target or a pass criterion."},
        "thresholds": {"minimum_cell43_effective_entries": minimum_effective_entries,
                       "maximum_ring_cell43_mean_q2": maximum_ring_mean},
        "inputs": [], "targets": {}, "failures": [],
    }
    for target in ("P", "N"):
        target_results = {}
        totals = {selection: {"cell43": Moments(), "low_w2": [Moments() for _ in range(7)]}
                  for selection in SELECTIONS}
        details = {selection: {} for selection in SELECTIONS}
        for helicity in HELICITIES:
            for order, order_sign in ORDERS.items():
                key = ("nominal", 0, 0, 1., target, helicity, order)
                jobs = groups.get(key, [])
                if len(jobs) != expected_shards:
                    raise ValueError(f"Incomplete shard matrix for {target}/{helicity}/{order}")
                if any(job.get("status") != "success" for job in jobs):
                    raise ValueError(f"Unfinished jobs for {target}/{helicity}/{order}")
                events = sum(int(job["events"]) for job in jobs)
                if events <= 0:
                    raise ValueError("Zero generated events")
                logical = {selection: {"cell43": Moments(), "low_w2": [Moments() for _ in range(7)]}
                           for selection in SELECTIONS}
                for job in jobs:
                    path = directory/str(job["output_yoda"])
                    objects = read_yoda_objects(path)
                    fraction = int(job["events"])/events
                    report["inputs"].append({"job": job["id"], "path": str(job["output_yoda"]),
                                             "events": int(job["events"]), "sha256": sha256(path)})
                    for selection, prefix in SELECTIONS.items():
                        cell = normalized_bin(objects, prefix+"BornSigmaQ2_X19", 0)
                        logical[selection]["cell43"].add(cell, fraction)
                        edges = objects[f"/RAW/{ANALYSIS}/{prefix}Accepted_W2Fine_Q2GT1"]["edges"]
                        if len(edges) != 29 or any(not math.isclose(value, 3.24+.1*i, abs_tol=1e-8)
                                                   for i, value in enumerate(edges)):
                            raise ValueError("Unexpected low-W2 diagnostic binning")
                        for i in range(7):
                            logical[selection]["low_w2"][i].add(
                                normalized_bin(objects, prefix+"Accepted_W2Fine_Q2GT1", i), fraction)
                for selection in SELECTIONS:
                    item = logical[selection]
                    totals[selection]["cell43"].add(item["cell43"], order_sign/4.)
                    for i, moment in enumerate(item["low_w2"]):
                        totals[selection]["low_w2"][i].add(moment, order_sign/4.)
                    details[selection][f"{helicity}/{order}"] = {
                        "events": events, "cell43": item["cell43"].report(),
                        "low_w2_cross_sections": [value.sum_w for value in item["low_w2"]],
                        "low_w2_raw_entries": [value.raw_entries for value in item["low_w2"]],
                    }
        for selection in SELECTIONS:
            total = totals[selection]
            cell = total["cell43"].report()
            low = [{"low": 3.24+.1*i, "high": 3.24+.1*(i+1),
                    "cross_section_pb": item.sum_w, "mc_error_pb": math.sqrt(item.variance_w),
                    "raw_entries": item.raw_entries} for i, item in enumerate(total["low_w2"])]
            if any(item["cross_section_pb"] <= 0 or item["raw_entries"] <= 0 for item in low):
                report["failures"].append(f"{target}/{selection}: low-W2 bins are not all populated with positive normalized UU cross section")
            if cell["effective_entries"] < minimum_effective_entries:
                report["failures"].append(f"{target}/{selection}: insufficient cell-43 effective entries ({cell['effective_entries']:.1f})")
            if selection == "ring_control" and (cell["mean_q2"] is None or cell["mean_q2"] >= maximum_ring_mean):
                report["failures"].append(f"{target}/{selection}: cell-43 mean Q2 has not moved below {maximum_ring_mean}")
            target_results[selection] = {"cell43": cell, "low_w2_bins": low,
                                         "contributions": details[selection]}
        report["targets"][target] = target_results
    report["passed"] = not report["failures"]
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--minimum-effective-entries", type=float, default=100.)
    parser.add_argument("--maximum-ring-mean", type=float, default=6.8)
    args = parser.parse_args()
    try:
        report = inspect_campaign(args.campaign, args.minimum_effective_entries, args.maximum_ring_mean)
    except (OSError, ValueError, KeyError, campaign.CampaignError) as exc:
        report = {"passed": False, "campaign": str(args.campaign), "failures": [str(exc)]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"passed": report["passed"], "failures": report["failures"],
                      "report": str(args.output)}))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
