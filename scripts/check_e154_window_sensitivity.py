#!/usr/bin/env python3
"""Report nested-window shifts after E154 postprocessing; this is not a chi2 test."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def compare(summary_csv: Path) -> dict:
    with summary_csv.open(newline="", encoding="utf-8") as stream:
        rows = {(row["selection"], int(row["bin"])): row for row in csv.DictReader(stream)}
    snapshot = json.loads((ROOT/"data/experimental/E154_1997_I443170/reference.json").read_text())
    records = []
    for dataset in snapshot["datasets"]:
        target = dataset["target"]
        for point in dataset["points"]:
            index = point["bin"]
            nominal = rows[(f"{target}_W050", index)]
            for window in ("W025", "W100"):
                alternative = rows[(f"{target}_{window}", index)]
                valid = bool(nominal["a1"] and alternative["a1"])
                delta = float(alternative["a1"])-float(nominal["a1"]) if valid else None
                data_error = math.hypot(point["stat"], point["systematic_combined"])
                records.append({
                    "target": target, "bin": index, "x_mean": point["x_mean"],
                    "window": window, "nominal_a1": float(nominal["a1"]) if nominal["a1"] else None,
                    "alternative_a1": float(alternative["a1"]) if alternative["a1"] else None,
                    "delta_a1": delta,
                    "delta_over_data_total_error": delta/data_error if delta is not None else None,
                    "nominal_mc_stat": float(nominal["a1_stat"]) if nominal["a1_stat"] else None,
                    "alternative_mc_stat": float(alternative["a1_stat"]) if alternative["a1_stat"] else None,
                })
    return {
        "measurement": "E154_1997_I443170", "source": str(summary_csv),
        "interpretation": "Changes relative to the 5-percent window, scaled by the experimental total error for context; not a significance or goodness-of-fit statistic.",
        "correlation": "Nested windows share events. MC errors are reported separately; no independent-error significance is calculated.",
        "records": records,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary_csv", type=Path, help="Campaign postprocess/summary.csv")
    args = parser.parse_args()
    print(json.dumps(compare(args.summary_csv), indent=2, allow_nan=False))
