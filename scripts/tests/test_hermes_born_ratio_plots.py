#!/usr/bin/env python3
"""Additive Born ratios preserve signed data, support masks and original plots."""
from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"scripts"))
import hermes_born_cell_plots as born
import hermes_born_ratio_plots as plots


def fixture():
    reference = json.loads(born.REFERENCE.read_text())
    rows = []
    for dataset in reference["datasets"]:
        for point in dataset["points"]:
            supported = point["supported_prediction"]
            rows.append({"selection": dataset["selection"], "bin": point["bin"],
                         "supported": supported, "bin_low": point["q2_low"],
                         "bin_high": point["q2_high"],
                         "a_parallel": .2 if supported else None,
                         "a_parallel_stat": .03 if supported else None})
    return {"measurement": born.MEASUREMENT, "tag": "ratio-test", "bins": rows}, reference


class BornRatioSnapshotTests(unittest.TestCase):
    def setUp(self):
        self.summary, self.reference = fixture()
        self.base = born.build_snapshot(self.summary, self.reference)
        self.snapshot = plots.with_ratio_metadata(self.base)

    def test_preserves_the_original_cells_and_every_published_coordinate(self):
        original = copy.deepcopy(self.base)
        summary, reference = copy.deepcopy(self.summary), copy.deepcopy(self.reference)
        result = plots.build_snapshot(self.summary, self.reference)
        self.assertEqual(self.base, original)
        self.assertEqual(summary, self.summary)
        self.assertEqual(reference, self.reference)
        for target in ("P", "D"):
            self.assertEqual(result["targets"][target]["cells"], original["targets"][target]["cells"])
            count = 0
            for panel in result["targets"][target]["panels"]:
                for cell, row in zip(panel["cells"], panel["ratio_rows"]):
                    count += 1
                    self.assertEqual([row[k] for k in ("low", "high", "marker")],
                                     [cell[k] for k in ("q2_low", "q2_high", "q2_mean")])
                    self.assertEqual(row["data"], cell["value"])
                    self.assertEqual(row["mc"], cell["mc_a_parallel"])
                    self.assertEqual(row["supported"], cell["supported_prediction"])
            self.assertEqual(count, 45)

    def test_negative_tiny_and_zero_compatible_data_are_retained(self):
        cell = copy.deepcopy(self.base["targets"]["D"]["cells"][5])
        cell.update(value=-1e-14, mc_a_parallel=.2, mc_stat=.03,
                    stat=.01, total_error=.02, mc_supported=True)
        row = plots.adapt_rows([cell])[0]
        self.assertEqual(row["ratio"], -2e13)
        self.assertEqual(row["ratio_mc_stat"], 3e12)
        self.assertTrue(row["denominator_zero_compatible"])
        cell["value"] = 0.
        self.assertIsNone(plots.adapt_rows([cell])[0]["ratio"])

    def test_missing_MC_is_distinct_from_unsupported_low_Q2(self):
        cell = copy.deepcopy(self.base["targets"]["P"]["cells"][5])
        cell.update(mc_a_parallel=None, mc_stat=None, mc_supported=False)
        missing = plots.adapt_rows([cell])[0]
        self.assertTrue(missing["supported"])
        self.assertEqual(missing["ratio_missing_reason"], "Missing MC prediction")
        low = plots.adapt_rows([self.base["targets"]["P"]["cells"][0]])[0]
        self.assertFalse(low["supported"])
        self.assertIsNone(low["ratio"])

    def test_full_range_limits_contain_all_supported_MC_and_data_errors(self):
        for block in self.snapshot["targets"].values():
            for panel in block["panels"]:
                low, high = panel["full_range_ratio_limits"]
                for row in panel["ratio_rows"]:
                    if row["ratio"] is None:
                        continue
                    for value in (row["ratio"]-row["ratio_mc_stat"],
                                  row["ratio"]+row["ratio_mc_stat"],
                                  1-row["data_relative_total"], 1+row["data_relative_total"]):
                        self.assertLess(low, value)
                        self.assertGreater(high, value)
        first = self.snapshot["targets"]["P"]["panels"][8]["full_range_ratio_limits"]
        last = self.snapshot["targets"]["D"]["panels"][17]["full_range_ratio_limits"]
        self.assertNotEqual(first, last)  # The full-range companion has independent scales.

    def test_all_new_outputs_are_distinct_from_existing_absolute_outputs(self):
        self.assertEqual(len(plots.REQUIRED_OUTPUTS), 84)
        self.assertFalse(plots.REQUIRED_OUTPUTS & born.REQUIRED_OUTPUTS)
        self.assertIn("AParallel_P_BornCells_Ratio_5x4.pdf", plots.REQUIRED_OUTPUTS)
        self.assertIn("AParallel_D_BornCells_Ratio_5x4_FullRange.png", plots.REQUIRED_OUTPUTS)
        for target in ("P", "D"):
            for index in range(1, 20):
                self.assertIn(f"AParallel_{target}_X{index:02d}_WithRatio.pdf", plots.REQUIRED_OUTPUTS)

    def test_render_dispatches_all_38_single_and_four_atlas_figures(self):
        try:
            import matplotlib
        except ImportError:
            self.skipTest("Matplotlib unavailable")
        def single(snapshot, target, panel, output):
            return [output/f"AParallel_{target}_X{panel['panel_index']+1:02d}_WithRatio.{suffix}"
                    for suffix in ("pdf", "png")]
        def atlas(snapshot, target, output, full_range=False):
            stem = plots.FIGURE_STEMS[target]+("_FullRange" if full_range else "")
            return [output/f"{stem}.{suffix}" for suffix in ("pdf", "png")]
        with tempfile.TemporaryDirectory() as temporary, \
             mock.patch.object(plots, "_render_single", side_effect=single) as singles, \
             mock.patch.object(plots, "_render_atlas", side_effect=atlas) as atlases:
            paths = plots.render(self.base, Path(temporary))
        self.assertEqual(singles.call_count, 38)
        self.assertEqual(atlases.call_count, 4)
        self.assertEqual({path.name for path in paths}, plots.REQUIRED_OUTPUTS)


if __name__ == "__main__":
    unittest.main()
