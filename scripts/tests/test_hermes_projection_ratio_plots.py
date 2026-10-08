"""Keep published A1 and model-assisted Born ratios paired to their own estimators."""
from __future__ import annotations

import copy
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import hermes_projection_ratio_plots as plots


def inputs():
    summary = {"measurement": plots.MEASUREMENT, "tag": "fixture", "bins": [],
               "acceptance": "rectangle_intersect_polar_ring"}
    reference = {"measurement": plots.MEASUREMENT, "datasets": []}
    integrated = {"measurement": plots.MEASUREMENT, "campaign_tag": "fixture",
                  "fit": {"cuts": {"acceptance": summary["acceptance"]}}, "targets": {}}
    for target, selection in (("P", "Q2GT1"), ("D", "D_Q2GT1")):
        points = [dict(bin=i, x_mean=x, value=value, stat=.03, systematic_combined=.04,
                       systematics={"experimental": .024, "parameterization": .032, "evolution": 0.})
                  for i, x, value in ((1, .15, -.02), (2, .3, .5))]
        reference["datasets"].append(dict(target=target, bin_edges=[.1, .2, .4], points=points))
        for i, low, high in ((1, .1, .2), (2, .2, .4)):
            summary["bins"].append(dict(selection=selection, bin=i, bin_low=low,
                                       bin_high=high, a1=.3, a1_stat=.02,
                                       a_parallel=99., supported=True))
            summary["bins"].append(dict(selection=selection.replace("1", "4"), bin=i,
                                       bin_low=low, bin_high=high, a1=123., a1_stat=456.))
        rows = []
        for projection in ("x", "q2"):
            for cut in ("Q2GT1", "Q2GT4"):
                rows.append(dict(id=f"{target}_{projection}_{cut}_1", bin=1,
                                 projection=projection, selection=cut, bin_low=.1, bin_high=.4,
                                 marker=.2, a_parallel=.2, stat=.06,
                                 systematic_upper_bound=.08, stat_plus_systematic_upper_bound=.1,
                                 mc_a_parallel=.12, mc_stat_independent_cell_approx=.01,
                                 supported=True, mc_supported=True))
        integrated["targets"][target] = {"bins": rows}
    return summary, reference, integrated


class ProjectionRatioTests(unittest.TestCase):
    def test_ten_comparisons_use_correct_observables_weights_and_published_means(self):
        summary, reference, integrated = inputs()
        before = copy.deepcopy((summary, reference, integrated))
        comparisons = plots.build_comparisons(summary, reference, integrated)
        self.assertEqual((summary, reference, integrated), before)
        self.assertEqual(len(comparisons), 10)
        self.assertEqual(tuple(item["id"] for item in comparisons), plots.COMPARISON_IDS)
        self.assertEqual(len(plots.REQUIRED_OUTPUTS), 40)
        a1 = [item for item in comparisons if item["kind"] == "a1"]
        self.assertEqual(len(a1), 2)
        self.assertTrue(all("Q2GT1" in item["id"] for item in a1))
        first = a1[0]["rows"][0]
        self.assertEqual((first["marker"], first["low"], first["high"]), (.15, .1, .2))
        self.assertAlmostEqual(first["total"], .05)
        self.assertEqual(first["mc"], .3)
        self.assertAlmostEqual(first["ratio"], -15.)
        self.assertTrue(first["denominator_zero_compatible"])
        for item in comparisons[2:]:
            row = item["rows"][0]
            self.assertEqual(row["mc"], .12)  # never the native a_parallel=99
            self.assertAlmostEqual(row["ratio"], .6)
            self.assertEqual(row["total"], .1)
            self.assertEqual(row["systematic_upper_bound"], .08)
            self.assertIn("conservative systematic correlation bound", " ".join(item["notes"]))

    def test_A1_total_matches_original_components_with_explicit_legacy_fallback(self):
        summary, reference, integrated = inputs()
        point = reference["datasets"][0]["points"][0]
        point["systematic_combined"] = .09  # deliberately rounded differently
        comparisons = plots.build_comparisons(summary, reference, integrated)
        self.assertAlmostEqual(comparisons[0]["rows"][0]["total"], .05)
        del point["systematics"]
        comparisons = plots.build_comparisons(summary, reference, integrated)
        self.assertAlmostEqual(comparisons[0]["rows"][0]["total"], (.03**2+.09**2)**.5)
        point["systematics"] = {"experimental": .024}
        with self.assertRaisesRegex(ValueError, "systematic components"):
            plots.build_comparisons(summary, reference, integrated)

    def test_A1_restores_original_Rivet_metadata_and_axis_range(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        comparisons = plots.build_comparisons(*inputs())
        for comparison, target, rivet_id in zip(comparisons[:2], ("p", "d"),
                                                ("d14-x01-y01", "d14-x01-y02")):
            self.assertEqual(comparison["rivet_id"], rivet_id)
            figure, axis = plt.subplots()
            try:
                plots._draw_absolute(axis, comparison)
                self.assertIn("HERMES", axis.get_title())
                self.assertIn("virtual-photon asymmetry", axis.get_title())
                self.assertEqual(axis.get_ylabel(),
                                 rf"$A_1^{{{target}}}$ ($\eta A_2$ neglected in prediction)")
                self.assertEqual(axis.get_ylim(), (-.15, 1.50))
                labels = [text.get_text() for text in axis.get_legend().get_texts()]
                self.assertEqual(labels, ["Data", "Data stat. + syst.",
                                          "NLO+PS (polarized; full spin)", "MC statistical"])
                self.assertFalse(any("PDF" in label or "scale" in label for label in labels))
            finally:
                plt.close(figure)

    def test_zero_data_and_missing_mc_are_explicit_without_masking_small_nonzero_data(self):
        summary, reference, integrated = inputs()
        reference["datasets"][0]["points"][0]["value"] = 0.
        summary["bins"][0]["a1"] = None
        summary["bins"][0]["a1_stat"] = None
        integrated["targets"]["P"]["bins"][0].update(mc_a_parallel=None,
                                                              mc_stat_independent_cell_approx=None,
                                                              mc_supported=False)
        result = plots.build_comparisons(summary, reference, integrated)
        self.assertIsNone(result[0]["rows"][0]["ratio"])
        self.assertIsNone(result[2]["rows"][0]["ratio"])
        self.assertEqual(result[2]["rows"][0]["data"], .2)
        self.assertTrue(result[2]["rows"][0]["supported"])
        self.assertFalse(result[2]["rows"][0]["mc_supported"])

    def test_empty_geometric_projection_keeps_nulls_without_inventing_data(self):
        summary, reference, integrated = inputs()
        row = integrated["targets"]["P"]["bins"][1]
        row.update(a_parallel=None, stat=None, stat_plus_systematic_upper_bound=None,
                   systematic_upper_bound=None, mc_a_parallel=None,
                   mc_stat_independent_cell_approx=None, supported=False, mc_supported=False)
        result = plots.build_comparisons(summary, reference, integrated)
        actual = next(item for item in result if item["id"] == "AParallel_P_vs_x_Q2GT4_reconstructed")["rows"][0]
        self.assertIsNone(actual["data"])
        self.assertIsNone(actual["ratio"])
        self.assertFalse(actual["supported"])
        self.assertEqual(actual["ratio_missing_reason"], "No physical projection support")

    def test_misaligned_campaign_acceptance_and_A1_bins_are_rejected(self):
        for case in ("campaign", "acceptance", "edges", "duplicates"):
            summary, reference, integrated = inputs()
            if case == "campaign": integrated["campaign_tag"] = "other"
            elif case == "acceptance": integrated["fit"]["cuts"]["acceptance"] = "ring"
            elif case == "edges": summary["bins"][0]["bin_high"] = .21
            else: summary["bins"].append(dict(summary["bins"][0]))
            with self.subTest(case=case), self.assertRaises(ValueError):
                plots.build_comparisons(summary, reference, integrated)

    def test_render_writes_both_ranges_with_embedded_pdf_fonts_and_300dpi_png(self):
        from PIL import Image
        comparisons = plots.build_comparisons(*inputs())
        with tempfile.TemporaryDirectory() as temporary:
            paths = plots.render([comparisons[0], comparisons[2]], Path(temporary))
            self.assertEqual(len(paths), 8)
            self.assertTrue(all(path.stat().st_size > 1000 for path in paths))
            self.assertTrue(any("fullrange" in path.name for path in paths))
            for path in paths:
                if path.suffix == ".pdf":
                    self.assertIn(b"/FontFile2", path.read_bytes())
                else:
                    with Image.open(path) as image:
                        self.assertAlmostEqual(image.info["dpi"][0], 300., delta=.1)
                        expected_width = 140. if path.name.startswith("A1_") else 90.
                        self.assertAlmostEqual(image.width / 300 * 25.4, expected_width, delta=.1)
                        self.assertIn("MC/data", image.info["Description"])


if __name__ == "__main__":
    unittest.main()
