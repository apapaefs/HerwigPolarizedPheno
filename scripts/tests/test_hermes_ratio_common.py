"""Signed ratio displays keep separate uncertainties and disclose clipping."""
from __future__ import annotations

import math
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import hermes_ratio_common as ratios


def display_row(data, stat, total, mc, mcstat, supported=True, **geometry):
    return {"low": .1, "high": .2, "marker": .15, "data": data,
            "supported": supported, **geometry,
            **ratios.ratio_record(data, stat, total, mc, mcstat, supported)}


class RatioRecordTests(unittest.TestCase):
    def test_signed_denominator_and_fixed_data_mc_error(self):
        # Holding the measured central value fixed, shift MC up/down by one
        # MC error. The resulting half-distance is positive for either sign.
        for data, mc in ((2., 4.), (-2., 4.), (2., -4.), (-2., -4.)):
            with self.subTest(data=data, mc=mc):
                low, high = sorted(((mc - .3) / data, (mc + .3) / data))
                result = ratios.ratio_record(data, .2, .5, mc, .3)
                self.assertAlmostEqual(result["ratio"], (high + low) / 2)
                self.assertAlmostEqual(result["ratio_mc_stat"], (high - low) / 2)
                self.assertEqual(result["data_relative_stat"], .1)
                self.assertEqual(result["data_relative_total"], .25)
                self.assertIsNone(result["ratio_missing_reason"])
        self.assertEqual(ratios.ratio_record(-2., .2, .5, 4., .3)["ratio"], -2.)

    def test_zero_mc_retains_nonzero_mc_uncertainty(self):
        result = ratios.ratio_record(-.25, .02, .03, 0., .1)
        self.assertEqual(result["ratio"], 0.)
        self.assertEqual(result["ratio_mc_stat"], .4)
        self.assertFalse(result["denominator_zero_compatible"])

    def test_exact_zero_data_is_explicitly_undefined(self):
        result = ratios.ratio_record(0., .01, .02, .1, .003)
        self.assertIsNone(result["ratio"])
        self.assertIsNone(result["ratio_mc_stat"])
        self.assertIsNone(result["data_relative_stat"])
        self.assertIsNone(result["data_relative_total"])
        self.assertEqual(result["ratio_missing_reason"], "Data central value is zero")
        self.assertTrue(result["denominator_zero_compatible"])

    def test_zero_compatible_data_still_has_signed_finite_ratio(self):
        for data, expected in ((.001, 11.), (-.001, -11.)):
            result = ratios.ratio_record(data, .009, .01, .011, .001)
            self.assertAlmostEqual(result["ratio"], expected)
            self.assertEqual(result["ratio_mc_stat"], 1.)
            self.assertEqual(result["data_relative_total"], 10.)
            self.assertTrue(result["denominator_zero_compatible"])
            self.assertIsNone(result["ratio_missing_reason"])
        self.assertTrue(ratios.ratio_record(.01, .005, .01, .01, .001)
                        ["denominator_zero_compatible"])

    def test_data_uncertainty_is_separate_from_mc_error(self):
        narrow = ratios.ratio_record(.2, .01, .02, .1, .004)
        wide = ratios.ratio_record(.2, .05, .3, .1, .004)
        self.assertEqual(narrow["ratio"], wide["ratio"])
        self.assertEqual(narrow["ratio_mc_stat"], wide["ratio_mc_stat"])
        self.assertEqual(wide["ratio_mc_stat"], .02)
        self.assertGreater(wide["data_relative_total"], narrow["data_relative_total"])

    def test_support_and_missing_prediction_are_distinct(self):
        for supported, reason in ((False, "No MC support"), (True, "Missing MC prediction")):
            record = ratios.ratio_record(.2, .01, .02, None, None, supported)
            self.assertIsNone(record["ratio"])
            self.assertIsNone(record["ratio_mc_stat"])
            self.assertEqual(record["ratio_missing_reason"], reason)
            self.assertAlmostEqual(record["data_relative_total"], .1)
        for mc, error in ((None, .01), (.2, None)):
            with self.assertRaisesRegex(ValueError, "present or missing together"):
                ratios.ratio_record(.2, .01, .02, mc, error)
        with self.assertRaisesRegex(ValueError, "Unsupported"):
            ratios.ratio_record(.2, .01, .02, .3, .01, False)
        with self.assertRaisesRegex(ValueError, "boolean"):
            ratios.ratio_record(.2, .01, .02, .3, .01, 1)

    def test_invalid_values_fail_explicitly(self):
        valid = [.2, .01, .02, .3, .04]
        invalid = ((0, None), (0, True), (0, math.nan), (0, math.inf),
                   (0, "not a number"), (1, -.001), (1, math.nan),
                   (2, .005), (2, -1e-13), (2, math.inf),
                   (3, math.nan), (3, False), (4, -.01), (4, math.inf))
        for index, value in invalid:
            arguments = list(valid)
            arguments[index] = value
            with self.subTest(index=index, value=value), self.assertRaises(ValueError):
                ratios.ratio_record(*arguments)
        # Nonnegative total remains necessary when stat=0, independently of
        # any tolerance permitted when comparing two nearly equal errors.
        with self.assertRaises(ValueError):
            ratios.ratio_record(.2, 0., -1e-13, .3, .04)


class RatioLimitsTests(unittest.TestCase):
    def test_focus_is_explicit_and_full_range_contains_every_finite_extent(self):
        rows = [display_row(.1, .02, .03, .13, .01),
                display_row(.0027, .0584, math.hypot(.0584, .0053), .0646, .0048),
                display_row(-.009, .01, math.hypot(.01, .0004), .0052, .001),
                display_row(.2, .01, .02, None, None)]
        self.assertEqual(ratios.ratio_limits(rows), [-1., 3.])
        low, high = ratios.ratio_limits(rows, focus=False)
        for row in rows:
            if row["ratio"] is None:
                continue
            for endpoint in (row["ratio"] - row["ratio_mc_stat"],
                             row["ratio"] + row["ratio_mc_stat"],
                             1 - row["data_relative_stat"], 1 + row["data_relative_stat"],
                             1 - row["data_relative_total"], 1 + row["data_relative_total"]):
                self.assertLess(low, endpoint)
                self.assertGreater(high, endpoint)
        self.assertLess(low, 0.)
        self.assertGreater(high, 2.)

    def test_empty_full_range_still_contains_unity_and_zero(self):
        low, high = ratios.ratio_limits([], focus=False)
        self.assertLess(low, 0.)
        self.assertGreater(high, 1.)


class RatioRenderingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError as exc:
            raise unittest.SkipTest("Matplotlib is unavailable") from exc
        cls.plt = plt

    def axis(self):
        figure, axis = self.plt.subplots()
        self.addCleanup(self.plt.close, figure)
        return axis

    def test_in_range_zero_compatible_ratio_has_hollow_marker(self):
        axis = self.axis()
        row = display_row(-.009, .01, .011, .0052, .001)
        ratios.draw_ratio(axis, [row])
        markers = [line for line in axis.lines if line.get_marker() == "o"]
        self.assertEqual(len(markers), 1)
        self.assertEqual(markers[0].get_markerfacecolor(), "white")
        self.assertAlmostEqual(markers[0].get_ydata()[0], row["ratio"])

    def test_every_mc_and_data_truncation_is_marked(self):
        axis = self.axis()
        # Center above the focus and MC error below its lower edge: both
        # directions must be marked, in addition to clipped data bands.
        row = display_row(1., 2., 5., 5., 10.)
        ratios.draw_ratio(axis, [row])
        orange = {line.get_marker() for line in axis.lines if line.get_color() == ratios.COLOR}
        gray = {line.get_marker() for line in axis.lines if line.get_color() == "0.48"}
        self.assertTrue({"^", "v"}.issubset(orange))
        self.assertEqual(gray, {"^", "v"})
        self.assertIn("5", [text.get_text() for text in axis.texts])
        for patch in axis.patches:
            self.assertGreaterEqual(patch.get_y(), ratios.FOCUS[0])
            self.assertLessEqual(patch.get_y() + patch.get_height(), ratios.FOCUS[1])

    def test_full_range_needs_no_clipping_arrows(self):
        axis = self.axis()
        rows = [display_row(.001, .009, .01, .011, .001),
                display_row(-.009, .01, .011, .0052, .001, low=.2, high=.3, marker=.25)]
        limits = ratios.ratio_limits(rows, focus=False)
        ratios.draw_ratio(axis, rows, limits)
        self.assertFalse(any(line.get_marker() in ("^", "v") for line in axis.lines))
        for patch in axis.patches:
            self.assertGreater(patch.get_y(), limits[0])
            self.assertLess(patch.get_y() + patch.get_height(), limits[1])

    def test_missing_and_exact_zero_ratios_are_labelled(self):
        axis = self.axis()
        rows = [display_row(0., .01, .02, .1, .003),
                display_row(.1, .01, .02, None, None, low=.2, high=.3, marker=.25)]
        ratios.draw_ratio(axis, rows)
        self.assertEqual({text.get_text() for text in axis.texts}, {"data = 0", "MC missing"})
        self.assertFalse(any(line.get_marker() == "o" for line in axis.lines))

    def test_invalid_display_geometry_and_limits_fail(self):
        for changes in ({"low": 0.}, {"high": .05}, {"marker": .21}):
            with self.subTest(changes=changes), self.assertRaisesRegex(ValueError, "geometry"):
                ratios.draw_ratio(self.axis(), [display_row(.1, .01, .02, .1, .01, **changes)])
        with self.assertRaisesRegex(ValueError, "axis limits"):
            ratios.draw_ratio(self.axis(), [], (1., 1.))


if __name__ == "__main__":
    unittest.main()
