"""Tests of the conditional fitted-UU projection, independent of fit numerics."""
from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import hermes_apar_integrated as projection


ROOT = Path(__file__).resolve().parents[2]


def fixtures():
    reference = json.loads((ROOT / "data/experimental/HERMES_2007_I726689/born-apar-reference.json").read_text())
    rows = []
    for dataset in reference["datasets"]:
        for point in dataset["points"]:
            rows.append({
                "selection": dataset["selection"], "bin": point["bin"],
                "bin_low": point["q2_low"], "bin_high": point["q2_high"],
                "supported": point["supported_prediction"],
                "a_parallel": 0.1 if point["supported_prediction"] else None,
                "a_parallel_stat": 0.02 if point["supported_prediction"] else None,
                # Differential MC density columns must not supply fitted weights.
                "sigma_uu_pb": 9999999., "sigma_ll_pb": 9999999.,
            })
    return reference, {"measurement": projection.MEASUREMENT, "tag": "toy", "bins": rows}


def area(target, xlow, xhigh, qlow, qhigh):
    return (xhigh - xlow) * (qhigh - qlow) * (1.0 if target == "P" else 2.0)


class HermesAparProjectionTests(unittest.TestCase):
    def setUp(self):
        self.reference, self.summary = fixtures()
        self.result = projection.build_projection_snapshot(self.reference, self.summary, area)

    def test_whole_cell_x_projection_uses_integrals_and_physical_edges(self):
        target = self.result["targets"]["P"]
        row = target["bins"][4]  # X09 has widths .505, .760, 17.735.
        weights = target["weights"][4]
        self.assertEqual((row["bin_low"], row["bin_high"]), (.05681, .07265))
        self.assertFalse(row["cell_constant_assumption"])
        self.assertEqual(row["source_global_bins"], [13, 14, 15])
        for source, width in ((12, .505), (13, .760), (14, 17.735)):
            self.assertAlmostEqual(weights[source], width / 19.)
        expected = sum(weights[p["global_bin"] - 1] * p["value"]
                       for dataset in self.reference["datasets"] if dataset["target"] == "P"
                       for p in dataset["points"])
        self.assertAlmostEqual(row["a_parallel"], expected)
        self.assertAlmostEqual(row["mc_a_parallel"], .1)

    def test_full_source_covariance_is_retained_and_split_bins_correlate(self):
        # Identity covariance isolates correlations induced by reusing a cell.
        reference = copy.deepcopy(self.reference)
        for target in ("P", "D"):
            reference["statistical_covariance"][target]["matrix"] = [
                [1. if i == j else 0. for j in range(45)] for i in range(45)]
        result = projection.build_projection_snapshot(reference, self.summary, area)
        block = result["targets"]["P"]
        left, right = 15, 16  # Q2 1–1.5 and 1.5–2 split many of the same cells.
        expected = sum(a * b for a, b in zip(block["weights"][left], block["weights"][right]))
        self.assertGreater(expected, 0.)
        self.assertAlmostEqual(block["statistical_covariance"][left][right], expected)
        self.assertTrue(block["bins"][left]["cell_constant_assumption"])
        # A separate off-diagonal input changes the projection variance.
        reference["statistical_covariance"]["P"]["matrix"][12][13] = .5
        reference["statistical_covariance"]["P"]["matrix"][13][12] = .5
        changed = projection.build_projection_snapshot(reference, self.summary, area)
        w = block["weights"][4]
        self.assertAlmostEqual(changed["targets"]["P"]["bins"][4]["stat"] ** 2,
                               block["bins"][4]["stat"] ** 2 + w[12] * w[13])

    def test_constant_asymmetry_is_preserved_and_unavailable_cells_excluded(self):
        reference = copy.deepcopy(self.reference)
        for dataset in reference["datasets"]:
            for cell in dataset["points"]:
                cell["value"] = .23 if cell["supported_prediction"] else -1000.
        result = projection.build_projection_snapshot(reference, self.summary, area)
        for target in ("P", "D"):
            block = result["targets"][target]
            unsupported = [cell["global_bin"] - 1 for cell in block["source_cells"]
                           if not cell["supported_prediction"]]
            self.assertEqual(len(unsupported), 8)
            for row, weights in zip(block["bins"], block["weights"]):
                self.assertAlmostEqual(sum(weights), 1.)
                self.assertAlmostEqual(row["a_parallel"], .23)
                self.assertTrue(all(weights[index] == 0 for index in unsupported))
                self.assertTrue(all(weight >= 0 for weight in weights))
                self.assertLessEqual(row["bin_low"], row["marker"])
                self.assertGreaterEqual(row["bin_high"], row["marker"])

    def test_systematic_bound_includes_published_normalization_once(self):
        for target in ("P", "D"):
            block = self.result["targets"][target]
            for row, weights in zip(block["bins"], block["weights"]):
                expected = sum(weight * cell["systematic_combined"]
                               for weight, cell in zip(weights, block["source_cells"]))
                self.assertAlmostEqual(row["systematic_upper_bound"], expected)
                self.assertGreaterEqual(row["systematic_upper_bound"], row["systematic_diagonal_diagnostic"])
            self.assertIn("neither added nor subtracted", block["normalization_policy"])

    def test_zero_fitted_yield_remains_masked_in_covariance(self):
        result = projection.build_projection_snapshot(self.reference, self.summary,
                                                       lambda *args: 0.)
        for target in ("P", "D"):
            block = result["targets"][target]
            self.assertTrue(all(not row["supported"] for row in block["bins"]))
            self.assertTrue(all(row["a_parallel"] is None for row in block["bins"]))
            self.assertTrue(all(row["stat"] is None for row in block["bins"]))
            self.assertTrue(all(value is None for row in block["statistical_covariance"] for value in row))
            self.assertTrue(all(value == 0 for row in block["weights"] for value in row))

    def test_theory_and_data_share_weights_and_controls_preserve_operators(self):
        summary = copy.deepcopy(self.summary)
        reference_values = {(d["selection"], p["bin"]): p["value"]
                            for d in self.reference["datasets"] for p in d["points"]}
        for row in summary["bins"]:
            if row["supported"]:
                row["a_parallel"] = reference_values[row["selection"], row["bin"]]
        result = projection.build_projection_snapshot(
            self.reference, summary, area,
            controls={"shape": lambda target, xl, xh, ql, qh: area(target, xl, xh, ql, qh) * (ql + qh)},
        )
        for target in ("P", "D"):
            for row in result["targets"][target]["bins"]:
                self.assertAlmostEqual(row["a_parallel"], row["mc_a_parallel"])
                self.assertAlmostEqual(row["weight_controls"]["shape"]["data_minus_mc_shift"], 0.)
            self.assertEqual(result["controls"]["shape"]["targets"][target]["projection_order"],
                             result["targets"][target]["projection_order"])

    def test_invalid_geometry_support_duplicates_and_integrals_rejected(self):
        bad = copy.deepcopy(self.summary)
        bad["bins"].append(copy.deepcopy(bad["bins"][0]))
        with self.assertRaisesRegex(projection.ProjectionError, "Repeated theory"):
            projection.build_projection_snapshot(self.reference, bad, area)
        bad = copy.deepcopy(self.summary)
        bad["bins"][5]["bin_low"] = .999
        with self.assertRaisesRegex(projection.ProjectionError, "boundaries differ"):
            projection.build_projection_snapshot(self.reference, bad, area)
        bad = copy.deepcopy(self.summary)
        bad["bins"][0]["supported"] = True
        with self.assertRaisesRegex(projection.ProjectionError, "support mismatch"):
            projection.build_projection_snapshot(self.reference, bad, area)
        for invalid in (-1., float("nan")):
            with self.assertRaises(projection.ProjectionError):
                projection.build_projection_snapshot(self.reference, self.summary, lambda *args: invalid)

    def test_numerical_outputs_do_not_overwrite_existing_reconstruction(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "derived"
            projection.write_outputs(self.result, output)
            saved = json.loads((output / "integrated.json").read_text())
            self.assertEqual(saved["targets"]["P"]["weights"], self.result["targets"]["P"]["weights"])
            self.assertEqual(len((output / "integrated.csv").read_text().splitlines()), 85)
            with self.assertRaisesRegex(projection.ProjectionError, "not empty"):
                projection.write_outputs(self.result, output)


if __name__ == "__main__":
    unittest.main()
