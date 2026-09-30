#!/usr/bin/env python3
"""Stored negative-order magnitudes require explicit subtraction, before ratios."""
from __future__ import annotations

import copy
import json
import math
from pathlib import Path
import sys
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import run_experimental_campaign as campaign

MEASUREMENT = "HERMES_2007_I726689"
SIGNS = {"PP": 1, "PM": -1, "MP": -1, "MM": 1}
MAGNITUDES = {
    "P": {"POSNLO": [14., 8., 8., 10.], "NEGNLO": [2., 3., 2., 1.]},
    "N": {"POSNLO": [8., 5., 5., 4.], "NEGNLO": [2., 1., 1., 1.]},
}


def descriptor():
    return copy.deepcopy(campaign.get_measurement(MEASUREMENT))


def magnitude_inputs(family_id="nominal", helicities=tuple(SIGNS), member=0):
    groups, histograms = {}, {}
    for component, orders in MAGNITUDES.items():
        for helicity in helicities:
            for order, values in orders.items():
                value = (sum(values)/4. if helicity == "00" else values[list(SIGNS).index(helicity)])
                # Replica input deliberately has a different normalization.
                value *= member+1
                output = f"{component}-{helicity}-{order}-p{member}.yoda"
                key = (family_id, member, 0, 1., component, helicity, order)
                groups[key] = [{"output_yoda": output, "events": 1}]
                inverse_d = 2. if component == "P" else 3.
                histograms[output] = (value, inverse_d)

    def read(path, histogram_path):
        value, inverse_d = histograms[Path(path).name]
        if "SigmaOverD_" in histogram_path:
            return campaign.BinSeries([.1, .3], [inverse_d*value], [inverse_d**2*value**2])
        if "CovarianceProxy_" in histogram_path:
            return campaign.BinSeries([.1, .3], [0.], [inverse_d*value**2])
        return campaign.BinSeries([.1, .3], [value], [value**2])

    return groups, read


def independent_ratio_variance(selection, a1):
    target_uu = {"P": 1., "N": 0.} if selection == "Q2GT1" else {"P": .5, "N": .5}
    target_ll = {"P": 1., "N": 0.} if selection == "Q2GT1" else {"P": .4625, "N": .4625}
    denominator, numerator = 0., 0.
    for component, orders in MAGNITUDES.items():
        inverse_d = (2. if component == "P" else 3.) if a1 else 1.
        for order, magnitudes in orders.items():
            order_sign = 1. if order == "POSNLO" else -1.
            for helicity, value in zip(SIGNS, magnitudes):
                denominator += target_uu[component]/4. * order_sign * value
                numerator += target_ll[component]/4. * SIGNS[helicity] * inverse_d * order_sign * value
    ratio = numerator/denominator
    variance, covariance = 0., 0.
    for component, orders in MAGNITUDES.items():
        inverse_d = (2. if component == "P" else 3.) if a1 else 1.
        for order, magnitudes in orders.items():
            order_sign = 1. if order == "POSNLO" else -1.
            for helicity, value in zip(SIGNS, magnitudes):
                a = target_uu[component]/4.
                b = target_ll[component]/4.*SIGNS[helicity]
                # Differentiate wrt each independent Poisson event count.
                gradient = (b*inverse_d*order_sign*value-ratio*a*order_sign*value)/denominator
                variance += gradient**2
                covariance += a*b*inverse_d*value**2  # Order sign squared is +1.
    return ratio, variance, covariance


class ExperimentalOrderCoefficientTests(unittest.TestCase):
    def test_optional_map_requires_exact_nominal_keys_and_finite_numbers(self):
        invalid = (None, {}, {"POSNLO": 1.}, {"POSNLO": 1., "NEGNLO": -1., "LO": 1.},
                   {"POSNLO": 1., "NEGNLO": float("nan")},
                   {"POSNLO": 1., "NEGNLO": float("inf")},
                   {"POSNLO": 1., "NEGNLO": True},
                   {"POSNLO": 1., "NEGNLO": "-1"},
                   {"POSNLO": 1., "NEGNLO": 10**400})
        for coefficients in invalid:
            with self.subTest(coefficients=coefficients):
                measurement = descriptor()
                measurement["combination"]["order_coefficients"] = coefficients
                with self.assertRaises(campaign.CampaignError):
                    campaign._validate_measurement(measurement, Path(f"{MEASUREMENT}.json"))
        measurement = descriptor()
        measurement["combination"].pop("order_coefficients")
        self.assertEqual(campaign._order_combination_coefficients(measurement, ("POSNLO", "NEGNLO")),
                         {"POSNLO": 1., "NEGNLO": 1.})

    def test_positive_order_magnitudes_subtract_after_normalizing_shards(self):
        groups = {
            ("nominal", 0, 0, 1., "P", "PP", "POSNLO"): [
                {"output_yoda": "positive1", "events": 100},
                {"output_yoda": "positive2", "events": 300},
            ],
            ("nominal", 0, 0, 1., "P", "PP", "NEGNLO"): [
                {"output_yoda": "negative", "events": 400},
            ],
        }
        series = {
            "positive1": campaign.BinSeries([0., 1.], [10.], [4.]),
            "positive2": campaign.BinSeries([0., 1.], [14.], [9.]),
            "negative": campaign.BinSeries([0., 1.], [2.], [1.]),
        }
        with mock.patch.object(campaign, "read_histogram_series", side_effect=lambda path, name: series[path.name]):
            result = campaign._load_component_series(groups, ROOT, MEASUREMENT, "ordinary",
                                                    component="P", helicities=("PP",),
                                                    order_coefficients={"POSNLO": 1., "NEGNLO": -1.})["PP"]
        self.assertEqual(result.values, [11.])  # .25*10+.75*14 - 2.
        self.assertEqual(result.variances, [.25**2*4. + .75**2*9. + 1.])

    def test_already_signed_input_keeps_legacy_sum_without_sign_inference(self):
        groups = {
            ("nominal", 0, 0, 1., "P", "PP", order): [{"output_yoda": order, "events": 1}]
            for order in ("POSNLO", "NEGNLO")
        }
        def read(path, name):
            return campaign.BinSeries([0., 1.], [10. if path.name == "POSNLO" else -2.], [1.])
        with mock.patch.object(campaign, "read_histogram_series", side_effect=read):
            legacy = campaign._load_component_series(groups, ROOT, MEASUREMENT, "ordinary",
                                                    component="P", helicities=("PP",))["PP"]
            # Negative moment bins can occur even for a magnitude-weight stream.
            # Explicit coefficients act linearly, rather than inspecting/abs'ing values.
            declared = campaign._load_component_series(groups, ROOT, MEASUREMENT, "ordinary",
                                                      component="P", helicities=("PP",),
                                                      order_coefficients={"POSNLO": 1., "NEGNLO": -1.})["PP"]
        self.assertEqual(legacy.values, [8.])
        self.assertEqual(declared.values, [12.])
        self.assertEqual(legacy.variances, declared.variances)

    def test_target_ratios_and_covariance_include_both_order_variances(self):
        measurement = descriptor()
        family = campaign.campaign_family_specs(measurement, False)["nominal"]
        groups, read = magnitude_inputs()
        with mock.patch.object(campaign, "read_histogram_series", side_effect=read):
            results, diagnostics = campaign._load_asymmetry_family_products(groups, ROOT, measurement, "nominal", family)
        self.assertEqual(results["Q2GT1"]["sigma_uu"].values, [8.])
        self.assertEqual(results["Q2GT1"]["sigma_ll"].values, [2.5])
        self.assertEqual(results["D_Q2GT1"]["sigma_uu"].values, [6.125])
        self.assertAlmostEqual(results["D_Q2GT1"]["sigma_ll"].values[0], 1.271875)
        for selection in ("Q2GT1", "D_Q2GT1"):
            for observable in ("apar", "a1"):
                ratio, variance, covariance = independent_ratio_variance(selection, observable == "a1")
                self.assertAlmostEqual(results[selection][f"{observable}_values"][0], ratio)
                self.assertAlmostEqual(results[selection][f"{observable}_errors"][0]**2, variance)
        self.assertEqual(diagnostics["Accepted_X_Q2GT1"].values, [8.])

    def test_direct00_comparison_and_replica_paths_use_same_order_signs(self):
        measurement = descriptor()
        families = campaign.campaign_family_specs(measurement, True)
        direct = families["unpolarized_nlo"]
        groups, read = magnitude_inputs("unpolarized_nlo", ("00",))
        with mock.patch.object(campaign, "read_histogram_series", side_effect=read):
            results, diagnostics = campaign._load_direct_unpolarized_products(
                groups, ROOT, measurement, "unpolarized_nlo", direct
            )
        self.assertEqual(results["Q2GT1"].values, [8.])
        self.assertEqual(results["D_Q2GT1"].values, [6.125])
        self.assertEqual(diagnostics["Accepted_X_Q2GT1"].values, [8.])
        for family_id, member in (("nominal", 1), ("no_real_spin_nlo", 0)):
            groups, read = magnitude_inputs(family_id, member=member)
            with mock.patch.object(campaign, "read_histogram_series", side_effect=read):
                results, _ = campaign._load_asymmetry_family_products(
                    groups, ROOT, measurement, family_id, families[family_id], member
                )
            self.assertEqual(results["Q2GT1"]["sigma_uu"].values, [8.*(member+1)])
            self.assertAlmostEqual(results["Q2GT1"]["a1_values"][0], .625)

    def test_LO_default_and_sign_annotations_are_explicit(self):
        measurement = descriptor()
        families = campaign.campaign_family_specs(measurement, True)
        self.assertEqual(campaign._order_combination_coefficients(measurement, ("LO",)), {"LO": 1.})
        groups = {("polarized_lo", 0, 0, 1., "P", "PP", "LO"): [{"output_yoda": "born", "events": 1}]}
        with mock.patch.object(campaign, "read_histogram_series", return_value=campaign.BinSeries([0., 1.], [3.], [4.])):
            result = campaign._load_component_series(groups, ROOT, MEASUREMENT, "ordinary", "polarized_lo", "P",
                                                    ("PP",), ("LO",), order_coefficients={"POSNLO": 1., "NEGNLO": -1.})["PP"]
        self.assertEqual(result.values, [3.])
        nominal = campaign._family_annotations(measurement, "nominal", families["nominal"])
        born = campaign._family_annotations(measurement, "polarized_lo", families["polarized_lo"])
        self.assertEqual(nominal["OrderCombination"], "POSNLO-NEGNLO")
        self.assertEqual(json.loads(nominal["OrderCoefficients"]), {"POSNLO": 1., "NEGNLO": -1.})
        self.assertIn("magnitudes", nominal["OrderInputConvention"])
        self.assertEqual(born["OrderCombination"], "LO")
        self.assertEqual(json.loads(born["OrderCoefficients"]), {"LO": 1.})

    def test_new_storage_sign_contract_rejects_old_campaign_signature(self):
        measurement = descriptor()
        old = copy.deepcopy(measurement)
        old["combination"].pop("order_coefficients")
        old["physics"].pop("order_input_convention")
        old_signature = campaign.measurement_signature(old)
        self.assertNotEqual(old_signature, campaign.measurement_signature(measurement))
        with self.assertRaises(campaign.CampaignError):
            campaign._assert_manifest_signatures_current(
                {"configuration": {"measurement_signature": old_signature}}, measurement
            )


if __name__ == "__main__":
    unittest.main()
