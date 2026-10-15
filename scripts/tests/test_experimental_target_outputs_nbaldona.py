"""Independent target combinations from shared proton/neutron samples."""
from __future__ import annotations

import copy
import contextlib
import io
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest.mock import Mock, patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import run_experimental_campaign as campaign


class TargetOutputTests(unittest.TestCase):
    def setUp(self):
        self.measurement = copy.deepcopy(campaign.get_measurement("COMPASS_2017_I1501480"))
        raw = self.measurement["raw_observables"]["Q2GT1"]
        outputs = self.measurement["outputs"]["Q2GT1"]
        self.measurement["raw_observables"] = {key: dict(raw) for key in ("He3", "Neutron")}
        self.measurement["outputs"] = {key: dict(outputs) for key in ("He3", "Neutron")}
        for selection, output in self.measurement["outputs"].items():
            for key in ("a1", "apar", "parity_pp_mm", "parity_pm_mp"):
                output[key] = f"/TEST/{key}_{selection}"
            output.update(
                observable="A1He3" if selection == "He3" else "A1n",
                target="helium-3" if selection == "He3" else "neutron",
                target_model="effective helium-3" if selection == "He3" else "free neutron",
            )
        self.measurement["diagnostics"] = ["Accepted_X"]
        self.helium = {
            "unpolarized": {"P": 2/3, "N": 1/3},
            "longitudinal": {"P": 2*(-0.027)/3, "N": 0.87/3},
        }
        self.measurement["combination"]["target_components"] = self.helium
        self.measurement["outputs"]["Neutron"]["target_components"] = {
            "unpolarized": {"P": 0, "N": 1}, "longitudinal": {"P": 0, "N": 1},
        }
        self.family = campaign.campaign_family_specs(self.measurement, False)["nominal"]
        self.proton_variance = 0.01

    def load(self, groups, directory, analysis, name, family, component, helicities, orders, *args):
        # sigma_h = U + sign(h)*L; a fixed D=1/2 makes the expected
        # A1 exactly 2L/U. The covariance proxy's variance is var(w)/D.
        u, ll = (10.0, 2.0) if component == "P" else (20.0, -3.0)
        variance = self.proton_variance if component == "P" else 0.02
        answer = {}
        for helicity in helicities:
            value = u if helicity == "00" else u + (ll if helicity in ("PP", "MM") else -ll)
            factor, variance_factor = 1.0, 1.0
            if name.startswith("SigmaOverD"):
                factor, variance_factor = 2.0, 4.0
            elif name.startswith("CovarianceProxy"):
                factor, variance_factor = 2**0.5, 2.0
            answer[helicity] = campaign.BinSeries([0.03, 0.04], [value*factor], [variance*variance_factor])
        return answer

    def products(self):
        with patch.object(campaign, "_load_component_series", side_effect=self.load):
            return campaign._load_asymmetry_family_products({}, Path("."), self.measurement, "nominal", self.family)

    def test_helium_and_free_neutron_share_samples_but_not_polarization_factors(self):
        results, diagnostics = self.products()
        # Compute the nuclear ratio from whole-nucleus cross sections, avoiding
        # the per-nucleon coefficients used in the implementation.
        expected_he = 2*(2*(-0.027)*2 + 0.87*(-3))/(2*10 + 20)
        self.assertAlmostEqual(results["He3"]["a1_values"][0], expected_he)
        self.assertAlmostEqual(results["Neutron"]["a1_values"][0], -0.3)
        self.assertAlmostEqual(results["Neutron"]["apar_values"][0], -0.15)
        self.assertAlmostEqual(diagnostics["Accepted_X"].values[0], (2*10+20)/3)

    def test_proton_noise_does_not_inflate_free_neutron_error(self):
        before, _ = self.products()
        self.proton_variance = 1e6
        after, _ = self.products()
        self.assertEqual(before["Neutron"]["a1_errors"], after["Neutron"]["a1_errors"])
        self.assertGreater(after["He3"]["a1_errors"][0], before["He3"]["a1_errors"][0])

    def test_diagnostics_do_not_depend_on_output_order(self):
        _, before = self.products()
        self.measurement["raw_observables"] = dict(reversed(list(self.measurement["raw_observables"].items())))
        _, after = self.products()
        self.assertEqual(before, after)

    def test_direct_unpolarized_comparison_uses_same_target_as_each_output(self):
        with patch.object(campaign, "_load_component_series", side_effect=self.load):
            results, diagnostics = campaign._load_direct_unpolarized_products(
                {}, Path("."), self.measurement, "unpolarized_nlo", self.family
            )
        self.assertAlmostEqual(results["He3"].values[0], (2*10+20)/3)
        self.assertEqual(results["Neutron"].values, [20.0])
        self.assertAlmostEqual(diagnostics["Accepted_X"].values[0], (2*10+20)/3)

    def test_original_deuteron_combination_is_retained_without_overrides(self):
        original = campaign.get_measurement("COMPASS_2017_I1501480")
        for selection in original["outputs"]:
            self.assertEqual(campaign._target_component_coefficients(original, self.family, "longitudinal", selection), {"P": 0.4625, "N": 0.4625})

    @staticmethod
    def estimate(yoda, edges, path, values, errors, annotations):
        # Isolate output metadata checks from the optional native YODA binding.
        return Mock(path=path, annotations=dict(annotations))

    def test_output_annotations_distinguish_helium_from_free_neutron(self):
        results, diagnostics = self.products()
        with patch.object(campaign, "_estimate_from_values", side_effect=self.estimate):
            objects = campaign._asymmetry_family_objects(
                None, self.measurement, {}, results, diagnostics,
                campaign._family_annotations(self.measurement, "nominal", self.family),
            )
        by_path = {obj.path: obj.annotations for obj in objects}
        helium, neutron = by_path["/TEST/a1_He3"], by_path["/TEST/a1_Neutron"]
        self.assertEqual(helium["Observable"], "A1He3")
        self.assertEqual(helium["TargetComponents"], "P,N")
        self.assertEqual(neutron["Observable"], "A1n")
        self.assertEqual(neutron["TargetComponents"], "N")
        self.assertEqual(neutron["TargetModel"], "free neutron")

    def test_e142_descriptor_produces_both_targets_in_every_window(self):
        self.measurement = campaign.get_measurement("E142_1996_I424108")
        self.family = campaign.campaign_family_specs(self.measurement, False)["nominal"]
        self.assertNotIn("selection_target_components", self.measurement["combination"])
        results, diagnostics = self.products()
        with patch.object(campaign, "_estimate_from_values", side_effect=self.estimate):
            objects = campaign._asymmetry_family_objects(
                None, self.measurement, {}, results, diagnostics,
                campaign._family_annotations(self.measurement, "nominal", self.family),
            )
        by_path = {obj.path: obj.annotations for obj in objects}
        self.assertEqual(len(results), 6)
        for selection, output in self.measurement["outputs"].items():
            neutron = selection.startswith("Neutron")
            self.assertAlmostEqual(results[selection]["a1_values"][0], -.3 if neutron else -.1359)
            annotations = by_path[output["a1"]]
            self.assertEqual(annotations["Observable"], "A1n" if neutron else "A1He3")
            self.assertEqual(annotations["TargetComponents"], "N" if neutron else "P,N")
            self.assertEqual(annotations["TargetModel"], output["target_model"])
            self.assertEqual(annotations["Acceptance"], output["acceptance"])

    def test_scale_variations_preserve_target_metadata_and_change_only_paths(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            descriptor = copy.deepcopy(self.measurement)
            descriptor["reference"]["snapshot"] = "reference.json"
            (root / "reference.json").write_text("{}")
            campaign.atomic_write_json(root / campaign.MANIFEST_NAME, {
                "configuration": {"variation_points": [[0, 0, 1.0], [0, 0, 2.0]]},
            })
            with patch.object(campaign, "DISPOL_ROOT", root), \
                    patch.object(campaign, "campaign_directory", return_value=root), \
                    patch.object(campaign, "_assert_manifest_signatures_current"), \
                    patch.object(campaign, "_require_complete_matrix", return_value={}), \
                    patch.object(campaign, "campaign_family_specs", return_value={"nominal": self.family}), \
                    patch.object(campaign, "_load_component_series", side_effect=self.load), \
                    patch.object(campaign, "_import_yoda", return_value=None), \
                    patch.object(campaign, "_estimate_from_values", side_effect=self.estimate), \
                    patch.object(campaign, "_asymmetry_family_objects", wraps=campaign._asymmetry_family_objects) as output_writer, \
                    contextlib.redirect_stdout(io.StringIO()):
                campaign.postprocess_campaign(types.SimpleNamespace(tag="test", dry_run=True), descriptor)
            variation = output_writer.call_args_list[1].args[1]
            for selection, original in descriptor["outputs"].items():
                varied = variation["outputs"][selection]
                metadata = set(original) - {"a1", "apar", "parity_pp_mm", "parity_pm_mp"}
                self.assertEqual({key: varied[key] for key in metadata},
                                 {key: original[key] for key in metadata})
                self.assertIn("/VARIATIONS/p000-u000-mu2/", varied["a1"])
                self.assertEqual(original["a1"], f"/TEST/a1_{selection}")

    def test_legacy_neutron_annotations_survive_the_merged_writer(self):
        for identifier in ("E154_1997_I443170", "HERMES_1997_I440904"):
            with self.subTest(measurement=identifier):
                self.measurement = campaign.get_measurement(identifier)
                self.family = campaign.campaign_family_specs(self.measurement, False)["nominal"]
                results, diagnostics = self.products()
                with patch.object(campaign, "_estimate_from_values", side_effect=self.estimate):
                    objects = campaign._asymmetry_family_objects(
                        None, self.measurement, {"observable": "A1"}, results, diagnostics,
                        campaign._family_annotations(self.measurement, "nominal", self.family),
                    )
                by_path = {obj.path: obj.annotations for obj in objects}
                for selection, output in self.measurement["outputs"].items():
                    actual = by_path[output["a1"]]
                    for key, value in self.measurement["output_annotations"][selection].items():
                        self.assertEqual(actual[key], value)
                    self.assertEqual(actual["TargetComponents"], "N")
                    self.assertEqual(by_path[output["apar"]]["Observable"], "A_parallel")

    def test_structured_output_metadata_overrides_legacy_annotations(self):
        self.measurement["output_annotations"] = {
            "He3": {"Observable": "legacy", "TargetModel": "legacy", "CustomNote": "retained"}
        }
        results, diagnostics = self.products()
        with patch.object(campaign, "_estimate_from_values", side_effect=self.estimate):
            objects = campaign._asymmetry_family_objects(
                None, self.measurement, {"observable": "A1"}, results, diagnostics,
                campaign._family_annotations(self.measurement, "nominal", self.family),
            )
        actual = {obj.path: obj.annotations for obj in objects}["/TEST/a1_He3"]
        self.assertEqual(actual["Observable"], "A1He3")
        self.assertEqual(actual["TargetModel"], "effective helium-3")
        self.assertEqual(actual["CustomNote"], "retained")

    def test_validation_rejects_missing_components_nonfinite_and_extra_channels(self):
        source = Path(self.measurement["id"] + ".json")
        campaign._validate_measurement(self.measurement, source)
        for bad in (
            {"unpolarized": {"N": 1}, "longitudinal": {"N": 1}},
            {"unpolarized": {"P": 0, "N": float("nan")}, "longitudinal": {"P": 0, "N": 1}},
            {**self.helium, "extra": {}},
        ):
            with self.subTest(bad=bad), self.assertRaises(campaign.CampaignError):
                measurement = copy.deepcopy(self.measurement)
                measurement["outputs"]["Neutron"]["target_components"] = bad
                campaign._validate_measurement(measurement, source)

    def test_validation_rejects_output_without_matching_raw_selection(self):
        self.measurement["outputs"]["typo"] = self.measurement["outputs"].pop("Neutron")
        with self.assertRaises(campaign.CampaignError):
            campaign._validate_measurement(self.measurement, Path(self.measurement["id"] + ".json"))


if __name__ == "__main__":
    unittest.main()
