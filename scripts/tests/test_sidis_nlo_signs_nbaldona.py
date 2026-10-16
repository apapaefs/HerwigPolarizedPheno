"""SIDIS NLO storage conventions, ratio errors, and campaign integration."""
from __future__ import annotations

import argparse
import contextlib
import io
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import run_experimental_campaign as experimental
import run_phenomenology_campaign as campaign
import polarized_sidis_postprocess as sidis


VALUES = {
    "P": {"POSNLO": [14., 8., 8., 10.], "NEGNLO": [2., 3., 2., 1.]},
    "N": {"POSNLO": [8., 5., 5., 4.], "NEGNLO": [2., 1., 1., 1.]},
}
IDS = ("COMPASS_2009_I820721", "COMPASS_2010_I862410", "HERMES_2019_I1698889")


def metadata(path, xsec=10., sumw=10.):
    path.write_text(f"""BEGIN YODA_ESTIMATE0D_V3 /_XSEC
Path: /_XSEC
---
{xsec} -0.1 0.1
END YODA_ESTIMATE0D_V3
BEGIN YODA_COUNTER_V3 /RAW/_EVTCOUNT
Path: /RAW/_EVTCOUNT
---
{sumw} 10 10
END YODA_COUNTER_V3
""")


class SIDISSignTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)
        self.args = argparse.Namespace(tag="test", dry_run=False, include_diagnostics=False)
        self.measurement = {
            "id": "TEST", "process_kind": "polarized_sidis",
            "analysis": {"name": "TEST"},
            "families": {"nominal": {"helicities": list(sidis.HELICITIES), "mpi": "off", "label": "test"}},
        }
        self.jobs = []
        for component, orders in VALUES.items():
            for order in orders:
                for helicity in sidis.HELICITIES:
                    name = f"{component}-{helicity}-{order}-s0001"
                    self.jobs.append({
                        "id": name, "family": "nominal", "channel": "sidis",
                        "target_component": component, "polarized_pdf_member": 0,
                        "unpolarized_pdf_member": 0, "scale": 1., "mpi": "off",
                        "helicity": helicity, "contribution": order, "order": order,
                        "shard": 1, "events": 10, "status": "success", "output_yoda": name+".yoda",
                    })
        self.set_encoding(False)
        self.manifest = {"jobs": self.jobs, "history": [], "configuration": {}}
        self.groups = campaign._logical_sidis_groups(self.manifest, self.directory, self.measurement)
        self.key = ("nominal", "sidis", 0, 0, 1., "off")
        self.manifest_path = self.directory / experimental.MANIFEST_NAME
        self.manifest_path.write_text(json.dumps(self.manifest))

    def set_encoding(self, signed):
        self.signed = signed
        for job in self.jobs:
            sign = -1 if signed and job["order"] == "NEGNLO" else 1
            metadata(self.directory/job["output_yoda"], sign*10.)

    def audit(self):
        return campaign._audit_sidis_nlo_signs(self.args, self.jobs, self.directory)

    def load(self, jobs, directory, analysis, object_name):
        job = jobs[0]
        value = VALUES[job["target_component"]][job["order"]][sidis.HELICITIES.index(job["helicity"])]
        sign = -1. if self.signed and job["order"] == "NEGNLO" else 1.
        factor = 2. if object_name == "inverse_d" else math.sqrt(2.) if object_name == "covariance" else 1.
        return experimental.BinSeries([0., 1.], [sign*factor*value], [(factor*value)**2])

    def files(self):
        return {str(p.relative_to(self.directory)): p.read_bytes()
                for p in self.directory.rglob("*") if p.is_file()}

    def test_check_is_read_only_and_stale_provenance_is_not_waived(self):
        before = self.files()
        for signed in (False, True):
            self.set_encoding(signed)
            before = self.files()
            with mock.patch.object(campaign, "_campaign_dir", return_value=self.directory), \
                    mock.patch.object(campaign, "_assert_manifest_signature_current", side_effect=campaign.CampaignError("stale")), \
                    contextlib.redirect_stdout(io.StringIO()):
                report = campaign.check_sidis_nlo_signs(self.args, self.measurement)
            self.assertEqual(report["files_checked"], 16)
            self.assertEqual(report["order_coefficients"], {"POSNLO": 1., "NEGNLO": 1. if signed else -1.})
            self.assertFalse(report["postprocessing_compatible"])
            self.assertEqual(before, self.files())

    def test_reject_conflicting_mode_mixed_shards_and_invalid_orders(self):
        self.args.nlo_sign_convention = "signed"
        with self.assertRaisesRegex(campaign.CampaignError, "contradicts"):
            self.audit()
        self.args.nlo_sign_convention = "auto"
        negative = next(job for job in self.jobs if job["order"] == "NEGNLO")
        metadata(self.directory/negative["output_yoda"], -10.)
        with self.assertRaisesRegex(campaign.CampaignError, "Mixed stored"):
            self.audit()
        self.set_encoding(False)
        negative["order"] = "POSNLO"
        with self.assertRaisesRegex(campaign.CampaignError, "Inconsistent"):
            self.audit()
        negative["order"] = "NEGNLO"
        with self.assertRaisesRegex(campaign.CampaignError, "requires POSNLO"):
            campaign._audit_sidis_nlo_signs(self.args, [], self.directory)

    def test_loader_cache_is_not_sign_flipped_in_place(self):
        group = next(iter(self.groups.values()))
        cache = {}
        coefficients = self.audit()["order_coefficients"]
        with mock.patch.object(campaign, "_load_series", side_effect=self.load) as loader:
            first = campaign._load_sidis_object(group, self.directory, "TEST", "ordinary", cache,
                                                order_coefficients=coefficients)
            second = campaign._load_sidis_object(group, self.directory, "TEST", "ordinary", cache,
                                                 order_coefficients=coefficients)
        self.assertEqual(first, second)
        self.assertEqual(loader.call_count, 8)
        self.assertEqual(first["PP"].values, [12.])
        self.assertEqual(first["PP"].variances, [14.**2+2.**2])
        self.assertTrue(all(s.values[0] > 0 for s in cache.values()))

    def test_compass_and_hermes_predictions_match_independent_ratio_errors(self):
        # Both proton and deuteron; same physical samples in either storage convention.
        for target, weights, longitudinal in (("proton", {"P": 1.}, 1.), ("deuteron", {"P": .5, "N": .5}, .925)):
            previous = None
            for signed in (False, True):
                self.set_encoding(signed)
                coefficients = self.audit()["order_coefficients"]
                hermes_snapshot = {"datasets": {"Aparallel": {"target": target, "species": "piplus", "projection": "x"}}}
                compass_measurement = {**self.measurement, "schema_version": 6, "postprocessor": "compass_sidis_a1",
                    "postprocess_config": {"target_outputs": {"T": weights}, "longitudinal_target_scale": longitudinal}}
                compass_snapshot = {"datasets": {"piplus": {"published_target": "T", "points": [{}],
                    "raw_objects": {"ordinary": "ordinary", "inverse_depolarization": "inverse_d", "covariance": "covariance"}}}}
                with mock.patch.object(campaign, "_load_series", side_effect=self.load), \
                        mock.patch("compass_multiplicity_fiducial.unsupported_cells", return_value=[]):
                    hermes = campaign._sidis_prediction_sets(self.groups, self.directory, self.measurement, hermes_snapshot,
                        False, order_coefficients=coefficients)[self.key]["Aparallel"]
                    compass = campaign._compass_sidis_prediction_sets(self.groups, self.directory, compass_measurement,
                        compass_snapshot, order_coefficients=coefficients)[self.key]["piplus"]
                pair = (hermes, compass)
                if previous is not None:
                    self.assertEqual(pair, previous)
                previous = pair
                # Independent derivative propagation over the POS and NEG random yields.
                for result, factor in ((hermes, 1.), (compass, 2.*longitudinal)):
                    numerator = denominator = 0.
                    terms = []
                    for component, weight in weights.items():
                        for order, values in VALUES[component].items():
                            order_sign = -1. if order == "NEGNLO" else 1.
                            for helicity, value in zip(sidis.HELICITIES, values):
                                a = .25*weight*order_sign
                                b = a*sidis.LONGITUDINAL_SIGN[helicity]*factor
                                numerator += b*value
                                denominator += a*value
                                terms.append((a,b,value**2))
                    ratio = numerator/denominator
                    variance = sum(((b-ratio*a)/denominator)**2*v for a,b,v in terms)
                    self.assertAlmostEqual(result["values"][0], ratio)
                    self.assertAlmostEqual(result["errors"][0]**2, variance)

    def test_signed_moments_and_covariance_proxies_are_not_classified_by_bin_sign(self):
        combined = sidis.add_nlo_components({
            "POSNLO": experimental.BinSeries([0.,1.], [-4.], [3.]),
            "NEGNLO": experimental.BinSeries([0.,1.], [-2.], [5.]),
        }, order_coefficients=self.audit()["order_coefficients"])
        self.assertEqual(combined.values, [-2.])
        self.assertEqual(combined.variances, [8.])
        for coefficients in ({}, {"POSNLO":1.,"NEGNLO":0.}, {"POSNLO":-1.,"NEGNLO":-1.}):
            with self.assertRaises(experimental.CampaignError):
                sidis.add_nlo_components({"POSNLO":combined,"NEGNLO":combined}, order_coefficients=coefficients)

    def test_all_postprocessors_stop_before_writing_on_conflicting_convention(self):
        self.args.nlo_sign_convention = "signed"
        before = self.files()
        for postprocess in (campaign.postprocess_sidis, campaign.postprocess_compass_sidis, campaign.postprocess_sidis_diagnostic):
            with mock.patch.object(campaign, "_campaign_dir", return_value=self.directory), \
                    mock.patch.object(campaign, "_assert_manifest_signature_current"), \
                    mock.patch.object(campaign, "_measurement_snapshot") as snapshot, \
                    self.assertRaisesRegex(campaign.CampaignError, "contradicts"):
                postprocess(self.args, self.measurement)
            snapshot.assert_not_called()
            self.assertEqual(before, self.files())

    def test_postprocess_persists_the_audit_and_does_not_modify_generation_configuration(self):
        # Use the real HERMES arithmetic/entry point, replacing only native YODA output.
        snapshot = {"datasets": {"Aparallel": {"target": "proton", "species": "piplus", "projection": "x", "points": []}}}
        with mock.patch.object(campaign, "_campaign_dir", return_value=self.directory), \
                mock.patch.object(campaign, "_assert_manifest_signature_current"), \
                mock.patch.object(campaign, "_measurement_snapshot", return_value=snapshot), \
                mock.patch.object(campaign, "_load_series", side_effect=self.load), \
                mock.patch.object(campaign, "aggregate_uncertainties", return_value={}), \
                mock.patch.object(campaign, "_reference_path", return_value=None), \
                mock.patch.object(campaign, "_estimate_with_bands", return_value=None) as estimate, \
                mock.patch.object(experimental, "_import_yoda", return_value=None), \
                mock.patch.object(experimental, "_write_yoda_objects", side_effect=lambda y,o,p:p.write_text("test output")), \
                contextlib.redirect_stdout(io.StringIO()):
            campaign.postprocess_sidis(self.args, self.measurement)
        summary = json.loads((self.directory/"postprocess/summary.json").read_text())
        manifest = json.loads(self.manifest_path.read_text())
        self.assertEqual(summary["order_sign_audit"]["files_checked"], 16)
        self.assertEqual(manifest["postprocess"]["order_sign_audit"], summary["order_sign_audit"])
        self.assertEqual(manifest["configuration"], {})
        self.assertEqual(manifest["jobs"], self.jobs)
        self.assertEqual(estimate.call_args_list[0].args[3]["NLOCombination"], "normalized POSNLO-NEGNLO bins")

    def test_cli_supports_all_three_campaigns_and_rejects_pp(self):
        for identifier in IDS:
            with mock.patch.object(campaign, "check_sidis_nlo_signs") as check:
                self.assertEqual(campaign.main(["check-nlo-signs", "--measurement", identifier, "--tag", "test"]),0)
                check.assert_called_once()
            args = campaign.make_parser().parse_args([
                "postprocess", "--measurement", identifier, "--tag", "test", "--nlo-sign-convention", "magnitudes"])
            self.assertEqual(args.nlo_sign_convention, "magnitudes")
            handler = "postprocess_sidis" if identifier.startswith("HERMES") else "postprocess_compass_sidis"
            with mock.patch.object(campaign, handler) as postprocess:
                self.assertEqual(campaign.main([
                    "postprocess", "--measurement", identifier, "--tag", "test",
                    "--nlo-sign-convention", "magnitudes",
                ]), 0)
                self.assertEqual(postprocess.call_args.args[0].nlo_sign_convention, "magnitudes")
            with mock.patch.object(campaign, "prepare_pp") as prepare:
                self.assertEqual(campaign.main([
                    "full", "--measurement", identifier, "--tag", "test",
                    "--nlo-sign-convention", "auto", "--dry-run",
                ]), 0)
                self.assertEqual(prepare.call_args.args[0].nlo_sign_convention, "auto")
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(campaign.main(["check-nlo-signs", "--measurement", "PHENIX_2023_I2033856", "--tag", "test"]),2)

    def test_sidis_still_rejects_dis_only_plot_controls(self):
        for options in (["--plot-scale", "2"], ["--hide-scale-uncertainty"], ["--ratiolabel", "test"]):
            with contextlib.redirect_stderr(io.StringIO()) as output, \
                    mock.patch.object(campaign, "prepare_pp") as prepare:
                self.assertEqual(campaign.main([
                    "full", "--measurement", "HERMES_2019_I1698889", "--tag", "test",
                    "--dry-run", "--nlo-sign-convention", "auto", *options,
                ]), 2)
                self.assertIn("only valid for fixed-target DIS", output.getvalue())
                prepare.assert_not_called()


if __name__ == "__main__":
    unittest.main()
