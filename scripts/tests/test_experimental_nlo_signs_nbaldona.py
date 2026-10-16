"""NLO storage checks and independent signed-bin/uncertainty regressions."""
from __future__ import annotations

import contextlib
import copy
import gzip
import io
import json
import math
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import run_experimental_campaign as campaign
from experimental_nlo_signs import audit_order_signs, read_sign_evidence


def write_metadata(path, xsec=10., sumw=10., sumw2=10., entries=10., version=3):
    kind = "YODA_ESTIMATE0D_V3" if version == 3 else "YODA_SCATTER1D_V2"
    text = f"""BEGIN {kind} /_XSEC
Path: /_XSEC
---
{xsec} -0.1 0.1
END {kind}
BEGIN YODA_COUNTER_V{version} /RAW/_EVTCOUNT
Path: /RAW/_EVTCOUNT
---
{sumw} {sumw2} {entries}
END YODA_COUNTER_V{version}
"""
    payload = text.encode()
    path.write_bytes(gzip.compress(payload) if path.suffix == ".gz" else payload)


class SignInspectionTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.jobs = [{"order": order, "output_yoda": f"{order}.yoda"}
                     for order in ("POSNLO", "NEGNLO")]
        for job in self.jobs:
            write_metadata(self.root/job["output_yoda"])

    def audit(self, **options):
        return audit_order_signs(self.jobs, self.root, **options)

    def test_magnitudes_are_subtracted_and_signed_inputs_are_added(self):
        report = self.audit()
        self.assertEqual(report["order_coefficients"], {"POSNLO": 1., "NEGNLO": -1.})
        self.assertEqual(report["files_checked"], 2)
        write_metadata(self.root/"NEGNLO.yoda", xsec=-10., sumw=-10.)
        signed = self.audit()
        self.assertEqual(signed["order_coefficients"], {"POSNLO": 1., "NEGNLO": 1.})
        self.assertNotEqual(report["inputs_sha256"], signed["inputs_sha256"])

    def test_normalized_cross_section_not_raw_weight_sign_controls_arithmetic(self):
        write_metadata(self.root/"NEGNLO.yoda", xsec=10., sumw=-10.)
        self.assertEqual(self.audit()["order_coefficients"]["NEGNLO"], -1.)
        write_metadata(self.root/"NEGNLO.yoda", xsec=-10., sumw=10.)
        self.assertEqual(self.audit()["order_coefficients"]["NEGNLO"], 1.)

    def test_explicit_and_upstream_descriptor_conventions_are_checked(self):
        self.assertEqual(self.audit(mode="magnitudes")["order_coefficients"]["NEGNLO"], -1.)
        self.assertEqual(self.audit(configured={"POSNLO": 1., "NEGNLO": -1.})["resolution"], "configured")
        for options in ({"mode": "signed"}, {"configured": {"POSNLO": 1., "NEGNLO": 1.}},
                        {"mode": "signed", "configured": {"POSNLO": 1., "NEGNLO": -1.}}):
            with self.subTest(options=options), self.assertRaises(ValueError):
                self.audit(**options)

    def test_weighted_or_mixed_event_streams_require_an_explicit_convention(self):
        write_metadata(self.root/"NEGNLO.yoda", sumw=8.)
        with self.assertRaisesRegex(ValueError, "Cannot infer"):
            self.audit()
        self.assertEqual(self.audit(mode="magnitudes")["non_unit_weight_files"], 1)

    def test_zero_missing_nonfinite_and_wrong_sign_metadata_fail(self):
        for values in ({"xsec": 0.}, {"xsec": math.nan}, {"sumw": 0.}, {"entries": 0.}):
            with self.subTest(values=values), self.assertRaises(ValueError):
                write_metadata(self.root/"NEGNLO.yoda", **values)
                self.audit()
        (self.root/"NEGNLO.yoda").write_text("no normalization metadata")
        with self.assertRaisesRegex(ValueError, "Missing or unsupported"):
            self.audit(mode="magnitudes")
        write_metadata(self.root/"NEGNLO.yoda")
        write_metadata(self.root/"POSNLO.yoda", xsec=-10.)
        with self.assertRaisesRegex(ValueError, "POSNLO"):
            self.audit()

    def test_mixed_conventions_between_shards_are_rejected(self):
        write_metadata(self.root/"other.yoda", xsec=-10.)
        self.jobs.append({"order": "NEGNLO", "output_yoda": "other.yoda"})
        with self.assertRaisesRegex(ValueError, "Mixed stored"):
            self.audit()

    def test_v2_gzip_and_lo_only_inputs(self):
        path = self.root/"v2.yoda.gz"
        write_metadata(path, version=2)
        self.assertEqual(read_sign_evidence(path)["cross_section"], 10.)
        self.assertEqual(audit_order_signs([{"order": "LO"}], self.root)["files_checked"], 0)

    def test_read_only_check_reports_stale_provenance_without_waiving_it(self):
        campaign.atomic_write_json(self.root/campaign.MANIFEST_NAME, {"jobs": self.jobs})
        before = {p.name: p.read_bytes() for p in self.root.iterdir()}
        with mock.patch.object(campaign, "campaign_directory", return_value=self.root), \
                mock.patch.object(campaign, "_assert_manifest_signatures_current", side_effect=campaign.CampaignError("stale")), \
                contextlib.redirect_stdout(io.StringIO()):
            report = campaign.check_campaign_nlo_signs(
                types.SimpleNamespace(tag="test"), {"id": "TEST", "combination": {}}
            )
        self.assertFalse(report["generation_signature_current"])
        self.assertEqual(before, {p.name: p.read_bytes() for p in self.root.iterdir()})


SIGNS = {"PP": 1, "PM": -1, "MP": -1, "MM": 1}
VALUES = {
    "P": {"POSNLO": [14., 8., 8., 10.], "NEGNLO": [2., 3., 2., 1.]},
    "N": {"POSNLO": [8., 5., 5., 4.], "NEGNLO": [2., 1., 1., 1.]},
}


class NLOArithmeticTests(unittest.TestCase):
    def setUp(self):
        self.measurement = campaign.get_measurement("E142_1996_I424108")
        self.measurement["combination"]["order_coefficients"] = {"POSNLO": 1., "NEGNLO": -1.}
        self.family = campaign.campaign_family_specs(self.measurement, False)["nominal"]
        self.groups = {}
        self.histograms = {}
        for component, orders in VALUES.items():
            for order, values in orders.items():
                for helicity, value in zip(SIGNS, values):
                    name = f"{component}-{helicity}-{order}.yoda"
                    self.groups[("nominal", 0, 0, 1., component, helicity, order)] = [
                        {"output_yoda": name, "events": 1, "order": order},
                    ]
                    self.histograms[name] = (value, 2. if component == "P" else 3., order)
        self.signed = False

    def read(self, path, histogram):
        value, inverse_d, order = self.histograms[path.name]
        factor = (inverse_d if "SigmaOverD" in histogram else
                  math.sqrt(inverse_d) if "CovarianceProxy" in histogram else 1.)
        sign = -1. if self.signed and order == "NEGNLO" else 1.
        return campaign.BinSeries([.1, .3], [sign*factor*value], [(factor*value)**2])

    def products(self):
        with mock.patch.object(campaign, "read_histogram_series", side_effect=self.read):
            return campaign._load_asymmetry_family_products(
                self.groups, Path("."), self.measurement, "nominal", self.family,
            )

    def test_magnitude_and_signed_samples_give_identical_results_and_errors(self):
        magnitudes = self.products()
        self.signed = True
        self.measurement["combination"]["order_coefficients"]["NEGNLO"] = 1.
        self.assertEqual(self.products(), magnitudes)

    def test_asymmetry_uncertainty_matches_independent_derivatives(self):
        results, _ = self.products()
        for selection, output in self.measurement["outputs"].items():
            uu = output["target_components"]["unpolarized"]
            ll = output["target_components"]["longitudinal"]
            numerator = denominator = 0.
            terms = []
            for component, orders in VALUES.items():
                inverse_d = 2. if component == "P" else 3.
                for order, values in orders.items():
                    order_sign = 1. if order == "POSNLO" else -1.
                    for helicity, value in zip(SIGNS, values):
                        a = uu[component]/4.*order_sign
                        b = ll[component]/4.*SIGNS[helicity]*inverse_d*order_sign
                        denominator += a*value
                        numerator += b*value
                        terms.append((a, b, value**2))
            ratio = numerator/denominator
            variance = sum(((b-ratio*a)/denominator)**2*v for a,b,v in terms)
            self.assertAlmostEqual(results[selection]["a1_values"][0], ratio)
            self.assertAlmostEqual(results[selection]["a1_errors"][0]**2, variance)

    def test_shards_are_averaged_before_order_subtraction_and_moment_sign_is_not_used(self):
        groups = {("nominal", 0, 0, 1., "", "PP", order): [
            {"output_yoda": f"{order}{i}", "events": n} for i,n in enumerate([1,3])
        ] for order in ("POSNLO", "NEGNLO")}
        inputs = {
            "POSNLO0": campaign.BinSeries([0.,1.], [10.], [4.]),
            "POSNLO1": campaign.BinSeries([0.,1.], [20.], [9.]),
            "NEGNLO0": campaign.BinSeries([0.,1.], [-2.], [1.]),
            "NEGNLO1": campaign.BinSeries([0.,1.], [-4.], [4.]),
        }
        with mock.patch.object(campaign, "read_histogram_series", side_effect=lambda p,h: inputs[p.name]):
            series = campaign._load_component_series(groups, Path("."), "TEST", "signed_moment",
                helicities=("PP",), order_coefficients={"POSNLO":1.,"NEGNLO":-1.})["PP"]
        self.assertEqual(series.values, [21.])
        self.assertAlmostEqual(series.variances[0], .25**2*(4.+1.)+.75**2*(9.+4.))

    def test_lo_and_annotation_formula(self):
        self.assertEqual(campaign._order_combination_coefficients(self.measurement, ["LO"]), {"LO":1.})
        annotations = campaign._family_annotations(self.measurement, "nominal", self.family)
        self.assertEqual(annotations["NLOCombination"], "normalized POSNLO-NEGNLO bins")
        self.assertEqual(json.loads(annotations["OrderCoefficients"]), {"POSNLO":1.,"NEGNLO":-1.})

    def test_postprocess_resolves_auto_signs_before_ratios_without_changing_descriptor(self):
        expected, _ = self.products()
        descriptor = copy.deepcopy(self.measurement)
        descriptor["combination"].pop("order_coefficients")
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            for filename in self.histograms:
                write_metadata(directory/filename)
            campaign.atomic_write_json(directory/campaign.MANIFEST_NAME, {"configuration": {}})
            with mock.patch.object(campaign, "campaign_directory", return_value=directory), \
                    mock.patch.object(campaign, "_assert_manifest_signatures_current"), \
                    mock.patch.object(campaign, "_require_complete_matrix", return_value=self.groups), \
                    mock.patch.object(campaign, "_import_yoda", return_value=None), \
                    mock.patch.object(campaign, "read_histogram_series", side_effect=self.read), \
                    mock.patch.object(campaign, "_asymmetry_family_objects", return_value=[]) as writer, \
                    contextlib.redirect_stdout(io.StringIO()) as output:
                campaign.postprocess_campaign(types.SimpleNamespace(tag="test", dry_run=True), descriptor)
            self.assertEqual(writer.call_args.args[3], expected)
            self.assertEqual(writer.call_args.args[5]["NLOCombination"], "normalized POSNLO-NEGNLO bins")
            self.assertIn('"files_checked": 16', output.getvalue())
            self.assertNotIn("order_coefficients", descriptor["combination"])
            self.assertNotIn("_nlo_sign_audit", descriptor)

    def test_unpolarized_comparison_also_subtracts_each_target_order(self):
        family = {**self.family, "helicities": ["00"], "postprocess": "direct_unpolarized"}
        groups = {}
        for component, orders in VALUES.items():
            for order, values in orders.items():
                filename = f"{component}-00-{order}.yoda"
                groups[("unpolarized_nlo",0,0,1.,component,"00",order)] = [
                    {"output_yoda":filename, "events":1}]
                self.histograms[filename] = (sum(values)/4., 1., order)
        with mock.patch.object(campaign, "read_histogram_series", side_effect=self.read):
            results, _ = campaign._load_direct_unpolarized_products(
                groups, Path("."), self.measurement, "unpolarized_nlo", family)
        proton = (sum(VALUES["P"]["POSNLO"])-sum(VALUES["P"]["NEGNLO"]))/4.
        neutron = (sum(VALUES["N"]["POSNLO"])-sum(VALUES["N"]["NEGNLO"]))/4.
        for selection in results:
            self.assertAlmostEqual(results[selection].values[0],
                neutron if selection.startswith("Neutron") else (2*proton+neutron)/3.)

    def test_cli_forwarding_for_both_campaign_entry_points(self):
        import run_phenomenology_campaign as primary
        for command in ("postprocess", "full", "check-nlo-signs"):
            argv = [command, "--measurement", "E142_1996_I424108", "--tag", "test",
                    "--nlo-sign-convention", "magnitudes"]
            args = primary.make_parser().parse_args(argv)
            forwarded = campaign.make_parser().parse_args(primary._legacy_arguments(args))
            self.assertEqual(forwarded.nlo_sign_convention, "magnitudes")
            self.assertEqual(forwarded.command, command)


if __name__ == "__main__":
    unittest.main()
