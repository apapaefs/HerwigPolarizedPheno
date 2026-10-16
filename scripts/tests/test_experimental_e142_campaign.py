from __future__ import annotations

import contextlib
import copy
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import run_experimental_campaign as campaign
import phenomenology_reference_data as reference
import run_E142_1996_I424108 as wrapper
import check_e142_window_sensitivity as sensitivity

ID = "E142_1996_I424108"


class E142Tests(unittest.TestCase):
    def setUp(self):
        self.measurement = campaign.get_measurement(ID)
        self.snapshot = reference.validate_vendored(ID)

    def test_reference_values_are_A1_not_g1_or_A2(self):
        helium, neutron = self.snapshot["datasets"]
        self.assertEqual([p["value"] for p in helium["points"]],
                         [-.0264, -.0238, -.0317, -.0447, -.0463, -.0333, -.0003, -.0145])
        self.assertEqual([p["value"] for p in neutron["points"]],
                         [-.092, -.082, -.109, -.162, -.170, -.113, .050, .006])
        self.assertEqual([p["stat"] for p in neutron["points"]],
                         [.061, .038, .031, .033, .041, .044, .083, .159])
        self.assertEqual([p["systematic_combined"] for p in helium["points"]],
                         [.0054, .0039, .0048, .0068, .0094, .0121, .022, .0199])
        self.assertEqual(helium["provenance"]["table_doi"], "10.17182/hepdata.22340.v1/t1")
        self.assertEqual(neutron["provenance"]["table_doi"], "10.17182/hepdata.22340.v1/t3")

    def test_reference_plot_coordinates_use_published_x_for_both_targets_and_aliases(self):
        for dataset in self.snapshot["datasets"]:
            for path in [dataset["rivet_path"], *dataset["rivet_aliases"]]:
                points = campaign._fixed_reference_overlay_points(self.measurement, self.snapshot, Path(path).name)
                self.assertEqual([p["plot_x"] for p in points], [p["x_mean"] for p in dataset["points"]])
        self.assertEqual(self.snapshot["datasets"][0]["points"][4]["x_mean"], .175)
        self.assertEqual(self.snapshot["datasets"][1]["points"][4]["x_mean"], .174)

    def test_legacy_multidataset_reference_still_plots_Q2(self):
        data = reference.validate_vendored("HERMES_2007_I726689_LEGACY")
        dataset = data["datasets"][0]
        points = campaign._fixed_reference_overlay_points({}, data, Path(dataset["rivet_path"]).name)
        self.assertEqual(points[0]["plot_x"], dataset["points"][0]["q2_mean"])

    def test_all_outputs_share_sixteen_logical_runs_and_have_distinct_paths(self):
        jobs = campaign.build_job_matrix(self.measurement, 3000000, 300000, 100, 424108)
        self.assertEqual(len(jobs), 1600)
        self.assertEqual(sum(job["events"] for job in jobs), 26400000)
        self.assertEqual(len({job["seed"] for job in jobs}), len(jobs))
        self.assertEqual(len({job["stem"] for job in jobs}), 16)
        self.assertEqual(len(self.measurement["outputs"]), 6)
        paths = [output[key] for output in self.measurement["outputs"].values()
                 for key in ("a1", "apar", "parity_pp_mm", "parity_pm_mp")]
        self.assertEqual(len(paths), len(set(paths)))
        for window in ("W025", "W050", "W100"):
            self.assertEqual(self.measurement["raw_observables"][f"He3_{window}"],
                             self.measurement["raw_observables"][f"Neutron_{window}"])

    def test_every_planned_card_exists_and_nominal_is_polarized_NLO(self):
        jobs = campaign.build_job_matrix(self.measurement, 1000, 100, 1, 10, True)
        self.assertEqual(len(jobs), 60)
        common = (ROOT/self.measurement["cards"]["directory"]/self.measurement["cards"]["common"]).read_text()
        self.assertIn("BeamA /Herwig/Particles/e-", common)
        self.assertIn("BeamEMaxA 25.51*GeV", common)
        self.assertIn("SubProcess:MatrixElements[0] /Herwig/MatrixElements/PowhegMEDISNCPol", common)
        self.assertIn("PowhegMEDISNCPol:UsePOWHEGRealSpinVertex Yes", common)
        self.assertIn("ShowerHandler:SpinCorrelations Yes", common)
        for job in jobs:
            text = (ROOT/job["card_source"]).read_text()
            self.assertIn("saverun " + job["stem"] + " EventGenerator", text)
            particle = "p+" if job["component"] == "P" else "n0"
            self.assertIn(f"TargetParticle /Herwig/Particles/{particle}", text)
            if job["family"] == "nominal":
                contribution = "PositiveNLO" if job["order"] == "POSNLO" else "NegativeNLO"
                self.assertIn("PowhegMEDISNCPol:Contribution " + contribution, text)

    def test_actual_target_coefficients_give_nuclear_ratio_and_free_neutron(self):
        family = campaign.campaign_family_specs(self.measurement, False)["nominal"]
        # Independent whole-nucleus expectation for sigma_U=(10,20),
        # sigma_LL/D=(4,-6). Overall per-nucleon normalization cancels.
        expected = {"He3": (2*(-.027)*4+.87*(-6))/(2*10+20), "Neutron": -6/20}
        for selection in self.measurement["outputs"]:
            uu = campaign._target_component_coefficients(self.measurement, family, "unpolarized", selection)
            ll = campaign._target_component_coefficients(self.measurement, family, "longitudinal", selection)
            ratio = (ll["P"]*4+ll["N"]*(-6))/(uu["P"]*10+uu["N"]*20)
            self.assertAlmostEqual(ratio, expected[selection.split("_")[0]])

    def test_yoda_reference_has_both_targets_and_window_aliases(self):
        try:
            import yoda
        except ImportError:
            self.skipTest("YODA is not active")
        objects = yoda.read(str(ROOT/self.measurement["analysis"]["reference_yoda"]))
        self.assertEqual(len(objects), 6)
        for dataset in self.snapshot["datasets"]:
            for path in [dataset["rivet_path"], *dataset["rivet_aliases"]]:
                obj = objects["/REF"+path]
                self.assertEqual(list(obj.xEdges()), dataset["bin_edges"])
                for bin_object, point in zip(obj.bins(), dataset["points"]):
                    self.assertAlmostEqual(bin_object.val(), point["value"])
                    self.assertAlmostEqual(bin_object.errAvg("stat"), point["stat"])
                    self.assertAlmostEqual(bin_object.errAvg("published_total"), point["systematic_combined"])

    def test_wrapper_dry_run_writes_nothing_and_uses_personal_output_root(self):
        with tempfile.TemporaryDirectory() as directory:
            with mock.patch.object(campaign, "CAMPAIGN_ROOT", Path(directory)), contextlib.redirect_stdout(io.StringIO()) as stream:
                status = wrapper.main(["--output-root", directory, "full", "--tag", "check", "--dry-run"])
            self.assertEqual(status, 0)
            plan = json.loads(stream.getvalue())
            self.assertEqual(plan["campaign_directory"], str(Path(directory)/ID/"check"))
            self.assertEqual(list(Path(directory).iterdir()), [])

    def test_postprocessing_roundtrip_and_sensitivity_report(self):
        try:
            import yoda
        except ImportError:
            self.skipTest("YODA is not active")
        family = campaign.campaign_family_specs(self.measurement, False)["nominal"]
        edges = self.snapshot["datasets"][0]["bin_edges"]
        def load(groups, directory, analysis, name, family_id, component, helicities, orders, *args):
            u, ll = (10., 2.) if component == "P" else (20., -3.)
            factor, variance_factor = 1., 1.
            if name.startswith("SigmaOverD"):
                factor, variance_factor = 2., 4.
            elif name.startswith("CovarianceProxy"):
                factor, variance_factor = 2**.5, 2.
            return {h: campaign.BinSeries(edges, [(u+(ll if h in ("PP","MM") else -ll))*factor]*8,
                                           [.01*variance_factor]*8) for h in helicities}
        with mock.patch.object(campaign, "_load_component_series", side_effect=load):
            results, diagnostics = campaign._load_asymmetry_family_products(
                {}, Path("."), self.measurement, "nominal", family)
        objects = campaign._asymmetry_family_objects(
            yoda, self.measurement, self.snapshot, results, diagnostics,
            campaign._family_annotations(self.measurement, "nominal", family))
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary)/"fixture.yoda"
            campaign._write_yoda_objects(yoda, objects, destination)
            loaded = yoda.read(str(destination))
            for selection, output in self.measurement["outputs"].items():
                obj = loaded[output["a1"]]
                neutron = selection.startswith("Neutron")
                self.assertAlmostEqual(obj.bin(1).val(), -.3 if neutron else -.1359)
                self.assertEqual(obj.annotation("Observable"), "A1n" if neutron else "A1He3")
                combination = obj.annotation("TargetComponentCombination")
                if isinstance(combination, str):
                    combination = json.loads(combination)
                self.assertEqual(combination["unpolarized"]["P"], 0. if neutron else 2/3)
            csv_path = Path(temporary)/"summary.csv"
            campaign._write_summary_csv(csv_path, [row for selection,result in results.items()
                                                   for row in campaign._summary_rows(selection,result)])
            report = sensitivity.compare(csv_path)
            self.assertEqual(len(report["records"]), 32)
            self.assertTrue(all(record["delta_a1"] == 0 for record in report["records"]))

    def test_publication_mean_disagreement_fails(self):
        tables = [json.loads((ROOT/source["raw"]).read_text()) for source in reference.SOURCES[ID]]
        tables[0]["values"][0]["x"][0]["value"] = "0.039"
        with self.assertRaises(reference.ReferenceDataError):
            reference.normalize_e142(*tables)


if __name__ == "__main__":
    unittest.main()
