"""Hard-spin comparison contracts, separate from the legacy azimuthal control."""
import copy
import json
import io
import os
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from test_mc_poljetshapes_campaign import arguments, series, helicities
import run_phenomenology_campaign as campaign
import run_mc_poljetshapes_lhe_campaign as wrapper


class HardProcessSpinCampaignTests(unittest.TestCase):
    def setUp(self):
        self.measurement = campaign.discover_pp_registry()["MC_POLJETSHAPES_LHE"]

    def test_descriptor_and_snapshot(self):
        m = self.measurement
        self.assertEqual(campaign._comparison_pair(m), ("nominal", "lhe_like_shower"))
        self.assertEqual(campaign._measurement_snapshot(m)["measurement"], m["id"])
        legacy = campaign.discover_pp_registry()["MC_POLJETSHAPES"]
        self.assertEqual(m["channels"], legacy["channels"])
        self.assertEqual(m["analysis"], legacy["analysis"])
        self.assertNotIn("comparison_pair", legacy)
        self.assertEqual(wrapper.pinned_arguments(["plot", "--tag", "x"])[:3],
                         ["plot", "--measurement", "MC_POLJETSHAPES_LHE"])

    def test_cards_preserve_polarized_hard_process(self):
        m = self.measurement
        args = arguments(tag="hard-spin-test", families="nominal,lhe_like_shower")
        plan = campaign._plan(m, args)
        jobs = [j for j in plan["jobs"] if j["helicity"] == "PP"]
        self.assertEqual(len(jobs), 2)
        cards = [campaign._card_text(m, job) for job in jobs]
        for card in cards:
            self.assertIn("ShowerHandler:SpinCorrelations Yes", card)
            self.assertIn("FirstLongitudinalPolarization 1", card)
            self.assertIn("read MC_POLJETSHAPES-Common.in", card)
        self.assertEqual(set(jobs[0]["analysis_instances"]), {"MC_POLJETSHAPES", "MC_POLDIJETS"})
        self.assertIn("HardProcessSpin Yes", cards[0])
        self.assertIn("HardProcessSpin No", cards[1])
        normalized = cards[1].replace(jobs[1]["stem"], jobs[0]["stem"])
        normalized = normalized.replace("HardProcessSpin No", "HardProcessSpin Yes")
        self.assertEqual(cards[0], normalized)
        config = campaign._manifest_configuration(m, args, plan)
        self.assertEqual(config["shower_spin_policy"]["lhe_like_shower"],
                         {"hard_process_spin": "off", "shower_spin_correlations": "on"})
        self.assertEqual(jobs[1]["hard_process_spin"], "off")
        changed = copy.deepcopy(m)
        changed["families"]["lhe_like_shower"]["hard_process_spin"] = "on"
        self.assertNotEqual(config, campaign._manifest_configuration(changed, args, plan))

    def test_new_control_ratio_and_sparse_support(self):
        m = self.measurement
        predictions = {}
        for family, value in (("nominal", 20), ("lhe_like_shower", 10)):
            key = campaign._variation_id((family, "combined", 0, 0, 1.0, "off"))
            predictions[key] = {"SigmaUU_jet3_pt": {
                "edges": [0, 1, 2], "values": [value, 0.1], "errors": [1, 1]}}
        ratios = campaign._mc_poljetshapes_plot_ratios(m, {"variations": predictions})
        self.assertEqual(ratios["SigmaUU_jet3_pt"]["values"], [2.0, None])
        key = ("nominal", "combined", 0, 0, 1.0, "off")
        outputs = {key: {"ALL_jet3_pt": {
            "edges": [0, 1, 2], "values": [0.0, 1.0], "errors": [0.1, 1.e-8]}}}
        raw = {key: {"jet3_pt": helicities(series([100, 1], [100, 1]))}}
        campaign._mask_spin_asymmetry_support(m, outputs, raw)
        self.assertEqual(outputs[key]["ALL_jet3_pt"]["values"], [0.0, None])

    def test_preflight_rejects_missing_switch_and_accepts_both_values(self):
        runtime = {"tools": {"Herwig": "/test/bin/Herwig"},
                   "hwshower_library": "/test/lib/HwShower.so",
                   "hwmehadron_library": "/test/lib/HwMEHadron.so"}
        with patch.object(campaign.subprocess, "run") as run:
            run.return_value = subprocess.CompletedProcess([], 1, "", "Unknown interface")
            with self.assertRaisesRegex(campaign.CampaignError, "HardProcessSpin"):
                campaign._preflight_hard_process_spin(self.measurement, runtime)
            run.return_value = subprocess.CompletedProcess([], 0, "0 [No]\n", "")
            with self.assertRaisesRegex(campaign.CampaignError, "actually loaded"):
                campaign._preflight_hard_process_spin(self.measurement, runtime)
            run.return_value = subprocess.CompletedProcess([], 0, "0 [No]\n",
                "calling init: /test/lib/HwShower.so\ncalling init: /test/lib/HwMEHadron.so\n")
            campaign._preflight_hard_process_spin(self.measurement, runtime)

    def test_conditional_moments_mask_each_fraction_bin_separately(self):
        measurement = copy.deepcopy(self.measurement)
        measurement["channels"]["combined"]["conditional_angular_observables"] = {
            "angle_vs_fraction": {"raw_observable": "flat_angle", "angle_bins": 4}
        }
        key = ("nominal", "combined", 0, 0, 1.0, "off")
        outputs = {key: {"C2LL_angle_vs_fraction": {
            "edges": [0, 1, 2], "values": [0.0, 1.0], "errors": [0.1, 0.1]}}}
        raw = {key: {"flat_angle": helicities(series([100]*4+[1]*4, [100]*4+[1]*4))}}
        campaign._mask_spin_asymmetry_support(measurement, outputs, raw)
        self.assertEqual(outputs[key]["C2LL_angle_vs_fraction"]["values"], [0.0, None])
        self.assertEqual(outputs[key]["C2LL_angle_vs_fraction"]["support_mask"], [True, False])

    def test_comparison_pair_validation(self):
        for pair in (["nominal", "nominal"], ["nominal", "missing"], ["bad"]):
            m = copy.deepcopy(self.measurement)
            m["comparison_pair"] = pair
            with self.assertRaises(campaign.CampaignError):
                campaign._comparison_pair(m)

    def test_incompatible_manifest_is_rejected_before_campaign_writes(self):
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory)
            manifest = destination / campaign.experimental.MANIFEST_NAME
            manifest.write_text(json.dumps({"configuration": {"legacy": True}}))
            original = manifest.read_bytes()
            args = arguments(tag="immutable-test",families="nominal,lhe_like_shower",dry_run=False)
            with patch.object(campaign,"_runtime",return_value={}), \
                 patch.object(campaign,"_preflight_hard_process_spin"), \
                 patch.object(campaign,"_campaign_dir",return_value=destination):
                with self.assertRaisesRegex(campaign.CampaignError,"incompatible manifest"):
                    campaign.prepare_pp(args,self.measurement)
            self.assertEqual(manifest.read_bytes(),original)
            self.assertEqual(list(destination.iterdir()),[manifest])

    def test_focus_labels_do_not_claim_internal_correlations_are_off(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "index.html").write_text("<h1>Plots</h1>")
            (root / "MC_POLJETSHAPES").mkdir()
            (root / "MC_POLJETSHAPES" / "SigmaUU_jet3_pt.png").write_bytes(b"fixture")
            page = campaign._write_mc_poljetshapes_focus_index(root, self.measurement)
            text = page.read_text()
            self.assertIn("full-spin showering and blue is LHE-like showering", text)
            self.assertIn("ordinary shower-generated spin correlations", text)
            self.assertIn("closure validation", text)

    def test_companion_spectra_receive_particle_level_axis_labels(self):
        with tempfile.TemporaryDirectory() as directory:
            script=Path(directory)/"SigmaUU_jet3_pt.py"
            script.write_text("ax.set_xlabel(ax_xLabel)\n")
            campaign._label_companion_jet_spectrum(script)
            self.assertIn("Jet 3 $p_T$ [GeV]",script.read_text())
            self.assertNotIn("parton",script.read_text())

    def test_parallel_renderer_preserves_order_and_reports_failures(self):
        with tempfile.TemporaryDirectory() as directory:
            paths=[Path(directory)/f"plot-{i}.py" for i in range(3)]
            for i,path in enumerate(paths):
                path.write_text(f"print('plot {i}')\n")
            log=io.StringIO()
            campaign._render_parallel_plot_scripts(paths,os.environ,log,2)
            self.assertLess(log.getvalue().index("plot 0"),log.getvalue().index("plot 2"))
            paths[1].write_text("raise SystemExit(3)\n")
            with self.assertRaisesRegex(campaign.CampaignError,"status 3"):
                campaign._render_parallel_plot_scripts(paths,os.environ,io.StringIO(),2)
            with self.assertRaisesRegex(campaign.CampaignError,"between 1 and 100"):
                campaign._render_parallel_plot_scripts(paths,os.environ,io.StringIO(),0)

    def test_plot_worker_option_is_presentation_only(self):
        parsed=campaign.make_parser().parse_args([
            "plot","--measurement","MC_POLJETSHAPES_LHE","--tag","x","--plot-jobs","16"])
        self.assertEqual(parsed.plot_jobs,16)


if __name__ == "__main__":
    unittest.main()
