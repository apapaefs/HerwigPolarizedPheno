"""Contracts for the compensated, particle-level harder-jet pilot."""
import copy
import math
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from test_mc_poljetshapes_campaign import arguments, series, helicities, ROOT
import run_phenomenology_campaign as campaign
import run_mc_poljetshapes_hard_campaign as wrapper
import check_harder_jet_pilot as audit


class HarderJetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.measurement = campaign.discover_pp_registry()["MC_POLJETSHAPES_HARD"]

    def test_descriptor_snapshot_and_windows(self):
        m = self.measurement
        s = campaign._measurement_snapshot(m)
        channel = m["channels"]["combined"]
        self.assertEqual(set(s["observables"]), set(channel["raw_objects"]))
        self.assertEqual(len(channel["angular_observables"]), 32)
        self.assertEqual(len(s["observables"]), 71)
        self.assertEqual(s["selection"]["hardness_windows_gev"], [[20,30],[30,45]])
        self.assertFalse(s["final_state"]["truth_information"])
        self.assertIn("each jet's own window", m["physics"]["hardness_windows"])
        self.assertEqual(wrapper.pinned_arguments(["full", "--tag", "x"])[:3],
                         ["full", "--measurement", m["id"]])
        self.assertEqual(campaign._comparison_pair(m), ("nominal", "lhe_like_shower"))

    def test_compensated_cards_and_single_physics_difference(self):
        m = self.measurement
        args = arguments(tag="hard-pilot-test", families="nominal,lhe_like_shower")
        plan = campaign._plan(m, args)
        common = (ROOT / m["cards"]["directory"] / m["cards"]["common"]).read_text()
        self.assertIn("SubProcess:Preweights 0 /Herwig/Weights/HardJetSampling", common)
        self.assertNotIn("SubProcess:Reweights", common)
        self.assertIn("HardJetSampling:Power 4", common)
        self.assertIn("HardJetSampling:Scale 20*GeV", common)
        self.assertIn("JetKtCut:MinKT 3.0*GeV", common)
        for h, polarizations in m["cards"]["helicities"].items():
            jobs = [j for j in plan["jobs"] if j["helicity"] == h]
            a, b = [campaign._card_text(m, job) for job in jobs]
            self.assertIn(f"FirstLongitudinalPolarization {polarizations[0]}", a)
            self.assertIn(f"SecondLongitudinalPolarization {polarizations[1]}", a)
            self.assertIn("SpinCorrelations Yes", a)
            self.assertEqual(a, b.replace(jobs[1]["stem"], jobs[0]["stem"])
                             .replace("HardProcessSpin No", "HardProcessSpin Yes"))
            self.assertEqual(set(jobs[0]["analysis_instances"]),
                             {"MC_POLJETSHAPES_HARD", "MC_POLDIJETS"})
        config = campaign._manifest_configuration(m, args, plan)
        changed = copy.deepcopy(m)
        changed["sampling"]["power"] = 2
        self.assertNotEqual(config, campaign._manifest_configuration(changed, args, plan))

    def test_windowed_rate_identities_and_missing_inputs(self):
        groups = self.measurement["channels"]["combined"]["rate_groups"]
        inputs = {}
        for group in groups.values():
            for role, values in (("denominator", [100,100]), ("ge3", [40,20]), ("ge4", [10,5])):
                inputs[group[role]] = helicities(series(values, values))
        predictions = campaign._mc_poljetshapes_prediction(inputs, (), rate_groups=groups)
        for w in groups:
            self.assertEqual(predictions["R32_UU_"+w]["values"], [.4,.2])
            self.assertEqual(predictions["R43_UU_"+w]["values"], [.25,.25])
            self.assertEqual(predictions["ThirdJetVeto_UU_"+w]["values"], [.6,.8])
        del inputs[groups["pt20_30"]["ge4"]]
        with self.assertRaisesRegex(campaign.CampaignError, "Missing rate inputs"):
            campaign._mc_poljetshapes_prediction(inputs, (), rate_groups=groups)

    def test_weighted_moments_neff_and_shape_closure(self):
        edges = [-math.pi+i*math.pi/12 for i in range(25)]
        weights = [1+i%3 for i in range(24)]
        sample = series(weights, [w*w for w in weights], edges)
        expected = sum(weights)**2/sum(w*w for w in weights)
        self.assertAlmostEqual(campaign._effective_entries(sample), expected)
        outputs = campaign._mc_poljetshapes_prediction({"angle": helicities(sample)}, ["angle"])
        shape = outputs["ShapeUU_angle"]
        self.assertAlmostEqual(sum(v*(b-a) for v,a,b in zip(shape["values"], edges, edges[1:])), 1)
        self.assertGreater(outputs["A2UU_angle"]["errors"][0], 0)
        self.assertAlmostEqual(outputs["A2LL_angle"]["values"][0], 0)

    def test_normalization_comparison_retains_both_errors(self):
        outputs = {f"Sigma{h}_jet{j}_pt": {"edges": [20,30,45],
                   "values": [100.,200.], "errors": [5.,10.]}
                   for h in ("UU", "PP", "PM", "MP", "MM") for j in (1,2)}
        pilot = {"variations": {"nominal": outputs}}
        reference = copy.deepcopy(pilot)
        rows = audit.normalization_comparison(pilot, reference)
        self.assertEqual(len(rows), 20)
        self.assertTrue(all(row["ratio"] == 1 for row in rows))
        self.assertAlmostEqual(rows[0]["independent_error"], math.sqrt(2)*0.05)

    def test_compensated_preflight_checks_loaded_sampler(self):
        runtime = {"tools": {"Herwig": "/test/bin/Herwig"},
                   "herwig_prefix": "/test", "hwshower_library": "/test/HwShower.so",
                   "hwmehadron_library": "/test/HwMEHadron.so"}
        output = "calling init: /test/HwShower.so\ncalling init: /test/HwMEHadron.so\n"
        with patch.object(campaign.subprocess, "run") as run, patch.object(
                campaign.experimental, "_find_runtime_library", return_value=Path("/test/ReweightMinPT.so")):
            run.return_value = subprocess.CompletedProcess([], 0, "No\n", output)
            with self.assertRaisesRegex(campaign.CampaignError, "ReweightMinPT"):
                campaign._preflight_hard_process_spin(self.measurement, runtime)
            run.return_value = subprocess.CompletedProcess([], 0, "No\n", output+"calling init: /test/ReweightMinPT.so\n")
            campaign._preflight_hard_process_spin(self.measurement, runtime)

    def test_missing_baseline_cannot_pass_projection(self):
        m = self.measurement
        name = m["statistics_projection"]["baseline_observables"][0]
        sample = series([100]*24, [.001]*24, [-math.pi+i*math.pi/12 for i in range(25)])
        raw = {name: helicities(sample)}
        prediction = campaign._mc_poljetshapes_prediction(raw, [name])
        keys = [(family, "combined", 0, 0, 1.0, "off") for family in m["comparison_pair"]]
        result, _, _ = campaign._mc_poljetshapes_assessment(
            m, {key: prediction for key in keys}, {key: raw for key in keys},
            {"configuration": {"families": m["comparison_pair"], "lo_events": 500000}})
        self.assertTrue(all(not tier["passes"] for tier in result["tiers"]))

    def test_focus_has_harder_spectra_and_ratios(self):
        m = self.measurement
        stems = {p[0] for section in m["focus_sections"] for p in section["plots"]}
        for w in ("pt20_30", "pt30_45"):
            for j in (3,4):
                self.assertIn(f"SigmaUU_{w}_jet{j}_pt", stems)
            self.assertIn("A2LL_"+w+"_dpsi12_j2_kt1", stems)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "index.html").write_text("<h1>Pilot</h1>")
            for stem in stems:
                path = root / m["analysis"]["name"] / (stem+".png")
                path.parent.mkdir(exist_ok=True)
                path.write_bytes(b"fixture")
            page = campaign._write_mc_poljetshapes_focus_index(root, m).read_text(encoding="utf-8")
            self.assertIn("20–30 GeV", page)
            self.assertIn("Full spin / LHE-like", page)
            self.assertIn("ordinary shower-generated spin correlations", page)

    def test_geometry_matches_frozen_analysis(self):
        frozen = (ROOT / "analyses/rivet/pp/MC_POLJETSHAPES.cc").read_text()
        new = (ROOT / "analyses/rivet/pp/PolJetShapesHard.hh").read_text()
        start = frozen.index("    constexpr double EPS")
        end = frozen.index("    struct JetAngles")
        self.assertIn(frozen[start:end], new)
        self.assertIn("pt >= 20.0 && pt < 30.0", new)
        self.assertIn("pt >= 30.0 && pt < 45.0", new)
        source = (ROOT / "analyses/rivet/pp/MC_POLJETSHAPES_HARD.cc").read_text()
        self.assertNotRegex(source, r"genParticle\(|ancestors\(|partons\(")
        self.assertIn("ptWindow(jets[j].pT()/GeV)", source)

    @unittest.skipUnless(shutil.which("rivet-config"), "Rivet environment not loaded")
    def test_compiled_geometry_fixture(self):
        import shlex
        compiler = shutil.which("g++-16") or shutil.which("g++")
        flags = shlex.split(subprocess.check_output(["rivet-config", "--cppflags", "--ldflags", "--libs"], text=True))
        with tempfile.TemporaryDirectory() as directory:
            binary = Path(directory)/"geometry"
            subprocess.run([compiler, "-std=c++17", "-I"+str(ROOT/"analyses/rivet/pp"),
                            str(ROOT/"scripts/tests/fixtures/harder_jet_geometry.cc"),
                            *flags, "-o", str(binary)], check=True, capture_output=True)
            subprocess.run([str(binary)], check=True)


if __name__ == "__main__":
    unittest.main()
