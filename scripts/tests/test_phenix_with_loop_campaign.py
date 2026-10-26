"""Three-process campaign budgets, hard cuts, and independent normalization."""
import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import run_phenomenology_campaign as campaign
import run_PHENIX_2023_I2033856_WithLoop as wrapper

BASE = "PHENIX_2023_I2033856"
BUDGETS = {
    "photon": [6500000, 4250000, 3250000, 1750000, 1750000],
    "fragmentation": [150000000, 54000000, 34000000, 22500000, 33500000],
    "gammajet_loop": [400000, 200000, 200000, 100000, 100000],
}


class WithLoopCampaignTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        registry = campaign.discover_pp_registry()
        cls.measurement = registry[wrapper.MEASUREMENT]
        cls.base = registry[BASE]
        cls.loop = registry[BASE + "_MEGammaJetLoop-Only"]

    def plan(self, *extra):
        argv = ["full", "--measurement", wrapper.MEASUREMENT, "--tag", "test",
                "--jobs", "63", "--shards", "252", "--seed-base", "426000000",
                "--profile", "central", "--scales", "central",
                "--hard-pt-edges", "4,8,12,16,20,inf", "--photon-enhancement", "12"]
        for channel, counts in BUDGETS.items():
            argv += ["--hard-pt-events", channel + "=" + ",".join(map(str, counts))]
        return campaign._plan(self.measurement, campaign.make_parser().parse_args(argv + list(extra)))

    def test_shared_analysis_reference_and_existing_variants_unchanged(self):
        for old in (self.base, self.loop):
            self.assertEqual(self.measurement["analysis"], old["analysis"])
            self.assertEqual(campaign._measurement_snapshot(self.measurement),
                             campaign._measurement_snapshot(old))
            common = lambda m: (ROOT / m["cards"]["directory"] / m["cards"]["common"]).read_bytes()
            self.assertEqual(common(self.measurement), common(old))
            self.assertEqual(self.measurement["pdf_ensembles"], old["pdf_ensembles"])
        self.assertEqual(set(self.base["channels"]), {"photon", "fragmentation"})
        self.assertEqual(set(self.loop["channels"]), {"photon"})

    def test_exact_budget_four_helicities_open_tail_and_unique_seeds(self):
        plan = self.plan()
        self.assertEqual(plan["logical_jobs"], 60)
        self.assertEqual(plan["shard_jobs"], 15120)
        self.assertEqual(sum(j["events"] for j in plan["jobs"]), 1250000000)
        self.assertEqual(plan["options"]["jobs"], 63)
        self.assertEqual(plan["options"]["hard_pt_ranges_gev"][-1], [20., None])
        for field in ("seed", "id", "output_yoda"):
            self.assertEqual(len({j[field] for j in plan["jobs"]}), 15120)
        self.assertEqual(len({j["run_file"] for j in plan["jobs"]}), 60)
        for channel, counts in BUDGETS.items():
            for helicity in ("PP", "PM", "MP", "MM"):
                for low, count in zip((4, 8, 12, 16, 20), counts):
                    jobs = [j for j in plan["jobs"] if j["channel"] == channel
                            and j["helicity"] == helicity and j["hard_pt_range_gev"][0] == low]
                    self.assertEqual(len(jobs), 252)
                    self.assertEqual(sum(j["events"] for j in jobs), count)

    def test_each_card_has_one_me_correct_cut_and_enhancement_only_on_dijets(self):
        matrix_elements = {"photon": "MEGammaJet", "fragmentation": "MEQCD2to2",
                           "gammajet_loop": "PHENIXGammaJetLoop"}
        jobs = {j["stem"]: j for j in self.plan()["jobs"]}.values()
        self.assertEqual(len(jobs), 60)
        for job in jobs:
            text = campaign._card_text(self.measurement, job)
            insertions = [line for line in text.splitlines()
                          if line.startswith("insert /Herwig/MatrixElements/SubProcess:MatrixElements")]
            self.assertEqual(insertions, ["insert /Herwig/MatrixElements/SubProcess:MatrixElements[0] "
                                         "/Herwig/MatrixElements/" + matrix_elements[job["channel"]]])
            cut = "JetKtCut" if job["channel"] == "fragmentation" else "PhotonKtCut"
            low, high = job["hard_pt_range_gev"]
            self.assertIn(f"{cut}:MinKT {low:g}*GeV", text)
            if high is None:
                self.assertNotIn(f"{cut}:MaxKT", text)
            else:
                self.assertIn(f"{cut}:MaxKT {high:g}*GeV", text)
            if job["channel"] == "fragmentation":
                self.assertIn("QtoQGammaSudakov:EnhancementFactor 12", text)
                self.assertIn("SplittingGenerator:Detuning 2", text)
                self.assertLess(text.index("EnhancementFactor"), text.index("saverun"))
            else:
                self.assertNotIn("EnhancementFactor", text)
            if job["channel"] == "gammajet_loop":
                self.assertIn("PHENIXGammaJetLoop:AmplitudeMode Massive", text)
            self.assertEqual(job["mpi"], "on")
            self.assertEqual(job["shower_spin_correlations"], "on")
            self.assertEqual(job["analysis_instance"], BASE)

    def test_add_all_three_normalized_components_before_asymmetry(self):
        def objects(values, variance):
            return {name: {h: campaign.experimental.BinSeries([6., 7.], [v], [variance])
                           for h, v in zip(("PP", "PM", "MP", "MM"), values)}
                    for name in ("inclusive_cross_section", "isolated_cross_section", "isolated_all")}
        key = lambda c: ("nominal", c, 0, 0, 1., "on")
        samples = {key("photon"): objects([12., 8., 8., 12.], 2.),
                   key("fragmentation"): objects([1., 9., 9., 1.], 3.),
                   key("gammajet_loop"): objects([4., 1., 1., 4.], .5)}
        summed, predictions = campaign._combine_phenix_channels(self.measurement, samples)
        # Per-helicity sums are [17,18,18,17], hence A=-2/70, not a mean of three A's.
        pred = predictions[key("photon")]
        self.assertAlmostEqual(pred["isolated_all"]["values"][0], -2. / 70.)
        self.assertAlmostEqual(pred["inclusive_cross_section"]["values"][0], 17.5)
        self.assertAlmostEqual(pred["inclusive_cross_section"]["errors"][0] ** 2, 5.5 / 4.)
        gradients = [(sign + 2. / 70.) / 70. for sign in (1., -1., -1., 1.)]
        self.assertAlmostEqual(pred["isolated_all"]["errors"][0] ** 2,
                               5.5 * sum(d * d for d in gradients))
        self.assertEqual(summed[key("photon")]["isolated_all"]["PP"].variances, [5.5])
        del samples[key("gammajet_loop")]
        with self.assertRaisesRegex(campaign.CampaignError, "Missing additive"):
            campaign._combine_phenix_channels(self.measurement, samples)

    def test_wrapper_dry_run_uses_separate_directory_without_writes(self):
        with tempfile.TemporaryDirectory() as output, \
                patch.object(campaign, "CAMPAIGN_ROOT", campaign.CAMPAIGN_ROOT):
            stream = io.StringIO()
            with contextlib.redirect_stdout(stream):
                code = wrapper.main(["--output-root", output, "full", "--tag", "test", "--dry-run"])
            self.assertEqual(code, 0)
            self.assertEqual(Path(json.loads(stream.getvalue())["campaign_directory"]),
                             Path(output) / wrapper.MEASUREMENT / "test")
            self.assertEqual(list(Path(output).iterdir()), [])

    def default_plan(self, *extra):
        args = campaign.make_parser().parse_args([
            "full", "--measurement", wrapper.MEASUREMENT, "--tag", "default-test", *extra])
        return campaign._plan(self.measurement, args)

    def test_default_sampling_and_explicit_overrides(self):
        plan = self.default_plan("--shards", "72")
        self.assertEqual(sum(j["events"] for j in plan["jobs"]), 125000000)
        self.assertEqual(plan["shard_jobs"], 4320)
        self.assertEqual(plan["options"]["hard_pt_ranges_gev"][-1], [20., None])
        self.assertEqual(plan["options"]["photon_shower_sampling"],
                         {"channel": "fragmentation", "enhancement_factor": 12., "detuning": 2.})
        budgets = plan["options"]["hard_pt_events_by_channel"]
        self.assertEqual(budgets, {c: [n // 10 for n in ns] for c, ns in BUDGETS.items()})
        uniform = self.default_plan("--lo-events", "600000", "--photon-enhancement", "1")
        self.assertNotIn("hard_pt_events_by_channel", uniform["options"])
        self.assertEqual(sum(j["events"] for j in uniform["jobs"]), 36000000)
        self.assertEqual(uniform["options"]["photon_shower_sampling"]["enhancement_factor"], 1.)
        custom = self.default_plan("--hard-pt-edges", "4,inf", "--lo-events", "1000")
        self.assertEqual(custom["logical_jobs"], 12)
        self.assertEqual(sum(j["events"] for j in custom["jobs"]), 12000)
        partial = self.default_plan("--hard-pt-events", "fragmentation=10,20,30,40,50")
        self.assertEqual(partial["options"]["hard_pt_events_by_channel"],
                         {"fragmentation": [10, 20, 30, 40, 50]})
        self.assertEqual(sum(j["events"] for j in partial["jobs"]), 8000600)

    def test_smoke_is_bounded_and_retains_all_components_and_enhancement(self):
        plan = self.default_plan("--smoke")
        self.assertEqual(plan["logical_jobs"], 60)
        self.assertEqual(sum(j["events"] for j in plan["jobs"]), 6000)
        self.assertNotIn("hard_pt_events_by_channel", plan["options"])
        for job in plan["jobs"]:
            text = campaign._card_text(self.measurement, job)
            self.assertEqual("EnhancementFactor 12" in text, job["channel"] == "fragmentation")
        # The two older variants do not silently inherit the new defaults.
        args = campaign.make_parser().parse_args([
            "full", "--measurement", BASE, "--tag", "old", "--smoke"])
        for measurement, count in ((self.base, 8), (self.loop, 4)):
            legacy = campaign._plan(measurement, args)
            self.assertEqual(legacy["logical_jobs"], count)
            self.assertNotIn("photon_shower_sampling", legacy["options"])
            self.assertNotIn("hard_pt_ranges_gev", legacy["options"])


if __name__ == "__main__":
    unittest.main()
