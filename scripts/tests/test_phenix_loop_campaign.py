"""Loop-only inputs must preserve the expanded PHENIX observable contract."""
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
import run_PHENIX_2023_I2033856_MEGammaJetLoop_Only as wrapper

BASE = "PHENIX_2023_I2033856"
BUDGETS = [1860000, 585000, 435000, 397500, 397500]


class LoopCampaignTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        registry = campaign.discover_pp_registry()
        cls.base = registry[BASE]
        cls.loop = registry[wrapper.MEASUREMENT]

    def plan(self, *extra):
        args = campaign.make_parser().parse_args([
            "full", "--measurement", wrapper.MEASUREMENT, "--tag", "test",
            "--jobs", "96", "--shards", "96", "--seed-base", "226000000",
            "--hard-pt-edges", "4,8,12,16,20,inf", "--hard-pt-events",
            "photon=" + ",".join(map(str, BUDGETS)), "--dry-run", *extra,
        ])
        return campaign._plan(self.loop, args)

    def test_same_analysis_reference_and_common_card(self):
        self.assertEqual(self.loop["analysis"], self.base["analysis"])
        for key in ("pdf_ensembles", "scales", "global_uncertainties"):
            self.assertEqual(self.loop[key], self.base[key])
        self.assertEqual(campaign._measurement_snapshot(self.loop),
                         campaign._measurement_snapshot(self.base))
        common = lambda m: (ROOT / m["cards"]["directory"] / m["cards"]["common"]).read_bytes()
        self.assertEqual(common(self.loop), common(self.base))
        for key in ("fiducial", "isolation", "counting", "sqrt_s_gev"):
            self.assertEqual(self.loop["physics"][key], self.base["physics"][key])
        self.assertNotEqual(campaign._signature(self.loop), campaign._signature(self.base))

    def test_budget_is_14p7_million_across_four_helicities_and_five_slices(self):
        plan = self.plan()
        self.assertEqual(plan["logical_jobs"], 20)
        self.assertEqual(plan["shard_jobs"], 1920)
        self.assertEqual(sum(j["events"] for j in plan["jobs"]), 14700000)
        self.assertEqual(set(self.loop["channels"]), {"photon"})
        for h in ("PP", "PM", "MP", "MM"):
            for low, count in zip((4, 8, 12, 16, 20), BUDGETS):
                jobs = [j for j in plan["jobs"] if j["helicity"] == h
                        and j["hard_pt_range_gev"][0] == low]
                self.assertEqual(sum(j["events"] for j in jobs), count)
                self.assertEqual(len(jobs), 96)
        for key in ("seed", "id", "output_yoda"):
            self.assertEqual(len({j[key] for j in plan["jobs"]}), 1920)
        self.assertEqual([j["seed"] for j in plan["jobs"]],
                         list(range(226000000, 226001920)))

    def test_only_loop_matrix_element_with_photon_slices_and_full_spin(self):
        jobs = {j["stem"]: j for j in self.plan()["jobs"]}.values()
        for job in jobs:
            card = campaign._card_text(self.loop, job)
            self.assertIn("create Herwig::MEPP2GammaJetLoop ", card)
            self.assertIn("PHENIXGammaJetLoop:AmplitudeMode Massive", card)
            insertions = [line for line in card.splitlines()
                          if line.startswith("insert /Herwig/MatrixElements/SubProcess:MatrixElements")]
            self.assertEqual(insertions, [
                "insert /Herwig/MatrixElements/SubProcess:MatrixElements[0] "
                "/Herwig/MatrixElements/PHENIXGammaJetLoop"])
            low, high = job["hard_pt_range_gev"]
            self.assertIn(f"PhotonKtCut:MinKT {low:g}*GeV", card)
            if high is None:
                self.assertNotIn("PhotonKtCut:MaxKT", card)
            else:
                self.assertIn(f"PhotonKtCut:MaxKT {high:g}*GeV", card)
            self.assertNotIn("EnhancementFactor", card)
            self.assertNotIn("Detuning", card)
            self.assertEqual(job["mpi"], "on")
            self.assertEqual(job["shower_spin_correlations"], "on")
            self.assertEqual(job["analysis_instance"], BASE)
            for member in ("polarized_pdf_member", "unpolarized_pdf_member"):
                self.assertEqual(job[member], 0)
            for particle, value in zip(("First", "Second"),
                                      self.loop["cards"]["helicities"][job["helicity"]]):
                self.assertIn(f"{particle}LongitudinalPolarization {value}", card)

    def test_scale_variations_address_loop_object(self):
        jobs = self.plan("--scales", "all")["jobs"]
        job = next(j for j in jobs if j["scale"] == 2)
        self.assertIn("PHENIXGammaJetLoop:ScalePreFactor 4",
                      campaign._card_text(self.loop, job))

    def test_fragmentation_enhancement_is_rejected_for_loop_only(self):
        with self.assertRaisesRegex(campaign.CampaignError, "fragmentation"):
            self.plan("--photon-enhancement", "12")

    def test_wrapper_keeps_new_directory_and_dry_run_writes_nothing(self):
        with tempfile.TemporaryDirectory() as output, \
                patch.object(campaign, "CAMPAIGN_ROOT", campaign.CAMPAIGN_ROOT):
            stream = io.StringIO()
            with contextlib.redirect_stdout(stream):
                code = wrapper.main(["--output-root", output, "full", "--tag", "loop",
                                     "--dry-run"])
            self.assertEqual(code, 0)
            plan = json.loads(stream.getvalue())
            self.assertEqual(Path(plan["campaign_directory"]),
                             Path(output) / wrapper.MEASUREMENT / "loop")
            self.assertEqual(list(Path(output).iterdir()), [])


if __name__ == "__main__":
    unittest.main()
