"""Photon sampling controls must be local, recorded, and safe to resume."""
import contextlib
import copy
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import run_experimental_campaign as experimental
import run_phenomenology_campaign as campaign
import run_PHENIX_2023_I2033856 as wrapper

ID = "PHENIX_2023_I2033856"


class PhotonEnhancementTests(unittest.TestCase):
    def setUp(self):
        self.measurement = campaign.discover_pp_registry()[ID]

    def plan(self, *extra):
        args = campaign.make_parser().parse_args([
            "full", "--measurement", ID, "--tag", "sampling-test",
            "--lo-events", "10", "--shards", "3", "--jobs", "2",
            "--seed-base", "91000000", "--hard-pt-edges", "4,8,12,16,20,inf",
            "--dry-run", *extra,
        ])
        return campaign._plan(self.measurement, args), args

    def test_factors_change_only_fragmentation_cards_before_saverun(self):
        original, _ = self.plan()
        for factor in (1, 2, 5, 10):
            with self.subTest(factor=factor):
                plan, args = self.plan("--photon-enhancement", str(factor))
                settings = {"channel": "fragmentation", "enhancement_factor": float(factor),
                            "detuning": 2.0}
                self.assertEqual(plan["options"]["photon_shower_sampling"], settings)
                self.assertEqual(plan["shard_jobs"], original["shard_jobs"])
                self.assertEqual(json.loads(json.dumps(plan)), plan)
                self.assertEqual(campaign._manifest_configuration(
                    self.measurement, args, plan)["photon_shower_sampling"], settings)
                for old_job, job in zip(original["jobs"], plan["jobs"]):
                    card = campaign._card_text(self.measurement, job)
                    if job["channel"] == "photon":
                        self.assertEqual(job, old_job)
                        self.assertEqual(card, campaign._card_text(self.measurement, old_job))
                    else:
                        self.assertEqual(job["photon_shower_sampling"], settings)
                        without_sampling = {k: v for k, v in job.items()
                                            if k != "photon_shower_sampling"}
                        self.assertEqual(without_sampling, old_job)
                        for line in (
                            f"set /Herwig/Shower/QtoQGammaSudakov:EnhancementFactor {factor}",
                            "set /Herwig/Shower/SplittingGenerator:Detuning 2",
                        ):
                            self.assertEqual(card.count(line), 1)
                            self.assertLess(card.index(line), card.index("saverun "))
                            card = card.replace(line + "\n", "")
                        self.assertEqual(card, campaign._card_text(self.measurement, old_job))

    def test_omission_preserves_legacy_configuration_and_cards(self):
        plan, args = self.plan()
        config = campaign._manifest_configuration(self.measurement, args, plan)
        self.assertNotIn("photon_shower_sampling", config)
        del args.photon_enhancement  # Namespace from the previous CLI version.
        legacy = campaign._plan(self.measurement, args)
        self.assertEqual(plan, legacy)
        self.assertEqual(config, campaign._manifest_configuration(self.measurement, args, legacy))
        for job in plan["jobs"]:
            self.assertNotIn("photon_shower_sampling", job)
            card = campaign._card_text(self.measurement, job)
            self.assertNotIn(":EnhancementFactor", card)
            self.assertNotIn(":Detuning", card)

    def test_control_is_explicit_and_settings_are_immutable_on_resume(self):
        plan, args = self.plan("--photon-enhancement", "1")
        configuration = campaign._manifest_configuration(self.measurement, args, plan)
        with tempfile.TemporaryDirectory() as directory, \
                patch.object(campaign, "CAMPAIGN_ROOT", Path(directory)), \
                patch.object(campaign, "_runtime", side_effect=AssertionError("No runtime")), \
                patch.object(experimental, "pending_jobs", return_value=([], [])) as pending:
            destination = Path(directory) / ID / args.tag
            destination.mkdir(parents=True)
            (destination / experimental.MANIFEST_NAME).write_text(json.dumps({
                "configuration": configuration, "jobs": plan["jobs"]}))
            with contextlib.redirect_stdout(io.StringIO()):
                campaign.run_pp(args, self.measurement)
            pending.assert_called_once()
            pending.reset_mock()
            for factor in (None, 2, 5, 10):
                changed = copy.copy(args)
                changed.photon_enhancement = factor
                with self.subTest(factor=factor), self.assertRaisesRegex(
                        campaign.CampaignError, "immutable manifest"):
                    campaign.run_pp(changed, self.measurement)
            pending.assert_not_called()

    def test_invalid_factors_and_other_measurements_fail_before_runtime(self):
        for factor in ("nan", "inf", "-inf", "0", "-2", "0.5", "100000001"):
            with self.subTest(factor=factor), self.assertRaises(campaign.CampaignError):
                self.plan(f"--photon-enhancement={factor}")
        for identifier in (ID + "_MEGammaJet-Only", "STAR_2019_I1708793", "COMPASS_2010_I843494"):
            with contextlib.redirect_stderr(io.StringIO()) as error, \
                    patch.object(experimental, "main") as legacy, \
                    patch.object(campaign, "_runtime") as runtime:
                code = campaign.main(["full", "--measurement", identifier, "--tag", "bad",
                                      "--photon-enhancement", "1", "--dry-run"])
            self.assertEqual(code, 2)
            self.assertIn("requires the expanded PHENIX", error.getvalue())
            legacy.assert_not_called()
            runtime.assert_not_called()

    def test_unsliced_smoke_and_comparison_families_receive_settings(self):
        args = campaign.make_parser().parse_args([
            "full", "--measurement", ID, "--tag", "smoke", "--smoke", "--shards", "1",
            "--families", "all", "--photon-enhancement", "5", "--dry-run"])
        plan = campaign._plan(self.measurement, args)
        for job in plan["jobs"]:
            self.assertEqual("photon_shower_sampling" in job, job["channel"] == "fragmentation")

    def test_pilot_dry_run_writes_nothing_and_has_requested_budgets(self):
        seeds = set()
        for factor, seed in ((1, 91401000), (2, 91411000), (5, 91421000), (10, 91431000)):
            with tempfile.TemporaryDirectory() as directory, \
                    patch.object(campaign, "CAMPAIGN_ROOT", campaign.CAMPAIGN_ROOT), \
                    patch.object(campaign, "_runtime", side_effect=AssertionError("No runtime")):
                output = io.StringIO()
                with contextlib.redirect_stdout(output):
                    code = wrapper.main([
                        "--output-root", directory, "full", "--tag", f"photon-C{factor}-D2",
                        "--jobs", "100", "--shards", "100", "--profile", "central",
                        "--families", "nominal", "--scales", "central", "--lo-events", "600000",
                        "--seed-base", str(seed), "--hard-pt-edges", "4,8,12,16,20,inf",
                        "--hard-pt-events", "fragmentation=6000000,6000000,6000000,6000000,6000000",
                        "--photon-enhancement", str(factor), "--dry-run",
                    ])
                plan = json.loads(output.getvalue())
                self.assertEqual(code, 0)
                self.assertEqual(plan["logical_jobs"], 40)
                self.assertEqual(plan["shard_jobs"], 4000)
                self.assertEqual(sum(j["events"] for j in plan["jobs"]), 132000000)
                for channel, per_shard in (("photon", 6000), ("fragmentation", 60000)):
                    self.assertEqual({j["events"] for j in plan["jobs"]
                                      if j["channel"] == channel}, {per_shard})
                current_seeds = {j["seed"] for j in plan["jobs"]}
                self.assertEqual(len(current_seeds), 4000)
                self.assertTrue(current_seeds.isdisjoint(seeds))
                seeds.update(current_seeds)
                self.assertEqual(list(Path(directory).iterdir()), [])


if __name__ == "__main__":
    unittest.main()
