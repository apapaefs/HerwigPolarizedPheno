#!/usr/bin/env python3
"""Changing worker concurrency preserves the prepared DIS campaign."""
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


class RuntimeConcurrencyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.measurement = campaign.get_measurement("HERMES_2007_I726689")
        cls.signature = campaign.measurement_signature(cls.measurement)

    def arguments(self, jobs=380, dry_run=False):
        command = [
            "campaign", "--measurement", self.measurement["id"],
            "--tag", "runtime-concurrency", "--smoke", "--shards", "1",
            "--seed-base", "5000", "--jobs", str(jobs),
            "--progress-interval", "-1",
        ]
        if dry_run:
            command.append("--dry-run")
        return campaign.make_parser().parse_args(command)

    def manifest(self, root, args):
        options = campaign._resolved_campaign_options(args, self.measurement)
        options["jobs"] = 100
        jobs = campaign.build_job_matrix(self.measurement, 100, 100, 1, 5000)
        # Existing outputs must remain complete; only the other saved jobs run.
        output = root / jobs[0]["output_yoda"]
        output.parent.mkdir(parents=True)
        output.write_text("existing completed shard", encoding="utf-8")
        jobs[0]["status"] = "success"
        jobs[0]["output_size"] = output.stat().st_size
        result = {
            "configuration": campaign._manifest_configuration(
                self.measurement, args.tag, options, self.signature
            ),
            "measurement": self.measurement["id"], "tag": args.tag,
            "runtime": {"toy": True}, "status": "prepared", "jobs": jobs,
            "history": [{"action": "prepared"}],
        }
        campaign.atomic_write_json(root / campaign.MANIFEST_NAME, result)
        return result

    def test_changed_or_missing_historical_jobs_are_compatible_without_mutation(self):
        options = campaign._resolved_campaign_options(self.arguments(), self.measurement)
        expected = campaign._manifest_configuration(
            self.measurement, "runtime-concurrency", options, self.signature
        )
        for recorded, requested in ((100, 380), (None, 380), (100, None)):
            with self.subTest(recorded=recorded, requested=requested):
                actual = copy.deepcopy(expected)
                desired = copy.deepcopy(expected)
                if recorded is None:
                    actual.pop("jobs")
                else:
                    actual["jobs"] = recorded
                if requested is None:
                    desired.pop("jobs")
                else:
                    desired["jobs"] = requested
                manifest = {"configuration": actual, "jobs": [{"seed": 5000}]}
                original = copy.deepcopy((manifest, desired))
                campaign._assert_manifest_compatible(manifest, desired)
                self.assertEqual((manifest, desired), original)

    def test_every_other_configuration_key_remains_locked(self):
        options = campaign._resolved_campaign_options(self.arguments(), self.measurement)
        options.update(
            comparisons=True, profile="paper", variation_points=[[0, 0, 1.0], [1, 0, 0.5]]
        )
        expected = campaign._manifest_configuration(
            self.measurement, "runtime-concurrency", options, self.signature
        )
        expected["families"] = ["nominal", "lo"]
        expected["profile_options"] = {"jobs": 100}
        actual = copy.deepcopy(expected)
        actual["jobs"] = 100
        for key, value in expected.items():
            if key == "jobs":
                continue
            with self.subTest(key=key):
                changed = copy.deepcopy(expected)
                if isinstance(value, bool):
                    changed[key] = not value
                elif isinstance(value, int):
                    changed[key] = value + 1
                elif isinstance(value, list):
                    changed[key] = value + ["changed"]
                elif isinstance(value, dict):
                    changed[key]["jobs"] = 380
                else:
                    changed[key] = value + "-changed"
                with self.assertRaisesRegex(campaign.CampaignError, key):
                    campaign._assert_manifest_compatible({"configuration": actual}, changed)
                missing = copy.deepcopy(expected)
                missing.pop(key)
                with self.assertRaisesRegex(campaign.CampaignError, key):
                    campaign._assert_manifest_compatible({"configuration": actual}, missing)
        with self.assertRaises(campaign.CampaignError):
            campaign._assert_manifest_compatible({}, expected)

    def test_existing_manifest_dry_run_schedules_saved_jobs_without_writes_or_generation(self):
        args = self.arguments(dry_run=True)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            original = self.manifest(root, args)
            path = root / campaign.MANIFEST_NAME
            original_bytes = path.read_bytes()
            original_files = {item.relative_to(root) for item in root.rglob("*")}
            output = io.StringIO()
            with mock.patch.object(campaign, "campaign_directory", return_value=root), \
                 mock.patch.object(campaign, "prepare_campaign") as prepare, \
                 mock.patch.object(campaign, "build_job_matrix") as build, \
                 mock.patch.object(campaign, "atomic_write_json") as write, \
                 mock.patch.object(campaign, "_run_tracked_campaign_job") as worker, \
                 contextlib.redirect_stdout(output):
                self.assertEqual(campaign.run_campaign(args, self.measurement), root)
            scheduled = json.loads(output.getvalue())
            expected_pending = copy.deepcopy(original["jobs"][1:])
            for job in expected_pending:
                job["status"] = "queued"
            self.assertEqual(scheduled["scheduled"], expected_pending)
            self.assertEqual(scheduled["already_complete"], 1)
            self.assertEqual(path.read_bytes(), original_bytes)
            self.assertEqual({item.relative_to(root) for item in root.rglob("*")}, original_files)
            for function in (prepare, build, write, worker):
                function.assert_not_called()

    def test_runtime_pool_and_history_use_actual_workers_preserving_prepared_configuration(self):
        executor_class = campaign.concurrent.futures.ThreadPoolExecutor
        for requested, effective in ((380, 380), (-4, 1)):
            with self.subTest(requested=requested), tempfile.TemporaryDirectory() as temporary:
                args = self.arguments(jobs=requested)
                root = Path(temporary)
                original = self.manifest(root, args)

                def fake_job(job, measurement, campaign_dir, runtime, campaign_tag, active, lock):
                    output = campaign_dir / job["output_yoda"]
                    output.parent.mkdir(parents=True, exist_ok=True)
                    output.write_text("toy completed shard", encoding="utf-8")
                    return campaign.JobResult(job["id"], True, 0, output.stat().st_size, "ok")

                with mock.patch.object(campaign, "campaign_directory", return_value=root), \
                     mock.patch.object(campaign, "prepare_campaign") as prepare, \
                     mock.patch.object(campaign, "build_job_matrix") as build, \
                     mock.patch.object(campaign, "_run_tracked_campaign_job", side_effect=fake_job) as worker, \
                     mock.patch.object(campaign.concurrent.futures, "ThreadPoolExecutor", wraps=executor_class) as pool, \
                     contextlib.redirect_stdout(io.StringIO()):
                    campaign.run_campaign(args, self.measurement)
                pool.assert_called_once_with(max_workers=effective)
                prepare.assert_not_called()
                build.assert_not_called()
                final = json.loads((root / campaign.MANIFEST_NAME).read_text())
                self.assertEqual(final["configuration"], original["configuration"])
                self.assertEqual(final["configuration"]["jobs"], 100)
                self.assertEqual(final["jobs"][0], original["jobs"][0])
                for before, after in zip(original["jobs"], final["jobs"]):
                    for key, value in before.items():
                        if key != "status":
                            self.assertEqual(after[key], value)
                self.assertEqual(
                    {call.args[0]["id"] for call in worker.call_args_list},
                    {job["id"] for job in original["jobs"][1:]},
                )
                start = next(entry for entry in final["history"] if entry["action"] == "campaign")
                self.assertEqual(start["max_workers"], effective)
                self.assertEqual(start["scheduled"], len(original["jobs"]) - 1)
                self.assertEqual(final["status"], "complete")


if __name__ == "__main__":
    unittest.main()
