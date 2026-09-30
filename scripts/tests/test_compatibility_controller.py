#!/usr/bin/env python3
"""Regression coverage for the corrected-comparison production handoff."""

from __future__ import annotations

import copy
import contextlib
import hashlib
import importlib.util
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


REPOSITORY = Path(__file__).resolve().parents[2]
SCRIPTS = REPOSITORY / "scripts"
CONTROL = (
    REPOSITORY
    / "campaigns"
    / "control"
    / "compatibility-corrected-20260819"
)
sys.path.insert(0, str(SCRIPTS))

import check_star_generator_cut_stability as stability  # noqa: E402
import package_corrected_comparison_plots as packager  # noqa: E402
import run_phenomenology_campaign as campaign  # noqa: E402
import runtime_provenance as provenance  # noqa: E402


def load_controller():
    path = CONTROL / "controller.py"
    specification = importlib.util.spec_from_file_location(
        "compatibility_controller", path
    )
    if specification is None or specification.loader is None:
        raise RuntimeError(f"Could not import {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


controller = load_controller()


class RuntimeProvenanceTests(unittest.TestCase):
    def test_pdf_inventory_hashes_every_file_and_member(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            pdf = root / "PinnedSet"
            pdf.mkdir()
            (pdf / "PinnedSet.info").write_text("metadata\n", encoding="utf-8")
            (pdf / "PinnedSet_0000.dat").write_text(
                "member zero\n", encoding="utf-8"
            )
            inventory = provenance.pdf_set_inventory(root, "PinnedSet")
            self.assertEqual(inventory["file_count"], 2)
            self.assertEqual(inventory["member_file_count"], 1)
            self.assertEqual(
                [item["path"] for item in inventory["files"]],
                ["PinnedSet.info", "PinnedSet_0000.dat"],
            )

            expected = hashlib.sha256()
            for item in inventory["files"]:
                expected.update(
                    (
                        f"{item['path']}\0{item['size']}\0"
                        f"{item['sha256']}\n"
                    ).encode("utf-8")
                )
            self.assertEqual(
                inventory["inventory_sha256"], expected.hexdigest()
            )

    def test_prepared_inventory_digest_is_order_independent(self) -> None:
        first = {
            "b.run": {"size": 2, "sha256": "bb"},
            "a.in": {"size": 1, "sha256": "aa"},
        }
        second = {
            "a.in": {"sha256": "aa", "size": 1},
            "b.run": {"sha256": "bb", "size": 2},
        }
        self.assertEqual(
            provenance.inventory_digest(first),
            provenance.inventory_digest(second),
        )


class ControllerConfigurationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = json.loads(
            (CONTROL / "campaigns.json").read_text(encoding="utf-8")
        )

    def test_exact_production_and_cut_scan_inventory(self) -> None:
        fixed = self.config["fixed"]
        star = self.config["star510"]
        scans = self.config["star_cut_scan"]
        self.assertEqual(sum(item["expected_shards"] for item in fixed), 5600)
        self.assertEqual(sum(item["expected_shards"] for item in star), 2000)
        self.assertEqual(
            sum(item["expected_shards"] for item in fixed + star), 7600
        )
        self.assertEqual(sum(item["expected_shards"] for item in scans), 1200)
        self.assertEqual(
            [item["jet_kt_min_gev"] for item in scans], [3.0, 4.0, 5.0]
        )
        self.assertTrue(
            all(item["lo_events"] == 50_000_000 for item in scans)
        )
        self.assertEqual(star[0]["lo_events"], 500_000_000)
        self.assertEqual(star[0]["shards"], 500)
        self.assertTrue(
            all(
                item["posnlo_events"] == 3_000_000
                and item["negnlo_events"] == 300_000
                and item["shards"] == 100
                for item in fixed
            )
        )
        production_measurements = {
            item["measurement"] for item in fixed + star
        }
        self.assertNotIn("HERMES_2007_I726689_LEGACY", production_measurements)
        self.assertNotIn("PHENIX_2023_I2033856", production_measurements)

    def test_unchanged_campaigns_retain_their_exact_historical_pending_matrix(self) -> None:
        for group in ("fixed", "star_cut_scan", "star510"):
            for item in self.config[group]:
                if item["measurement"] == "HERMES_2007_I726689":
                    # This frozen controller records the earlier proton-only
                    # definition. Its dual-target incompatibility is checked
                    # separately without changing archived event budgets.
                    continue
                argv = ["prepare", *controller.generation_arguments(item), "--dry-run"]
                output = io.StringIO()
                with contextlib.redirect_stdout(output):
                    self.assertEqual(campaign.main(argv), 0)
                plan = json.loads(output.getvalue())
                self.assertEqual(plan["shard_jobs"], item["expected_shards"])
                self.assertEqual(len(plan["jobs"]), item["expected_shards"])
                self.assertTrue(
                    all(job["status"] == "planned" for job in plan["jobs"])
                )
                self.assertEqual(
                    {job["family"] for job in plan["jobs"]}, {"nominal"}
                )

    def test_historical_hermes_guard_rejects_the_current_dual_target_matrix(self) -> None:
        item = next(
            spec for spec in self.config["fixed"]
            if spec["measurement"] == "HERMES_2007_I726689"
        )
        self.assertEqual(item["tag"], "compat-central-3m-20260819-v1")
        self.assertEqual(item["expected_shards"], 800)
        argv = ["prepare", *controller.generation_arguments(item), "--dry-run"]
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            self.assertEqual(campaign.main(argv), 0)
        plan = json.loads(output.getvalue())
        self.assertEqual(plan["shard_jobs"], 1600)
        self.assertEqual(len(plan["jobs"]), 1600)
        self.assertEqual({job["component"] for job in plan["jobs"]}, {"P", "N"})
        manifest = {"status": "prepared", "jobs": plan["jobs"]}
        with (
            mock.patch.object(controller, "_manifest", return_value=manifest),
            mock.patch.object(controller, "verify_manifest_provenance") as verify,
        ):
            with self.assertRaisesRegex(controller.ControllerError, "has 1600 shards, expected 800"):
                controller.assert_prepared(item, "a" * 40)
        verify.assert_not_called()

    def test_prepare_action_has_an_enforced_execution_hard_stop(self) -> None:
        verification = {"repository": {"commit": "a" * 40}}
        handoff = {
            "production_shards": 7600,
            "cut_scan_shards": 1200,
        }
        with (
            mock.patch.object(controller, "verify", return_value=verification),
            mock.patch.object(controller, "require_validation") as required,
            mock.patch.object(controller, "run") as run,
            mock.patch.object(
                controller, "assert_handoff", return_value=handoff
            ),
        ):
            with contextlib.redirect_stdout(io.StringIO()):
                controller.action_prepare()
        required.assert_called_once_with("a" * 40)
        commands = [call.args[0] for call in run.call_args_list]
        self.assertEqual(len(commands), 9)
        self.assertTrue(all(command[2] == "prepare" for command in commands))
        self.assertFalse(
            any("campaign" in command[2:] for command in commands)
        )

    def test_runtime_lock_pins_all_requested_binary_inputs(self) -> None:
        lock = json.loads(
            (CONTROL / "runtime-lock.json").read_text(encoding="utf-8")
        )
        self.assertEqual(
            set(lock["herwig"]["artifacts"]),
            {
                "Herwig",
                "HerwigDefaults.rpo",
                "HwMEDIS",
                "HwMEHadron",
                "HwShower",
                "FixedTargetLuminosity",
            },
        )
        self.assertEqual(
            {
                name: (item["file_count"], item["member_file_count"])
                for name, item in lock["pdf_sets"].items()
            },
            {
                "NNPDF40_nlo_pch_as_01180": (102, 101),
                "NNPDFpol20_nlo_as_01180": (102, 101),
            },
        )


def fake_prediction() -> dict[str, dict[str, list[float]]]:
    result: dict[str, dict[str, list[float]]] = {}
    for observable in stability.PRIMARY_OBSERVABLES:
        result[f"SigmaUU_{observable}"] = {
            "edges": [float(index) for index in range(7)],
            "values": [100.0, 80.0, 60.0, 50.0, 40.0, 30.0],
            "errors": [0.0] * 6,
        }
        result[observable] = {
            "edges": [float(index) for index in range(7)],
            "values": [0.01, 0.02, 0.03, 0.04, 0.05, 0.06],
            "errors": [0.0] * 6,
        }
    return result


def fake_policy_inputs() -> tuple[dict[str, object], dict[str, object]]:
    measurement = {
        "comparison_policy": {
            "primary_bin_masks": {
                "inclusive": {
                    "first_bin": 5,
                    "minimum_analysis_low_gev": 13.1,
                    "excluded_treatment": "diagnostic_only",
                }
            },
            "covariance_points": 59,
        }
    }
    snapshot = {
        "datasets": {
            observable: {
                "bin_edges": (
                    [7.0, 8.2, 9.6, 11.2, 13.1, 15.3, 17.9]
                    if observable == "inclusive"
                    else [float(index) for index in range(7)]
                )
            }
            for observable in stability.PRIMARY_OBSERVABLES
        }
    }
    return measurement, snapshot


def fake_campaign(cut: float, seed: int) -> dict[str, object]:
    return {
        "directory": f"/campaign/{cut:g}",
        "cut_gev": cut,
        "prediction": fake_prediction(),
        "source_commit": "b" * 40,
        "manifest_sha256": f"manifest-{cut:g}",
        "summary_sha256": f"summary-{cut:g}",
        "initial_seeds": {seed},
    }


class StarGeneratorCutStabilityTests(unittest.TestCase):
    def test_first_two_finite_bins_form_the_only_gate(self) -> None:
        reference = fake_campaign(3.0, 1)
        alternate = fake_campaign(4.0, 2)
        measurement, snapshot = fake_policy_inputs()
        report = stability.compare_campaigns(
            reference, alternate, gate=True,
            measurement=measurement, snapshot=snapshot,
        )
        self.assertTrue(report["passed"])
        gated = [row for row in report["rows"] if row["gated"]]
        self.assertEqual(len(gated), 20)
        inclusive = [row for row in gated if row["observable"] == "inclusive"]
        dijets = [row for row in gated if row["observable"] != "inclusive"]
        self.assertEqual({row["bin"] for row in inclusive}, {5, 6})
        self.assertEqual({row["bin"] for row in dijets}, {1, 2})

        failed = copy.deepcopy(alternate)
        failed["prediction"]["SigmaUU_inclusive"]["values"][0] = 200.0
        report = stability.compare_campaigns(
            reference, failed, gate=True,
            measurement=measurement, snapshot=snapshot,
        )
        self.assertTrue(report["passed"])
        failed["prediction"]["SigmaUU_inclusive"]["values"][4] = 80.0
        report = stability.compare_campaigns(
            reference, failed, gate=True,
            measurement=measurement, snapshot=snapshot,
        )
        self.assertFalse(report["passed"])
        self.assertIn("inclusive:sigma_uu:bin5", report["failures"])

    def test_report_requires_one_commit_and_disjoint_seeds(self) -> None:
        campaigns = {
            3.0: fake_campaign(3.0, 1),
            4.0: fake_campaign(4.0, 2),
            5.0: fake_campaign(5.0, 3),
        }
        measurement, snapshot = fake_policy_inputs()
        report = stability.build_report(campaigns, measurement, snapshot)
        self.assertTrue(report["gate_passed"])
        self.assertEqual(report["schema_version"], 2)
        self.assertEqual(report["campaign_source_commit"], "b" * 40)
        self.assertFalse(report["four_vs_five_gev_stress"]["gate"])
        campaigns[5.0]["initial_seeds"] = {2}
        with self.assertRaisesRegex(stability.StabilityError, "reuse"):
            stability.build_report(campaigns, measurement, snapshot)


class CorrectedPlotPackageTests(unittest.TestCase):
    def test_curated_selection_is_exactly_four_plus_twenty_four_plus_five(self) -> None:
        config = json.loads(
            (CONTROL / "campaigns.json").read_text(encoding="utf-8")
        )
        fixed = packager.fixed_plot_specs(config)
        sidis = packager.sidis_plot_specs(config)
        star = packager.star_plot_specs(config)
        self.assertEqual((len(fixed), len(sidis), len(star)), (4, 24, 5))
        destinations = [
            str(item["destination"]) for item in fixed + sidis + star
        ]
        self.assertEqual(len(destinations), len(set(destinations)))

    def test_historical_hermes_package_selects_the_proton_dataset_by_id(self) -> None:
        config = {"fixed": [{"measurement": "HERMES_2007_I726689"}]}
        snapshot = {
            "schema_version": 2,
            "datasets": [
                {"id": "D", "rivet_path": "/REF/HERMES_2007_I726689/d14-x01-y02"},
                {"id": "P", "rivet_path": "/REF/HERMES_2007_I726689/d14-x01-y01"},
            ],
        }
        with mock.patch.object(packager, "load_json", return_value=snapshot):
            entries = packager.fixed_plot_specs(config)
        self.assertEqual(len(entries), 1)
        self.assertEqual(entries[0]["stem"], "d14-x01-y01")
        self.assertEqual(str(entries[0]["destination"]), "dis/hermes-a1")
        for datasets in ([snapshot["datasets"][0]], [snapshot["datasets"][1]] * 2):
            with self.subTest(datasets=datasets):
                with mock.patch.object(packager, "load_json", return_value={"datasets": datasets}):
                    with self.assertRaisesRegex(packager.PackageError, "exactly one proton P"):
                        packager.fixed_plot_specs(config)

    def test_star_only_high_resolution_selection_is_explicit(self) -> None:
        arguments = packager.make_parser().parse_args(
            [
                "--campaign-config",
                "campaigns.json",
                "--output-root",
                "figures",
                "--source-commit",
                "a" * 40,
                "--selection",
                "star510",
            ]
        )
        self.assertEqual(arguments.selection, "star510")


if __name__ == "__main__":
    unittest.main()
