#!/usr/bin/env python3
"""Regression coverage for the projected-spin STAR production controller."""

from __future__ import annotations

import contextlib
import importlib.util
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock


REPOSITORY = Path(__file__).resolve().parents[2]
CONTROL = (
    REPOSITORY
    / "campaigns"
    / "control"
    / "star510-pt13p1-bloch-20260822"
)


def load_controller():
    path = CONTROL / "controller.py"
    specification = importlib.util.spec_from_file_location(
        "star510_projection_controller", path
    )
    if specification is None or specification.loader is None:
        raise RuntimeError(f"Could not import {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


controller = load_controller()


class Star510ProjectionControllerTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = json.loads(
            (CONTROL / "campaigns.json").read_text(encoding="utf-8")
        )
        self.lock = json.loads(
            (CONTROL / "runtime-lock.json").read_text(encoding="utf-8")
        )

    def test_campaign_is_exactly_500m_per_helicity_and_2000_shards(self) -> None:
        production = self.config["star510"]
        self.assertEqual(len(production), 1)
        spec = production[0]
        self.assertEqual(spec["lo_events"], 500_000_000)
        self.assertEqual(spec["jet_kt_min_gev"], 4.0)
        self.assertEqual(spec["shards"], 500)
        self.assertEqual(spec["expected_shards"], 2_000)

    def test_completed_cut_scans_are_declared_read_only(self) -> None:
        scans = self.config["star_cut_scan"]
        self.assertEqual(
            [item["jet_kt_min_gev"] for item in scans], [3.0, 4.0, 5.0]
        )
        self.assertTrue(all(item["access"] == "read_only" for item in scans))
        self.assertEqual(sum(item["expected_shards"] for item in scans), 1_200)
        self.assertEqual(
            self.lock["star_cut_scan"]["campaign_source_commit"],
            "60da223dcc9ad2c78a5109e2bb1b1513fb131ef9",
        )

    def test_prepare_only_creates_the_new_nominal_star_campaign(self) -> None:
        verification = {"repository": {"commit": "a" * 40}}
        handoff = {
            "production_shards": 2_000,
            "read_only_cut_scan_shards": 1_200,
        }
        with (
            mock.patch.object(controller, "verify", return_value=verification),
            mock.patch.object(controller, "require_validation") as required,
            mock.patch.object(controller, "run") as run,
            mock.patch.object(controller, "assert_handoff", return_value=handoff),
            contextlib.redirect_stdout(io.StringIO()),
        ):
            controller.action_prepare()
        required.assert_called_once_with("a" * 40)
        self.assertEqual(run.call_count, 1)
        command = run.call_args.args[0]
        self.assertEqual(command[2], "prepare")
        self.assertIn("star510-pt13p1-bloch-500m-20260822-v1", command)
        self.assertNotIn("compat-star510-gencut-3gev-50m-20260819-v1", command)

    def test_lock_pins_source_overlay_guard_and_nonveto_policy(self) -> None:
        self.assertEqual(
            self.lock["herwig_source"]["commit"],
            "ad0c5486ad137e0e34c6274959c1fbf8323c5937",
        )
        self.assertIn("libThePEG.so", self.lock["thepeg_overlay"]["library"]["path"])
        self.assertIn("HerwigCore", self.lock["herwig"]["artifacts"])
        policy = self.lock["spin_density_policy"]
        self.assertEqual(policy["name"], "radial_bloch_ball_projection")
        self.assertTrue(policy["strict_negative_isr_guard"])
        self.assertFalse(policy["event_or_branching_veto"])
        self.assertEqual(set(self.lock["projection_evidence"]), {
            "three_seed_projection_audit", "counterfactual_replay"
        })

    def test_package_action_requests_only_star_high_resolution_plots(self) -> None:
        verification = {"repository": {"commit": "b" * 40}}
        with (
            mock.patch.object(controller, "verify", return_value=verification),
            mock.patch.object(controller, "require_validation"),
            mock.patch.object(controller, "run") as run,
        ):
            controller.action_package()
        command = run.call_args.args[0]
        self.assertEqual(command[-2:], ["--selection", "star510"])

    def test_smoke_yoda_validation_requires_finite_numeric_content(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "smoke.yoda"
            path.write_text("nan inf -inf 1.0 -2.5e-3\n", encoding="utf-8")
            self.assertEqual(controller._finite_yoda_count(path), 2)


if __name__ == "__main__":
    unittest.main()
