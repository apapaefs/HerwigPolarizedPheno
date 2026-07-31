#!/usr/bin/env python3
"""Regression tests for Herwig-local experimental campaign filenames."""

from __future__ import annotations

import os
import sys
import unittest
from pathlib import Path


DISPOL_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(DISPOL_ROOT / "scripts"))

import run_experimental_campaign as campaign  # noqa: E402


class ExperimentalOutputFilenameTests(unittest.TestCase):
    def setUp(self) -> None:
        self.job = {
            "stem": (
                "HERMES_2019_I1698889_sidis-targetP-levelstable_particle-"
                "nominal-PP-POSNLO-p000-u000-mu0p5-mpioff"
            ),
            "seed": 1699529,
            "shard": 1,
            "shards": 10,
            "attempt": 2,
        }
        self.measurement = {"id": "HERMES_2019_I1698889"}
        self.campaign_tag = "hermes_sidis_scales_300k_20260723"

    def test_attempt_tag_is_compact_and_unique(self) -> None:
        tag = campaign._job_attempt_tag(
            self.job, self.measurement, self.campaign_tag
        )
        self.assertEqual(tag, "s001-of-010-a02")

        next_shard = dict(self.job, shard=2)
        next_attempt = dict(self.job, attempt=3)
        self.assertNotEqual(
            tag,
            campaign._job_attempt_tag(
                next_shard, self.measurement, self.campaign_tag
            ),
        )
        self.assertNotEqual(
            tag,
            campaign._job_attempt_tag(
                next_attempt, self.measurement, self.campaign_tag
            ),
        )

    def test_failing_sidis_case_stays_below_name_max(self) -> None:
        tag = campaign._job_attempt_tag(
            self.job, self.measurement, self.campaign_tag
        )
        output_stem = f"{self.job['stem']}-S{self.job['seed']}-{tag}"
        for suffix in (".yoda", "-EvtGen.log"):
            filename = os.fsencode(output_stem + suffix)
            self.assertLessEqual(len(filename), 255)


if __name__ == "__main__":
    unittest.main()
