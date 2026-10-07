#!/usr/bin/env python3
"""Actual-YODA support validation with unequal shard statistics and order signs."""
from __future__ import annotations

import json
import math
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"scripts"))
import check_hermes_generation_support as support


def yoda_block(name, edges, raw_rows, scale):
    path = f"/{support.ANALYSIS}/{name}"
    rows = [[0.]*5, *raw_rows, [0.]*5]
    raw = [f"BEGIN YODA_HISTO1D_V3 /RAW{path}", f"Edges(A1): {json.dumps(edges)}", "---"]
    raw += [" ".join(map(str, row)) for row in rows]
    raw += ["END YODA_HISTO1D_V3"]
    estimate = [f"BEGIN YODA_ESTIMATE1D_V3 {path}", f"ScaledBy: {scale:.17g}",
                f"Edges(A1): {json.dumps(edges)}", "---", "nan --- ---"]
    for i, row in enumerate(raw_rows):
        width = edges[i+1]-edges[i]
        value, error = scale*row[0]/width, abs(scale)*math.sqrt(row[1])/width
        estimate.append(f"{value:.16g} {-error:.16g} {error:.16g}")
    estimate += ["nan --- ---", "END YODA_ESTIMATE1D_V3"]
    return "\n".join(raw+estimate)+"\n"


def make_campaign(directory, missing_low_w2=False):
    jobs = []
    expected = {selection: support.Moments() for selection in support.SELECTIONS}
    for target in ("P", "N"):
        for hi, helicity in enumerate(support.HELICITIES):
            for order, sign in support.ORDERS.items():
                for shard, events in ((1, 1000), (2, 3000)):
                    name = f"{target}-{helicity}-{order}-{shard}"
                    scale = (1+hi*.1)*(2. if shard == 1 else .5)
                    count = (100 if shard == 1 else 300) if order == "POSNLO" else 10
                    q = (6. if shard == 1 else 6.4) if order == "POSNLO" else 6.7
                    text = ""
                    for selection, prefix in support.SELECTIONS.items():
                        n = count if selection == "ring_control" else count*.7
                        row = [n, n, n*q, n*(q*q+.04), n]
                        text += yoda_block(prefix+"BornSigmaQ2_X19", [1., 7.655, 12.528, 20.],
                                           [row, [0.]*5, [0.]*5], scale)
                        low_n = 0. if missing_low_w2 else n
                        low_rows = [[low_n, low_n, low_n*4., low_n*16.01, low_n] for _ in range(28)]
                        text += yoda_block(prefix+"Accepted_W2Fine_Q2GT1",
                                           [3.24+.1*i for i in range(29)], low_rows, scale)
                        if target == "P":
                            coefficient = events/4000*sign/4*scale
                            expected[selection].add(support.Moments(n, n*q, n, n*q,
                                                                    n*(q*q+.04), n), coefficient)
                    (directory/f"{name}.yoda").write_text(text)
                    jobs.append({"id": name, "component": target, "helicity": helicity,
                                 "order": order, "shard": shard, "events": events,
                                 "status": "success", "output_yoda": f"{name}.yoda"})
    manifest = {"configuration": {"measurement": support.ANALYSIS, "shards": 2}, "jobs": jobs}
    (directory/"manifest.json").write_text(json.dumps(manifest))
    return expected


class ActualCampaignSupportTests(unittest.TestCase):
    def test_weighted_shards_signed_orders_and_four_helicities(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            expected = make_campaign(directory)
            report = support.inspect_campaign(directory)
            self.assertTrue(report["passed"], report["failures"])
            self.assertEqual(len(report["inputs"]), 32)
            for selection in support.SELECTIONS:
                actual = report["targets"]["P"][selection]["cell43"]
                correct = expected[selection].report()
                for key in correct:
                    self.assertAlmostEqual(actual[key], correct[key], places=10, msg=key)
                # A raw accumulator merge or an equal-shard average would fail.
                self.assertGreater(abs(actual["mean_q2"] - 6.2), .05)

    def test_missing_low_mass_support_and_insufficient_statistics_fail(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            make_campaign(directory, missing_low_w2=True)
            report = support.inspect_campaign(directory, minimum_effective_entries=1e9)
            self.assertFalse(report["passed"])
            self.assertTrue(any("low-W2 bins" in failure for failure in report["failures"]))
            self.assertTrue(any("insufficient" in failure for failure in report["failures"]))

    def test_incomplete_manifest_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            make_campaign(directory)
            path = directory/"manifest.json"
            manifest = json.loads(path.read_text())
            manifest["jobs"][0]["status"] = "running"
            path.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, "Unfinished"):
                support.inspect_campaign(directory)

    def test_nonconstant_raw_weights_do_not_receive_invented_mean_errors(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)/"fixture.yoda"
            # Two weights 1 and 2: Cauchy equality is false (9 != 2*5).
            path.write_text(yoda_block("BornSigmaQ2_X19", [1., 7.655],
                                       [[3., 5., 18., 108., 2.]], 2.))
            objects = support.read_yoda_objects(path)
            with self.assertRaisesRegex(ValueError, "Nonconstant"):
                support.normalized_bin(objects, "BornSigmaQ2_X19", 0)

    def test_finalized_normalization_must_match_raw_moments(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)/"fixture.yoda"
            path.write_text(yoda_block("BornSigmaQ2_X19", [1., 7.655],
                                       [[3., 3., 18., 108., 3.]], 2.))
            objects = support.read_yoda_objects(path)
            objects[f"/{support.ANALYSIS}/BornSigmaQ2_X19"]["scale"] = 3.
            with self.assertRaisesRegex(ValueError, "normalization"):
                support.normalized_bin(objects, "BornSigmaQ2_X19", 0)


if __name__ == "__main__":
    unittest.main()
