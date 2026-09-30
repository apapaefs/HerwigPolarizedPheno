#!/usr/bin/env python3
from __future__ import annotations

import copy
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import hermes_born_reference_data as reference


class HermesBornReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.snapshot = reference.validate_vendored()

    def test_full_published_born_columns_and_target_identities(self) -> None:
        self.assertEqual(len(self.snapshot["datasets"]), 38)
        self.assertEqual(len({d["id"] for d in self.snapshot["datasets"]}), 38)
        for target, table in (("P", 7), ("D", 8)):
            panels = [d for d in self.snapshot["datasets"] if d["target"] == target]
            points = [p for d in panels for p in d["points"]]
            raw = json.loads((ROOT / reference.DATA / f"raw-born-apar-table-{table}.json").read_text())
            self.assertEqual(len(points), 45)
            self.assertEqual([p["global_bin"] for p in points], list(range(1, 46)))
            for point, row in zip(points, raw["values"]):
                measured, born = row["y"]
                self.assertEqual(point["value"], float(born["value"]))
                self.assertEqual(point["stat"], float(born["errors"][0]["symerror"]))
                self.assertEqual(point["systematic_combined"], float(born["errors"][1]["symerror"]))
                self.assertEqual(point["y_mean"], float(row["x"][1]["value"]))
                self.assertEqual(point["q2_mean"], float(row["x"][2]["value"]))
            self.assertNotEqual([p["value"] for p in points],
                                [float(row["y"][0]["value"]) for row in raw["values"]])
            self.assertEqual([p["value"] for p in (points[0], points[-1])],
                             [0.0205, 0.5274] if target == "P" else [0.0086, 0.5164])
            self.assertTrue(all(d["rivet_path"].startswith(f"/REF/HERMES_2007_I726689/d{table:02d}-")
                                for d in panels))

    def test_low_q2_support_comes_from_whole_cell_bounds(self) -> None:
        for target in ("P", "D"):
            panels = [d for d in self.snapshot["datasets"] if d["target"] == target]
            self.assertEqual([len(d["points"]) for d in panels], [1]*4 + [2]*4 + [3]*11)
            supported = [p["global_bin"] for d in panels for p in d["points"] if p["supported_prediction"]]
            self.assertEqual(supported, [6, 8, 10, 12] + list(range(13, 46)))
            self.assertEqual(len(supported), 37)
            for d in panels:
                for p in d["points"]:
                    self.assertEqual(p["supported_prediction"], p["q2_low"] >= 1)
                    self.assertLess(p["q2_low"], p["q2_mean"])
                    self.assertGreater(p["q2_high"], p["q2_mean"])

    def test_numerical_thesis_internal_edges_match_published_eps_geometry(self) -> None:
        # This checks independently drawn physical cuts, not midpoints of data
        # means.  Axis ticks at Q2=1,10 are y=842,1375 in the EPS integer frame.
        eps = gzip.decompress((ROOT / reference.DATA / "raw-born-apar-fig4.eps.gz").read_bytes()).decode()
        lines = [(int(x), int(y)) for x, y, width in
                 re.findall(r"(\d+)\s+(\d+)\s+m\s+(\d+)\s+X", eps)
                 if int(x) in {1340,1413,1485,1558,1630,1703,1775,1848,1920,1993,2065}]
        xs = [1340,1413,1485,1558,1630,1703,1775,1848,1920,1993,2065]
        for x, edges in zip(xs, reference.Q2_EDGES[8:]):
            observed = [y for xx, y in lines if xx == x and y > 850]
            self.assertEqual(len(observed), 2)
            for value, y in zip(edges[1:-1], sorted(observed)):
                transformed = 842.5 + (1600/3)*math.log10(value)
                self.assertLess(abs(y-transformed), 1.0)
        self.assertEqual(reference.Q2_EDGES[8], [1, 1.505, 2.265, 20])
        self.assertEqual(reference.Q2_EDGES[-1], [1, 7.655, 12.528, 20])
        # Published rounded x limits are retained separately from physical cuts.
        self.assertEqual(self.snapshot["datasets"][8]["x_low"], 0.05681)
        self.assertEqual(self.snapshot["datasets"][8]["display_x_low"], 0.0568)

    def test_missing_hepdata_means_are_paper_values_not_bin_midpoints(self) -> None:
        for target in ("P", "D"):
            points = [p for d in self.snapshot["datasets"] if d["target"] == target for p in d["points"]]
            for index, value in {4: 0.0190, 10: 0.0403, 12: 0.0506}.items():
                self.assertEqual(points[index-1]["x_mean"], value)
                self.assertEqual(points[index-1]["x_mean_source"], "paper Table XI/XII")

    def test_full_statistical_covariance_uses_same_global_born_bin_order(self) -> None:
        for target, table in (("P", 18), ("D", 19)):
            points = [p for d in self.snapshot["datasets"] if d["target"] == target for p in d["points"]]
            covariance = self.snapshot["statistical_covariance"][target]
            matrix, correlation = covariance["matrix"], covariance["correlation"]
            self.assertEqual(covariance["global_bin_order"], list(range(1, 46)))
            raw = json.loads((ROOT / reference.DATA / f"raw-born-apar-table-{table}.json").read_text())
            for row in raw["values"]:
                i,j = [int(float(x["value"]))-1 for x in row["x"]]
                self.assertEqual(correlation[i][j], float(row["y"][0]["value"]))
                self.assertAlmostEqual(matrix[i][j], points[i]["stat"]*points[j]["stat"]*correlation[i][j], places=17)
            for i in range(45):
                self.assertAlmostEqual(matrix[i][i], points[i]["stat"]**2, places=17)
                for j in range(45):
                    self.assertEqual(matrix[i][j], matrix[j][i])
            # The supplied anti-correlation must affect an actual two-cell
            # linear-combination variance; a diagonal approximation is different.
            correlated = matrix[0][0] + matrix[1][1] + 2*matrix[0][1]
            self.assertLess(correlated, matrix[0][0] + matrix[1][1])

    def test_rounded_published_correlations_remain_positive_definite(self) -> None:
        for target in ("P", "D"):
            matrix = self.snapshot["statistical_covariance"][target]["correlation"]
            lower = [[0.0]*45 for _ in range(45)]
            for i in range(45):
                for j in range(i+1):
                    residual = matrix[i][j] - sum(lower[i][k]*lower[j][k] for k in range(j))
                    if i == j:
                        self.assertGreater(residual, 0)
                        lower[i][j] = math.sqrt(residual)
                    else:
                        lower[i][j] = residual/lower[j][j]

    def test_no_reference_corrections_or_fabricated_integrated_data(self) -> None:
        for dataset in self.snapshot["datasets"]:
            provenance = dataset["provenance"]
            self.assertIn("retained unchanged", provenance["correction_convention"])
            self.assertEqual(provenance["normalization_uncertainty"], 0.052 if dataset["target"] == "P" else 0.05)
            self.assertEqual(provenance["normalization_uncertainty_included_in"], "published systematic uncertainty")
            self.assertNotIn("depolarization", dataset)
            self.assertNotIn("d_state_factor", dataset)
        self.assertIn("unpolarized yield/cross-section weights", self.snapshot["projection_reference_policy"])
        self.assertEqual(set(self.snapshot["statistical_covariance"]), {"P", "D"})

    def test_all_vendored_source_files_have_pinned_checksums(self) -> None:
        self.assertEqual(len(self.snapshot["source_files"]), 5)
        for source in self.snapshot["source_files"]:
            self.assertEqual(hashlib.sha256((ROOT / source["path"]).read_bytes()).hexdigest(), source["sha256"])

    def test_source_and_normalized_mutation_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for source in self.snapshot["source_files"]:
                dest = root/source["path"]
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(ROOT/source["path"], dest)
            normalized = root/reference.REFERENCE_PATH
            normalized.write_text(json.dumps(self.snapshot))
            self.assertEqual(reference.validate_vendored(root), self.snapshot)
            altered = copy.deepcopy(self.snapshot)
            altered["datasets"][0]["points"][0]["value"] *= 0.925
            normalized.write_text(json.dumps(altered))
            with self.assertRaisesRegex(reference.BornReferenceError, "differs"):
                reference.validate_vendored(root)
            raw = root/self.snapshot["source_files"][0]["path"]
            raw.write_bytes(raw.read_bytes()+b"\n")
            with self.assertRaisesRegex(reference.BornReferenceError, "Checksum"):
                reference.build_snapshot(root)


if __name__ == "__main__":
    unittest.main()
