"""Published-number fidelity and the open-charm production boundary."""
from __future__ import annotations

import contextlib
import io
import json
import math
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import compass_open_charm_reference as charm
import build_results_browser as browser
import run_phenomenology_campaign as campaign


class OpenCharmReferenceTests(unittest.TestCase):
    def test_published_cells_uncertainties_and_means(self):
        snapshot = charm.validate(ROOT)
        groups = snapshot["datasets"]
        self.assertEqual(len(groups), 9)
        self.assertEqual(sum(len(d["points"]) for d in groups), 45)
        # Independent anchors from Tables 5, 6 and 7, including estimates outside ±1.
        self.assertEqual(groups[0]["points"][0]["value"], -0.90)
        self.assertEqual(groups[2]["points"][2]["value"], 1.49)
        self.assertEqual(groups[5]["points"][0]["value"], -2.55)
        extreme = groups[6]["points"][0]
        self.assertEqual((extreme["value"], extreme["stat"], extreme["syst"]), (7.03, 4.74, 0.71))
        self.assertEqual(extreme["means"], {"y": .46, "Q2_GeV2": .38, "pt_GeV": .22, "E_GeV": 27.7, "D": .58})
        for dataset in groups:
            self.assertIsNone(dataset["points"][-1]["pt_bin_GeV"][1])
            for point in dataset["points"]:
                self.assertAlmostEqual(point["total"], math.hypot(point["stat"], point["syst"]))
        self.assertIsNone(groups[-1]["energy_bin_GeV"][1])

    @unittest.skipUnless(shutil.which("pdftotext"), "PDF verification requires Poppler pdftotext")
    def test_every_cell_and_mean_matches_independently_extracted_pdf(self):
        charm.validate(ROOT, verify_pdf_text=True)

    def temporary_repo(self):
        temp = tempfile.TemporaryDirectory(prefix="open-charm-reference-test-")
        self.addCleanup(temp.cleanup)
        repo = Path(temp.name)
        shutil.copytree(ROOT / charm.DATA, repo / charm.DATA)
        (repo / "config/reference").mkdir(parents=True)
        shutil.copyfile(ROOT / "config/reference" / f"{charm.ID}.json",
                        repo / "config/reference" / f"{charm.ID}.json")
        return repo

    def test_source_and_normalized_changes_fail_closed(self):
        for name in ["tables-5-7.csv", "arXiv-1211.6849v2.pdf", "reference.json", "source-manifest.json"]:
            with self.subTest(name=name):
                repo = self.temporary_repo()
                path = repo / charm.DATA / name
                path.write_bytes(path.read_bytes() + b"\n")
                with self.assertRaises(ValueError):
                    charm.validate(repo)

    def test_duplicate_missing_nonfinite_and_invented_upper_edges_rejected(self):
        raw = (ROOT / charm.DATA / "tables-5-7.csv").read_text()
        for changed in [raw + raw.splitlines()[1] + "\n", raw.replace("-0.90", "nan", 1),
                        raw.replace("5,1.5,,50,", "5,1.5,3,50,", 1),
                        raw.replace(raw.splitlines()[2] + "\n", "", 1)]:
            with self.assertRaises(ValueError):
                charm.parse_csv(changed.encode())

    def test_reference_registration_cannot_silently_enable_simulation(self):
        for change in [{"simulation": {"enabled": True}}, {"cards": {}},
                       {"reference": {"snapshot": "elsewhere.json"}}]:
            repo = self.temporary_repo()
            path = repo / "config/reference" / f"{charm.ID}.json"
            value = json.loads(path.read_text())
            for key, update in change.items():
                value.setdefault(key, {}).update(update)
            path.write_text(json.dumps(value))
            with self.assertRaises(ValueError):
                charm.descriptor(repo)

    def test_campaign_entrypoints_reject_before_dispatch_or_file_creation(self):
        self.assertIn(charm.ID, campaign.discover_all())
        self.assertNotIn(charm.ID, campaign.discover_pp_registry())
        with tempfile.TemporaryDirectory() as tmp, patch.object(campaign, "CAMPAIGN_ROOT", Path(tmp)), \
                patch.object(campaign, "prepare_pp") as prepare, patch.object(campaign, "run_pp") as run, \
                patch.object(campaign, "postprocess_pp") as post, patch.object(campaign, "plot_pp") as plot:
            for command in ["prepare", "campaign", "full", "postprocess", "plot"]:
                output = io.StringIO()
                with contextlib.redirect_stderr(output):
                    code = campaign.main([command, "--measurement", charm.ID, "--tag", "must-not-run"])
                self.assertEqual(code, 2)
                self.assertIn("massive-charm", output.getvalue())
            for function in [prepare, run, post, plot]:
                function.assert_not_called()
            self.assertEqual(list(Path(tmp).iterdir()), [])
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(campaign.main(["fetch-data", "--measurement", charm.ID]), 0)

    def test_browser_exports_data_without_claiming_a_campaign(self):
        value = charm.descriptor(ROOT)
        observation = {"available": False, "controllers": [], "worker_directories": []}
        item, files = browser.inspect_analysis(ROOT, [], charm.ID, "reference", value, observation, "both", False, None)
        self.assertEqual((item["availability"], item["reference_entries"]), ("reference", 45))
        self.assertEqual(len(item["plots"]), 9)
        self.assertIsNone(item["selected"])
        self.assertEqual(item["history"], [])
        self.assertFalse(item["simulation"]["enabled"])
        self.assertEqual(len(files), 2)
        self.assertTrue(all("data only" in p["labels"]["Title"] for p in item["plots"]))
        self.assertEqual(set(charm.plot_records(charm.validate(ROOT), "png")[0]["files"]), {"png"})
        with self.assertRaises(browser.common.ExportError):
            browser.inspect_analysis(ROOT, [], charm.ID, "reference", value, observation, "both", False, "fake")


if __name__ == "__main__":
    unittest.main()
