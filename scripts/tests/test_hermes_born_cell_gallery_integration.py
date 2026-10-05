"""Direct Born panels belong to both portable galleries and only real HERMES plots."""
from __future__ import annotations

import contextlib
import html
import io
import json
import re
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import hermes_born_cell_plots as panels
import run_experimental_campaign as campaign

MEASUREMENT = "HERMES_2007_I726689"


class BornCellGalleryIntegrationTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.campaign = self.root / "campaign"
        self.output = self.campaign / "plots"
        self.analysis = self.output / MEASUREMENT
        self.analysis.mkdir(parents=True)

    def fixture(self):
        gallery = self.analysis / "born-cell-panels" / "revision"
        gallery.mkdir(parents=True)
        for target in ("P", "D"):
            for suffix in ("pdf", "png"):
                (gallery / f"AParallel_{target}_BornCells_5x4.{suffix}").write_bytes(b"rendered")
        for name in ("index.html", "born-cells.json", "born-cells.csv"):
            (gallery / name).write_text("recorded data\n", encoding="utf-8")
        script = self.analysis / "d07-x09-y02.py"
        script.write_text("# Existing Rivet script\n", encoding="utf-8")
        (script.with_suffix(".png")).write_bytes(b"Original Rivet PNG")
        (script.with_suffix(".pdf")).write_bytes(b"Original Rivet PDF")
        return {"directory": gallery.relative_to(self.campaign).as_posix()}, script

    def assert_portable(self, root):
        for index in root.rglob("*.html"):
            for link in re.findall(r'(?:href|src)="([^"]+)"', index.read_text(encoding="utf-8")):
                path = (index.parent / html.unescape(link)).resolve()
                self.assertTrue(path.is_relative_to(root.resolve()), (index, link))
                self.assertTrue(path.is_file(), (index, link))

    def test_both_indexes_embed_separate_figures_without_changing_original_plots(self):
        info, script = self.fixture()
        originals = {path: path.read_bytes() for path in self.analysis.glob("d07-x09-y02.*")}
        reconstructed = self.analysis / "reconstructed-born" / "original"
        reconstructed.mkdir(parents=True)
        for selection in ("Q2GT1", "Q2GT4"):
            for axis in ("x", "q2"):
                for target in ("P", "D"):
                    for suffix in ("pdf", "png"):
                        (reconstructed / f"AParallel_{target}_vs_{axis}_{selection}_reconstructed.{suffix}").write_bytes(b"Existing reconstruction")
        for name in ("index.html", "integrated.csv", "integrated.json",
                     "fit-domain-sensitivity-Q2GT1.pdf", "integrated-overview-Q2GT1.pdf",
                     "integrated-overview-Q2GT1.png"):
            (reconstructed / name).write_bytes(b"Existing reconstruction")
        with mock.patch("hermes_apar_integrated_campaign.publish_gallery", return_value=reconstructed):
            index = campaign.write_plot_indexes(self.output, {"id": MEASUREMENT}, [script],
                                              supplemental={"directory": "unused"}, born_cells=info)
        for page in (index, self.analysis / "index.html"):
            text = page.read_text(encoding="utf-8")
            for target in ("P", "D"):
                self.assertIn(f"AParallel_{target}_BornCells_5x4.png", text)
                self.assertIn(f"AParallel_{target}_BornCells_5x4.pdf", text)
            self.assertIn("Reconstructed Born A_parallel projections", text)
            self.assertIn("no unpolarized-fit averaging", text)
        for path, payload in originals.items():
            self.assertEqual(path.read_bytes(), payload)
        for source, name in ((self.output, "copied-plots"), (self.analysis, "copied-analysis")):
            destination = self.root / name
            shutil.copytree(source, destination)
            self.assert_portable(destination)

    def test_runner_gates_other_measurements_and_reports_render_failure(self):
        with mock.patch.object(panels, "ensure_campaign_plots") as ensure:
            self.assertIsNone(campaign._hermes_born_cell_plot_gallery(
                self.campaign, self.output, {"id": "COMPASS_2017_I1501480"}))
        ensure.assert_not_called()
        sentinel = {"directory": "plots/HERMES/born-cell-panels/key"}
        with mock.patch.object(panels, "ensure_campaign_plots", return_value=sentinel) as ensure:
            self.assertIs(campaign._hermes_born_cell_plot_gallery(
                self.campaign, self.output, {"id": MEASUREMENT}), sentinel)
        ensure.assert_called_once_with(self.campaign, self.output)
        with mock.patch.object(panels, "ensure_campaign_plots", side_effect=RuntimeError("controlled failure")):
            with self.assertRaisesRegex(campaign.CampaignError, "HERMES Born cell.*controlled failure"):
                campaign._hermes_born_cell_plot_gallery(self.campaign, self.output, {"id": MEASUREMENT})

    def test_dry_run_does_not_render_either_supplement_or_change_manifest(self):
        measurement = campaign.get_measurement(MEASUREMENT)
        postprocess = self.campaign / "postprocess"
        postprocess.mkdir()
        (postprocess / "prediction.yoda").write_text("Normalized prediction\n", encoding="utf-8")
        manifest = {
            "configuration": {"measurement_signature": campaign.measurement_signature(measurement)},
            "postprocess": {"yoda": "postprocess/prediction.yoda"},
            "runtime": {"tools": {"rivet-mkhtml": "/test/rivet-mkhtml"}},
        }
        manifest_path = self.campaign / "manifest.json"
        campaign.atomic_write_json(manifest_path, manifest)
        original = manifest_path.read_bytes()
        args = campaign.make_parser().parse_args([
            "plot", "--measurement", MEASUREMENT, "--tag", "toy", "--dry-run",
        ])
        with mock.patch.object(campaign, "campaign_directory", return_value=self.campaign), \
                mock.patch.object(campaign, "_hermes_born_cell_plot_gallery") as born, \
                mock.patch.object(campaign, "_hermes_reconstructed_plot_gallery") as integrated, \
                contextlib.redirect_stdout(io.StringIO()):
            campaign.plot_campaign(args, measurement)
        born.assert_not_called()
        integrated.assert_not_called()
        self.assertEqual(manifest_path.read_bytes(), original)

    def test_gallery_refuses_missing_outputs_and_outside_paths(self):
        info, script = self.fixture()
        wrong = {"directory": "../outside"}
        with self.assertRaisesRegex(campaign.CampaignError, "outside"):
            campaign.write_plot_indexes(self.output, {"id": MEASUREMENT}, [script], born_cells=wrong)
        (self.campaign / info["directory"] / "AParallel_D_BornCells_5x4.png").unlink()
        with self.assertRaisesRegex(campaign.CampaignError, "missing or empty"):
            campaign.write_plot_indexes(self.output, {"id": MEASUREMENT}, [script], born_cells=info)


if __name__ == "__main__":
    unittest.main()
