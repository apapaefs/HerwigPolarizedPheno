#!/usr/bin/env python3
"""Automatic derived galleries use complete, immutable input-addressed caches."""
from __future__ import annotations

import contextlib
import copy
import html
import io
import json
import re
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import hermes_apar_integrated as projection  # noqa: E402
import hermes_apar_integrated_campaign as automatic  # noqa: E402
import run_experimental_campaign as campaign  # noqa: E402


class HermesAutomaticProjectionTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.campaign = Path(self.temporary.name) / "campaign"
        postprocess = self.campaign / "postprocess"
        postprocess.mkdir(parents=True)
        self.summary = postprocess / "summary.json"
        self.summary.write_text(json.dumps({
            "measurement": "HERMES_2007_I726689", "tag": "toy", "bins": [],
        }) + "\n", encoding="utf-8")
        self.manifest = self.campaign / "manifest.json"
        self.manifest.write_text('{"generation_signature":"unchanged"}\n', encoding="utf-8")
        old = self.campaign / "derived-integrated-data" / "integrated.json"
        old.parent.mkdir()
        old.write_text("prior standalone reconstruction\n", encoding="utf-8")
        self.old_standalone = old

    def fake_reconstruct(self, summary, output, **options):
        output = Path(output)
        output.mkdir(parents=True, exist_ok=True)
        files = {name: "fake nonempty rendered artifact\n"
                 for name in automatic.REQUIRED_OUTPUTS}
        files.update({
            "integrated.json": '{"measurement":"HERMES_2007_I726689"}\n',
            "integrated.csv": "id,a_parallel\nP_X_01,0.1\n",
            "index.html": '<h1>Model-assisted Born projections</h1><a href="integrated.csv">CSV</a>\n',
            "integrated-overview-Q2GT1.png": "fake png bytes\n",
            "integrated-overview-Q2GT1.pdf": "fake pdf bytes\n",
        })
        for name, content in files.items():
            (output / name).write_text(content, encoding="utf-8")
        return {"measurement": "HERMES_2007_I726689"}

    def ensure(self, side_effect=None):
        return mock.patch.object(automatic, "reconstruct",
                                 side_effect=side_effect or self.fake_reconstruct)

    def test_first_generation_is_separate_and_repeat_reuses_complete_cache(self):
        originals = {path: path.read_bytes() for path in
                     (self.summary, self.manifest, self.old_standalone)}
        with self.ensure() as reconstruct:
            first = automatic.ensure_campaign_plots(self.campaign)
            second = automatic.ensure_campaign_plots(self.campaign)
        self.assertTrue(first["created"])
        self.assertFalse(second["created"])
        self.assertEqual(first["cache_key"], second["cache_key"])
        self.assertEqual(first["index"], second["index"])
        reconstruct.assert_called_once()
        self.assertEqual(reconstruct.call_args.kwargs, automatic.SETTINGS)
        self.assertTrue(reconstruct.call_args.kwargs["controls"])
        index = self.campaign / first["index"]
        self.assertTrue(index.is_file())
        self.assertIn("automatic", index.parts)
        for path, original in originals.items():
            self.assertEqual(path.read_bytes(), original)

    def test_changed_summary_creates_new_revision_and_preserves_prior_results(self):
        with self.ensure() as reconstruct:
            first = automatic.ensure_campaign_plots(self.campaign)
            old_index = self.campaign / first["index"]
            old_contents = old_index.read_bytes()
            self.summary.write_text(self.summary.read_text() + " \n", encoding="utf-8")
            second = automatic.ensure_campaign_plots(self.campaign)
        self.assertEqual(reconstruct.call_count, 2)
        self.assertTrue(second["created"])
        self.assertNotEqual(first["cache_key"], second["cache_key"])
        self.assertNotEqual(first["directory"], second["directory"])
        self.assertEqual(old_index.read_bytes(), old_contents)

    def test_damaged_output_is_preserved_and_gets_a_fresh_revision(self):
        with self.ensure() as reconstruct:
            first = automatic.ensure_campaign_plots(self.campaign)
            damaged = self.campaign / first["directory"] / "integrated.csv"
            damaged.write_text("damaged reconstruction\n", encoding="utf-8")
            second = automatic.ensure_campaign_plots(self.campaign)
        self.assertEqual(reconstruct.call_count, 2)
        self.assertTrue(second["created"])
        self.assertEqual(first["cache_key"], second["cache_key"])
        self.assertNotEqual(first["directory"], second["directory"])
        self.assertEqual(damaged.read_text(), "damaged reconstruction\n")

    def test_failed_reconstruction_never_becomes_reusable(self):
        def failed(summary, output, **options):
            output = Path(output)
            output.mkdir(parents=True, exist_ok=True)
            (output / "partial.txt").write_text("retain for diagnosis\n", encoding="utf-8")
            raise projection.ProjectionError("controlled rendering failure")

        with self.ensure(failed):
            with self.assertRaisesRegex(projection.ProjectionError, "controlled"):
                automatic.ensure_campaign_plots(self.campaign)
        partials = list(self.campaign.rglob("partial.txt"))
        self.assertEqual(len(partials), 1)
        with self.ensure() as reconstruct:
            recovered = automatic.ensure_campaign_plots(self.campaign)
            reused = automatic.ensure_campaign_plots(self.campaign)
        reconstruct.assert_called_once()
        self.assertTrue(recovered["created"])
        self.assertFalse(reused["created"])
        self.assertEqual(partials[0].read_text(), "retain for diagnosis\n")
        self.assertNotEqual(partials[0].parent,
                            self.campaign / recovered["directory"])

    def test_changed_fit_or_helper_source_gets_a_new_input_key(self):
        source = Path(self.temporary.name) / "fit-or-helper-source.txt"
        source.write_text("source revision one\n", encoding="utf-8")
        with mock.patch.object(automatic, "SOURCE_INPUTS", (source,)), self.ensure() as reconstruct:
            first = automatic.ensure_campaign_plots(self.campaign)
            source.write_text("source revision two\n", encoding="utf-8")
            second = automatic.ensure_campaign_plots(self.campaign)
        self.assertEqual(reconstruct.call_count, 2)
        self.assertNotEqual(first["cache_key"], second["cache_key"])
        record = json.loads((self.campaign / second["directory"] /
                             automatic.CACHE_MANIFEST).read_text())
        self.assertIn(str(source.resolve()), record["inputs"])
        self.assertEqual(record["settings"], automatic.SETTINGS)
        self.assertIn("fit-domain-sensitivity-Q2GT1.pdf", record["outputs"])

    def test_source_change_during_render_does_not_publish_a_completion_record(self):
        def changing(summary, output, **options):
            result = self.fake_reconstruct(summary, output, **options)
            self.summary.write_text(self.summary.read_text() + " \n", encoding="utf-8")
            return result

        with self.ensure(changing):
            with self.assertRaisesRegex(automatic.CacheError, "changed during"):
                automatic.ensure_campaign_plots(self.campaign)
        self.assertFalse(list(self.campaign.rglob(automatic.CACHE_MANIFEST)))
        self.assertTrue(list(self.campaign.rglob("integrated.csv")))

    def test_incomplete_renderer_outputs_do_not_publish_completion(self):
        def incomplete(summary, output, **options):
            self.fake_reconstruct(summary, output, **options)
            (Path(output) / "fit-domain-sensitivity-Q2GT1.pdf").unlink()

        with self.ensure(incomplete):
            with self.assertRaisesRegex(automatic.CacheError, "Incomplete"):
                automatic.ensure_campaign_plots(self.campaign)
        self.assertFalse(list(self.campaign.rglob(automatic.CACHE_MANIFEST)))

    def test_fiducial_change_requires_review_instead_of_new_wrong_acceptance_cache(self):
        source = Path(self.temporary.name) / "analysis.cc"
        source.write_text("different cuts\n", encoding="utf-8")
        with mock.patch.object(automatic, "FIDUCIAL_SOURCES", {source: "old pin"}), self.ensure() as reconstruct:
            with self.assertRaisesRegex(automatic.CacheError, "fiducial source changed"):
                automatic.ensure_campaign_plots(self.campaign)
        reconstruct.assert_not_called()
        self.assertFalse((self.old_standalone.parent / "automatic").exists())

    def test_missing_and_foreign_summary_are_refused_before_rendering(self):
        with self.ensure() as reconstruct:
            self.summary.unlink()
            with self.assertRaisesRegex(automatic.CacheError, "postprocess/summary.json"):
                automatic.ensure_campaign_plots(self.campaign)
            self.summary.write_text('{"measurement":"COMPASS_2017_I1501480"}\n', encoding="utf-8")
            with self.assertRaisesRegex(automatic.CacheError, "requires a HERMES"):
                automatic.ensure_campaign_plots(self.campaign)
        reconstruct.assert_not_called()

    def test_runner_gates_other_measurements_and_wraps_reconstruction_failure(self):
        with mock.patch.object(automatic, "ensure_campaign_plots", side_effect=AssertionError("unexpected")) as ensure:
            result = campaign._hermes_reconstructed_plot_gallery(
                self.campaign, {"id": "COMPASS_2017_I1501480"})
        self.assertIsNone(result)
        ensure.assert_not_called()
        with mock.patch.object(automatic, "ensure_campaign_plots", side_effect=automatic.CacheError("controlled")):
            with self.assertRaisesRegex(campaign.CampaignError, "HERMES Born.*controlled"):
                campaign._hermes_reconstructed_plot_gallery(
                    self.campaign, {"id": "HERMES_2007_I726689"})
        sentinel = {"index": "derived/index.html"}
        with mock.patch.object(automatic, "ensure_campaign_plots", return_value=sentinel) as ensure:
            self.assertIs(campaign._hermes_reconstructed_plot_gallery(
                self.campaign, {"id": "HERMES_2007_I726689"}), sentinel)
        ensure.assert_called_once_with(self.campaign)

    def test_dry_run_does_not_create_automatic_products_or_change_manifest(self):
        measurement = campaign.get_measurement("HERMES_2007_I726689")
        prediction = self.campaign / "postprocess/prediction.yoda"
        prediction.write_text("existing normalized prediction\n", encoding="utf-8")
        manifest = {
            "configuration": {"measurement_signature": campaign.measurement_signature(measurement)},
            "postprocess": {"yoda": "postprocess/prediction.yoda"},
            "runtime": {"tools": {"rivet-mkhtml": "/test/rivet-mkhtml"}},
        }
        campaign.atomic_write_json(self.manifest, manifest)
        original = self.manifest.read_bytes()
        args = campaign.make_parser().parse_args([
            "plot", "--measurement", measurement["id"], "--tag", "toy", "--dry-run",
        ])
        with mock.patch.object(campaign, "campaign_directory", return_value=self.campaign), mock.patch.object(
            automatic, "ensure_campaign_plots", side_effect=AssertionError("dry run generated plots")
        ) as ensure, contextlib.redirect_stdout(io.StringIO()):
            campaign.plot_campaign(args, measurement)
        ensure.assert_not_called()
        self.assertEqual(self.manifest.read_bytes(), original)
        self.assertFalse((self.old_standalone.parent / "automatic").exists())

    def test_root_and_analysis_gallery_links_resolve_and_distinguish_estimators(self):
        with self.ensure():
            supplemental = automatic.ensure_campaign_plots(self.campaign)
        output = self.campaign / "plots"
        analysis = output / "HERMES_2007_I726689"
        analysis.mkdir(parents=True)
        script = analysis / "AParallelQ2_Q2GT1.py"
        script.write_text("# original event-level Rivet curve\n", encoding="utf-8")
        (analysis / f"{script.stem}.png").write_bytes(b"original png")
        (analysis / f"{script.stem}.pdf").write_bytes(b"original pdf")
        root = campaign.write_plot_indexes(
            output, {"id": "HERMES_2007_I726689"}, [script], supplemental=supplemental)
        for index in (root, analysis / "index.html"):
            text = index.read_text()
            self.assertIn("Model-assisted experimental data", text)
            self.assertIn("nominal Herwig cell predictions", text)
            self.assertIn("event-level Monte Carlo integrals", text)
            self.assertIn("Errors are conditional on the central fit", text)
            self.assertIn("Fit-domain sensitivity", text)
            for target in re.findall(r'(?:href|src)="([^"]+)"', text):
                self.assertTrue((index.parent / html.unescape(target)).is_file(), target)
        self.assertIn(f"{script.stem}.png", (analysis / "index.html").read_text())
        self.assertEqual((analysis / f"{script.stem}.png").read_bytes(), b"original png")
        plain = campaign.write_plot_indexes(output, {"id": "other"}, [script])
        self.assertNotIn("Reconstructed Born", plain.read_text())

    def test_sparse_mc_masks_predictions_without_changing_data_weights_or_covariance(self):
        reference = json.loads(projection.REFERENCE.read_text())
        rows = []
        missing_identity = None
        for dataset in reference["datasets"]:
            for cell in dataset["points"]:
                rows.append({
                    "selection": dataset["selection"], "bin": cell["bin"],
                    "bin_low": cell["q2_low"], "bin_high": cell["q2_high"],
                    "supported": cell["supported_prediction"],
                    "a_parallel": .1 if cell["supported_prediction"] else None,
                    "a_parallel_stat": .02 if cell["supported_prediction"] else None,
                })
                if dataset["target"] == "P" and cell["global_bin"] == 13:
                    missing_identity = (dataset["selection"], cell["bin"])
        summary = {"measurement": projection.MEASUREMENT, "tag": "sparse", "bins": rows}

        def area(target, xl, xh, ql, qh):
            return (xh - xl) * (qh - ql)

        controls = {"shape": lambda target, xl, xh, ql, qh:
                    area(target, xl, xh, ql, qh) * (ql + qh)}
        complete = projection.build_projection_snapshot(reference, summary, area, controls=controls)
        sparse = copy.deepcopy(summary)
        missing = next(row for row in sparse["bins"]
                       if (row["selection"], row["bin"]) == missing_identity)
        missing["a_parallel"] = missing["a_parallel_stat"] = None
        result = projection.build_projection_snapshot(reference, sparse, area, controls=controls)
        missing_projections = available_projections = 0
        for target in ("P", "D"):
            expected = complete["targets"][target]
            actual = result["targets"][target]
            self.assertEqual(actual["weights"], expected["weights"])
            self.assertEqual(actual["statistical_covariance"], expected["statistical_covariance"])
            for row, normal, weights in zip(actual["bins"], expected["bins"], actual["weights"]):
                self.assertEqual(row["a_parallel"], normal["a_parallel"])
                self.assertEqual(row["stat"], normal["stat"])
                self.assertEqual(row["systematic_upper_bound"], normal["systematic_upper_bound"])
                if target == "P" and weights[12] > 0:
                    missing_projections += 1
                    self.assertFalse(row["mc_supported"])
                    self.assertEqual(row["mc_missing_global_bins"], [13])
                    self.assertIsNone(row["mc_a_parallel"])
                    self.assertIsNone(row["mc_stat_independent_cell_approx"])
                    self.assertIsNone(row["weight_controls"]["shape"]["mc_shift"])
                    self.assertIsNone(row["weight_controls"]["shape"]["data_minus_mc_shift"])
                    self.assertIsNotNone(row["weight_controls"]["shape"]["data_shift"])
                else:
                    available_projections += 1
                    self.assertTrue(row["mc_supported"])
                    self.assertEqual(row["mc_missing_global_bins"], [])
                    self.assertEqual(row["mc_a_parallel"], normal["mc_a_parallel"])
                    self.assertIsNotNone(row["weight_controls"]["shape"]["mc_shift"])
        self.assertGreater(missing_projections, 0)
        self.assertGreater(available_projections, 0)
        for value, error in ((None, .02), (.1, None), (float("inf"), .02)):
            invalid = copy.deepcopy(sparse)
            row = next(row for row in invalid["bins"]
                       if (row["selection"], row["bin"]) == missing_identity)
            row["a_parallel"], row["a_parallel_stat"] = value, error
            with self.assertRaises(projection.ProjectionError):
                projection.build_projection_snapshot(reference, invalid, area)


if __name__ == "__main__":
    unittest.main()
