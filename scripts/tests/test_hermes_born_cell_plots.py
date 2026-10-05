#!/usr/bin/env python3
"""Direct HERMES Born-cell atlases preserve data, masks and portable caches."""
from __future__ import annotations

import base64
import contextlib
import copy
import csv
import html
import json
import math
from pathlib import Path
import re
import shutil
import sys
import tempfile
import types
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import hermes_born_cell_plots as atlas  # noqa: E402


def fixtures():
    """Use the pinned table, without depending on a generated campaign."""
    reference = json.loads((ROOT / "data/experimental/HERMES_2007_I726689/"
                            "born-apar-reference.json").read_text(encoding="utf-8"))
    rows = []
    for dataset in reference["datasets"]:
        for point in dataset["points"]:
            supported = point["supported_prediction"]
            rows.append({
                "selection": dataset["selection"], "bin": point["bin"],
                "axis": "q2", "bin_low": point["q2_low"], "bin_high": point["q2_high"],
                "q2_low": point["q2_low"], "q2_high": point["q2_high"],
                "supported": supported,
                "a_parallel": ((.1 if dataset["target"] == "P" else .2)
                               + .002 * point["global_bin"] if supported else None),
                "a_parallel_stat": (.007 + .0001 * point["global_bin"] if supported else None),
                # Direct plots must use ratios already in the summary. These
                # arbitrary density columns cannot act as GD11 or other weights.
                "sigma_uu_pb": 9999999., "sigma_ll_pb": -9999999.,
            })
    return reference, {"measurement": atlas.MEASUREMENT, "tag": "toy", "bins": rows}


class HermesBornCellSnapshotTests(unittest.TestCase):
    def setUp(self):
        self.reference, self.summary = fixtures()
        self.snapshot = atlas.build_snapshot(self.summary, self.reference)

    def test_exact_published_cells_and_fixed_x_layout_for_both_targets(self):
        self.assertEqual(set(self.snapshot["targets"]), {"P", "D"})
        for target in ("P", "D"):
            block = self.snapshot["targets"][target]
            self.assertEqual(len(block["panels"]), 19)
            self.assertEqual(len(block["cells"]), 45)
            self.assertEqual([cell["global_bin"] for cell in block["cells"]], list(range(1, 46)))
            self.assertEqual([len(panel["cells"]) for panel in block["panels"]],
                             [1] * 4 + [2] * 4 + [3] * 11)
            self.assertEqual([panel["panel_index"] for panel in block["panels"]], list(range(19)))
            self.assertEqual([(panel["row_index"], panel["column_index"])
                              for panel in block["panels"]],
                             [(index // 5, index % 5) for index in range(19)])
            self.assertEqual([panel["selection"] for panel in block["panels"]],
                             [f"Born{target}_X{index:02d}" for index in range(1, 20)])
            self.assertEqual([cell["global_bin"] for panel in block["panels"]
                              for cell in panel["cells"]], list(range(1, 46)))

    def test_37_predictions_and_eight_low_q2_data_only_cells_are_preserved(self):
        expected = [6, 8, 10, 12] + list(range(13, 46))
        for target in ("P", "D"):
            cells = self.snapshot["targets"][target]["cells"]
            self.assertEqual([cell["global_bin"] for cell in cells
                              if cell["supported_prediction"]], expected)
            self.assertEqual(sum(cell["mc_supported"] for cell in cells), 37)
            unsupported = [cell for cell in cells if not cell["supported_prediction"]]
            self.assertEqual(len(unsupported), 8)
            for cell in unsupported:
                self.assertFalse(cell["mc_supported"])
                self.assertIsNone(cell["mc_a_parallel"])
                self.assertIsNone(cell["mc_stat"])
                self.assertIsNotNone(cell["value"])
                self.assertIsNotNone(cell["stat"])
                self.assertIsNotNone(cell["systematic_combined"])
                self.assertLess(cell["q2_low"], 1.)

    def test_published_values_errors_marker_coordinates_and_edges_are_unchanged(self):
        theory = {(row["selection"], row["bin"]): row for row in self.summary["bins"]}
        for dataset in self.reference["datasets"]:
            target = dataset["target"]
            cells = self.snapshot["targets"][target]["cells"]
            panel = self.snapshot["targets"][target]["panels"][int(dataset["id"][-2:]) - 1]
            for key in ("x_low", "x_high", "display_x_low", "display_x_high"):
                self.assertEqual(panel[key], dataset[key])
            for point in dataset["points"]:
                cell = cells[point["global_bin"] - 1]
                for key in ("value", "stat", "systematic_combined", "q2_mean", "q2_low",
                            "q2_high", "x_mean", "supported_prediction", "bin", "global_bin"):
                    self.assertEqual(cell[key], point[key])
                self.assertAlmostEqual(cell["total_error"],
                                       math.hypot(point["stat"], point["systematic_combined"]))
                mc = theory[(dataset["selection"], point["bin"])]
                self.assertEqual(cell["mc_a_parallel"], mc["a_parallel"])
                self.assertEqual(cell["mc_stat"], mc["a_parallel_stat"])
        # This explicit deuteron comparison catches a second .925 correction.
        self.assertAlmostEqual(self.snapshot["targets"]["D"]["cells"][5]["mc_a_parallel"], .212)

    def test_direct_snapshot_never_calls_gd11_and_does_not_mutate_inputs(self):
        reference, summary = copy.deepcopy(self.reference), copy.deepcopy(self.summary)
        with mock.patch.dict(sys.modules, {"hermes_unpolarized_fit": None}):
            result = atlas.build_snapshot(summary, reference)
        self.assertEqual(reference, self.reference)
        self.assertEqual(summary, self.summary)
        self.assertEqual(result, self.snapshot)
        self.assertNotIn("weights", result)
        self.assertNotIn("statistical_covariance", result["targets"]["P"])

    def test_non_born_rows_and_their_large_values_do_not_change_the_atlas(self):
        self.summary["bins"].append({"selection": "D_Q2GT1", "bin": 1,
                                     "a_parallel": 9999., "a_parallel_stat": 9999.})
        self.assertEqual(atlas.build_snapshot(self.summary, self.reference), self.snapshot)

    def test_theory_join_uses_cell_identity_instead_of_summary_order(self):
        summary = copy.deepcopy(self.summary)
        summary["bins"].reverse()
        self.assertEqual(atlas.build_snapshot(summary, self.reference), self.snapshot)

    def test_shared_row_limits_contain_zero_all_data_total_bars_and_mc_bars(self):
        self.assertEqual(len(self.snapshot["row_y_limits"]), 4)
        for row_index, (low, high) in enumerate(self.snapshot["row_y_limits"]):
            self.assertLess(low, 0.)
            self.assertGreater(high, 0.)
            for target in ("P", "D"):
                for panel in self.snapshot["targets"][target]["panels"]:
                    if panel["row_index"] != row_index:
                        continue
                    for cell in panel["cells"]:
                        self.assertLess(low, cell["value"] - cell["total_error"])
                        self.assertGreater(high, cell["value"] + cell["total_error"])
                        if cell["mc_supported"]:
                            self.assertLess(low, cell["mc_a_parallel"] - cell["mc_stat"])
                            self.assertGreater(high, cell["mc_a_parallel"] + cell["mc_stat"])

    def test_row_limit_responds_to_deuteron_data_and_proton_mc_errors_together(self):
        reference, summary = copy.deepcopy(self.reference), copy.deepcopy(self.summary)
        deuteron = next(dataset for dataset in reference["datasets"]
                        if dataset["selection"] == "BornD_X05")["points"][0]
        deuteron.update(value=-.5, stat=1., systematic_combined=2.)
        proton = next(row for row in summary["bins"]
                      if row["selection"] == "BornP_X05" and row["supported"])
        proton.update(a_parallel=2., a_parallel_stat=3.)
        changed = atlas.build_snapshot(summary, reference)
        low, high = changed["row_y_limits"][0]
        self.assertLess(low, -.5 - math.hypot(1., 2.))
        self.assertGreater(high, 5.)
        self.assertEqual(changed["row_y_limits"][1:], self.snapshot["row_y_limits"][1:])

    def test_sparse_supported_mc_cell_is_masked_without_changing_data_or_geometry(self):
        summary = copy.deepcopy(self.summary)
        missing = next(row for row in summary["bins"]
                       if row["selection"] == "BornP_X09" and row["bin"] == 2)
        missing.update(a_parallel=None, a_parallel_stat=None)
        snapshot = atlas.build_snapshot(summary, self.reference)
        cells = snapshot["targets"]["P"]["cells"]
        cell = next(cell for cell in cells if cell["selection"] == "BornP_X09" and cell["bin"] == 2)
        original = self.snapshot["targets"]["P"]["cells"][cell["global_bin"] - 1]
        self.assertTrue(cell["supported_prediction"])
        self.assertFalse(cell["mc_supported"])
        self.assertIsNone(cell["mc_a_parallel"])
        self.assertIsNone(cell["mc_stat"])
        self.assertTrue(cell["mc_missing_reason"])
        for key in ("value", "stat", "systematic_combined", "total_error",
                    "q2_mean", "q2_low", "q2_high", "x_low", "x_high"):
            self.assertEqual(cell[key], original[key])
        self.assertEqual(sum(cell["mc_supported"] for cell in cells), 36)
        self.assertEqual(sum(not cell["supported_prediction"] for cell in cells), 8)
        self.assertEqual(snapshot["targets"]["D"], self.snapshot["targets"]["D"])

    def test_nonfinite_inconsistent_and_negative_mc_uncertainties_are_rejected(self):
        for value, error in ((float("nan"), .1), (float("inf"), .1),
                             (.1, float("nan")), (.1, float("inf")),
                             (None, .1), (.1, None), (.1, -.1)):
            with self.subTest(value=value, error=error):
                summary = copy.deepcopy(self.summary)
                row = next(row for row in summary["bins"] if row["supported"])
                row.update(a_parallel=value, a_parallel_stat=error)
                with self.assertRaises(atlas.BornCellError):
                    atlas.build_snapshot(summary, self.reference)

    def test_wrong_measurement_missing_duplicate_and_mismatched_theory_cells_are_rejected(self):
        cases = []
        summary = copy.deepcopy(self.summary)
        summary["measurement"] = "OTHER"
        cases.append(("measurement", summary))
        summary = copy.deepcopy(self.summary)
        summary["bins"].pop(5)
        cases.append(("missing", summary))
        summary = copy.deepcopy(self.summary)
        summary["bins"].append(copy.deepcopy(summary["bins"][5]))
        cases.append(("duplicate", summary))
        summary = copy.deepcopy(self.summary)
        summary["bins"][5]["supported"] = False
        cases.append(("support", summary))
        summary = copy.deepcopy(self.summary)
        summary["bins"][5]["bin_high"] += .001
        cases.append(("physical boundary", summary))
        for reason, summary in cases:
            with self.subTest(reason=reason), self.assertRaises(atlas.BornCellError):
                atlas.build_snapshot(summary, self.reference)

    def test_invalid_reference_cell_count_support_and_uncertainties_are_rejected(self):
        for reason in ("missing", "global_bin", "support", "nonfinite", "negative"):
            with self.subTest(reason=reason):
                reference = copy.deepcopy(self.reference)
                point = reference["datasets"][0]["points"][0]
                if reason == "missing":
                    reference["datasets"][0]["points"].clear()
                elif reason == "global_bin":
                    point["global_bin"] = 2
                elif reason == "support":
                    point["supported_prediction"] = True
                elif reason == "nonfinite":
                    point["systematic_combined"] = float("nan")
                else:
                    point["stat"] = -.1
                with self.assertRaises(atlas.BornCellError):
                    atlas.build_snapshot(self.summary, reference)


class HermesBornCellCacheTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.campaign = Path(self.temporary.name) / "campaign"
        self.output = self.campaign / "plots"
        postprocess = self.campaign / "postprocess"
        postprocess.mkdir(parents=True)
        self.reference, summary = fixtures()
        self.summary = postprocess / "summary.json"
        self.summary.write_text(json.dumps(summary) + "\n", encoding="utf-8")
        self.manifest = self.campaign / "manifest.json"
        self.manifest.write_text('{"generation_signature":"unchanged"}\n', encoding="utf-8")

    @staticmethod
    def minimal_pdf():
        """Make a real one-page PDF for cheap cache tests without matplotlib."""
        objects = [b"<< /Type /Catalog /Pages 2 0 R >>",
                   b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
                   b"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 1 1] >>"]
        output, offsets = bytearray(b"%PDF-1.4\n"), [0]
        for number, contents in enumerate(objects, 1):
            offsets.append(len(output))
            output.extend(f"{number} 0 obj\n".encode() + contents + b"\nendobj\n")
        position = len(output)
        output.extend(b"xref\n0 4\n0000000000 65535 f \n")
        for offset in offsets[1:]:
            output.extend(f"{offset:010d} 00000 n \n".encode())
        output.extend(f"trailer\n<< /Size 4 /Root 1 0 R >>\nstartxref\n{position}\n%%EOF\n".encode())
        return bytes(output)

    def fake_render(self, snapshot, output):
        """Keep the real snapshot/cache code; replace only costly drawing."""
        output = Path(output)
        output.mkdir(parents=True, exist_ok=True)
        png = base64.b64decode("iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO+ip1sAAAAASUVORK5CYII=")
        for name in atlas.REQUIRED_OUTPUTS:
            path = output / name
            if path.suffix == ".pdf":
                path.write_bytes(self.minimal_pdf())
            elif path.suffix == ".png":
                path.write_bytes(png)
            elif path.suffix == ".json":
                path.write_text(json.dumps(snapshot) + "\n", encoding="utf-8")
            elif path.suffix == ".csv":
                with path.open("w", encoding="utf-8", newline="") as stream:
                    writer = csv.writer(stream)
                    writer.writerow(["target", "global_bin", "value", "mc_a_parallel"])
                    for target in ("P", "D"):
                        for cell in snapshot["targets"][target]["cells"]:
                            writer.writerow([target, cell["global_bin"], cell["value"], cell["mc_a_parallel"]])
            elif name == "index.html":
                links = [f'<a href="{html.escape(item)}">{html.escape(item)}</a>'
                         for item in sorted(atlas.REQUIRED_OUTPUTS) if item != "index.html"]
                path.write_text('<h1>Born A_parallel — published cells</h1>\n' + "\n".join(links),
                                encoding="utf-8")
        return [output / name for name in sorted(atlas.REQUIRED_OUTPUTS)]

    def stub_plotting(self):
        """Exercise real exports and HTML while replacing only drawing APIs."""
        matplotlib = types.ModuleType("matplotlib")
        matplotlib.use = mock.Mock()
        pyplot = types.ModuleType("matplotlib.pyplot")
        pyplot.rc_context = lambda options: contextlib.nullcontext()
        pyplot.close = mock.Mock()

        def subplots(*args, **kwargs):
            figure = mock.Mock()

            def savefig(path):
                path = Path(path)
                path.write_bytes(self.minimal_pdf() if path.suffix == ".pdf" else
                                 base64.b64decode("iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO+ip1sAAAAASUVORK5CYII="))

            figure.savefig.side_effect = savefig
            return figure, types.SimpleNamespace(flat=[mock.Mock() for _ in range(20)])

        pyplot.subplots = subplots
        lines = types.ModuleType("matplotlib.lines")
        lines.Line2D = mock.Mock()
        patches = types.ModuleType("matplotlib.patches")
        patches.Patch = mock.Mock()
        matplotlib.pyplot = pyplot
        return {"matplotlib": matplotlib, "matplotlib.pyplot": pyplot,
                "matplotlib.lines": lines, "matplotlib.patches": patches}

    def ensure(self):
        return atlas.ensure_campaign_plots(self.campaign, self.output)

    def source_bytes(self, directory):
        return {path.relative_to(directory).as_posix(): path.read_bytes()
                for path in directory.rglob("*") if path.is_file()}

    def assert_gallery_is_self_contained(self, directory):
        root = directory.resolve()
        self.assertTrue(list(directory.rglob("*.html")))
        for index in directory.rglob("*.html"):
            for link in re.findall(r'(?:href|src)="([^"]+)"',
                                   index.read_text(encoding="utf-8")):
                if link.startswith("#"):
                    continue
                self.assertNotRegex(link, r"^[a-zA-Z][a-zA-Z0-9+.-]*:")
                target = (index.parent / html.unescape(link).split("#", 1)[0]).resolve()
                self.assertTrue(target.is_file(), f"{index}: {link}")
                try:
                    target.relative_to(root)
                except ValueError:
                    self.fail(f"Gallery link leaves its copied directory: {index}: {link}")

    def test_initial_render_and_intact_reuse_preserve_campaign_numerical_inputs(self):
        originals = {path: path.read_bytes() for path in (self.summary, self.manifest)}
        with mock.patch.object(atlas, "render", side_effect=self.fake_render) as render:
            first, second = self.ensure(), self.ensure()
        render.assert_called_once()
        self.assertTrue(first["created"])
        self.assertFalse(second["created"])
        self.assertEqual(first["cache_key"], second["cache_key"])
        self.assertEqual(first["directory"], second["directory"])
        self.assertEqual(first["index"], second["index"])
        self.assertEqual(first["theory_family"], "nominal")
        self.assertEqual({figure["target"] for figure in first["figures"]}, {"P", "D"})
        for figure in first["figures"]:
            for suffix in ("pdf", "png"):
                path = self.campaign / figure[suffix]
                self.assertTrue(path.is_file())
                self.assertEqual(path.name, f"AParallel_{figure['target']}_BornCells_5x4.{suffix}")
        for path, original in originals.items():
            self.assertEqual(path.read_bytes(), original)
        record = json.loads((self.campaign / first["directory"] / atlas.CACHE_MANIFEST)
                            .read_text(encoding="utf-8"))
        self.assertTrue(atlas.REQUIRED_OUTPUTS.issubset(record["outputs"]))

    def test_changed_summary_creates_new_revision_without_overwriting_previous_plots(self):
        with mock.patch.object(atlas, "render", side_effect=self.fake_render) as render:
            first = self.ensure()
            previous = self.source_bytes(self.campaign / first["directory"])
            self.summary.write_text(self.summary.read_text(encoding="utf-8") + " \n", encoding="utf-8")
            second = self.ensure()
        self.assertEqual(render.call_count, 2)
        self.assertTrue(second["created"])
        self.assertNotEqual(first["cache_key"], second["cache_key"])
        self.assertNotEqual(first["directory"], second["directory"])
        self.assertEqual(self.source_bytes(self.campaign / first["directory"]), previous)

    def test_tampered_cache_is_retained_and_recovered_in_new_revision(self):
        with mock.patch.object(atlas, "render", side_effect=self.fake_render) as render:
            first = self.ensure()
            damaged = self.campaign / first["directory"] / "born-cells.csv"
            damaged.write_text("retain damaged revision\n", encoding="utf-8")
            second, third = self.ensure(), self.ensure()
        self.assertTrue(second["created"])
        self.assertFalse(third["created"])
        self.assertEqual(render.call_count, 2)
        self.assertEqual(first["cache_key"], second["cache_key"])
        self.assertNotEqual(first["directory"], second["directory"])
        self.assertEqual(second["directory"], third["directory"])
        self.assertEqual(damaged.read_text(encoding="utf-8"), "retain damaged revision\n")

    def test_malformed_cache_inventory_is_retained_and_rebuilt(self):
        with mock.patch.object(atlas, "render", side_effect=self.fake_render) as render:
            first = self.ensure()
            record_path = self.campaign / first["directory"] / atlas.CACHE_MANIFEST
            record = json.loads(record_path.read_text(encoding="utf-8"))
            record["outputs"] = sorted(atlas.REQUIRED_OUTPUTS)
            damaged = json.dumps(record) + "\n"
            record_path.write_text(damaged, encoding="utf-8")
            recovered = self.ensure()
        self.assertEqual(render.call_count, 2)
        self.assertTrue(recovered["created"])
        self.assertEqual(first["cache_key"], recovered["cache_key"])
        self.assertNotEqual(first["directory"], recovered["directory"])
        self.assertEqual(record_path.read_text(encoding="utf-8"), damaged)

    def test_modified_manifest_provenance_does_not_reuse_the_original_key(self):
        for field in ("inputs", "settings"):
            with self.subTest(field=field), mock.patch.object(atlas, "render", side_effect=self.fake_render) as render:
                self.output = self.campaign / f"plots-{field}"
                first = self.ensure()
                record_path = self.campaign / first["directory"] / atlas.CACHE_MANIFEST
                record = json.loads(record_path.read_text(encoding="utf-8"))
                record[field]["foreign-source"] = "unvalidated provenance"
                damaged = json.dumps(record) + "\n"
                record_path.write_text(damaged, encoding="utf-8")
                recovered = self.ensure()
                self.assertTrue(recovered["created"])
                self.assertEqual(render.call_count, 2)
                self.assertEqual(first["cache_key"], recovered["cache_key"])
                self.assertNotEqual(first["directory"], recovered["directory"])
                self.assertEqual(record_path.read_text(encoding="utf-8"), damaged)

    def test_input_change_during_cache_reuse_refuses_stale_revision_and_preserves_it(self):
        with mock.patch.object(atlas, "render", side_effect=self.fake_render) as render:
            first = self.ensure()
            cache = self.campaign / first["directory"]
            previous = self.source_bytes(cache)
            original_validation = atlas._valid_cache

            def changing_validation(directory, key):
                valid = original_validation(directory, key)
                if valid:
                    self.summary.write_text(self.summary.read_text(encoding="utf-8") + " \n",
                                            encoding="utf-8")
                return valid

            with mock.patch.object(atlas, "_valid_cache", side_effect=changing_validation):
                with self.assertRaisesRegex(atlas.BornCellError, "changed|retry"):
                    self.ensure()
            self.assertEqual(self.source_bytes(cache), previous)
            recovered = self.ensure()
        self.assertEqual(render.call_count, 2)
        self.assertTrue(recovered["created"])
        self.assertNotEqual(first["cache_key"], recovered["cache_key"])
        self.assertEqual(self.source_bytes(cache), previous)

    def test_failed_render_has_no_completion_manifest_and_retry_preserves_partial_output(self):
        def failed(snapshot, output):
            path = Path(output)
            path.mkdir(parents=True, exist_ok=True)
            (path / "partial.txt").write_text("retain rendering failure\n", encoding="utf-8")
            raise RuntimeError("controlled rendering failure")

        with mock.patch.object(atlas, "render", side_effect=failed):
            with self.assertRaisesRegex(RuntimeError, "controlled rendering failure"):
                self.ensure()
        partials = list(self.campaign.rglob("partial.txt"))
        self.assertEqual(len(partials), 1)
        self.assertFalse((partials[0].parent / atlas.CACHE_MANIFEST).exists())
        with mock.patch.object(atlas, "render", side_effect=self.fake_render) as render:
            recovered, reused = self.ensure(), self.ensure()
        render.assert_called_once()
        self.assertTrue(recovered["created"])
        self.assertFalse(reused["created"])
        self.assertEqual(partials[0].read_text(encoding="utf-8"), "retain rendering failure\n")
        self.assertNotEqual(partials[0].parent, self.campaign / recovered["directory"])

    def test_missing_required_output_never_becomes_reusable(self):
        def incomplete(snapshot, output):
            result = self.fake_render(snapshot, output)
            # Keep a nonempty diagnostic artifact instead of the requested PDF.
            missing = Path(output) / "AParallel_D_BornCells_5x4.pdf"
            missing.rename(missing.with_suffix(".partial"))
            return result

        with mock.patch.object(atlas, "render", side_effect=incomplete):
            with self.assertRaises(atlas.BornCellError):
                self.ensure()
        self.assertFalse(list(self.campaign.rglob(atlas.CACHE_MANIFEST)))
        retained = list(self.campaign.rglob("*.partial"))
        self.assertEqual(len(retained), 1)
        with mock.patch.object(atlas, "render", side_effect=self.fake_render):
            recovered = self.ensure()
        self.assertTrue(recovered["created"])
        self.assertTrue(retained[0].is_file())

    def test_gallery_is_portable_when_plots_and_analysis_are_copied_independently(self):
        with mock.patch.object(atlas, "render", side_effect=self.fake_render):
            info = self.ensure()
        index = self.campaign / info["index"]
        index.relative_to(self.output / atlas.MEASUREMENT)
        self.assert_gallery_is_self_contained(self.output)
        for directory in (self.output, self.output / atlas.MEASUREMENT):
            copied = Path(self.temporary.name) / f"copy-{directory.name}"
            shutil.copytree(directory, copied)
            self.assert_gallery_is_self_contained(copied)

    def test_real_exports_and_html_keep_all_cells_masks_and_shared_limits(self):
        with mock.patch.dict(sys.modules, self.stub_plotting()), mock.patch.object(atlas, "_panel") as panel:
            info = self.ensure()
        directory = self.campaign / info["directory"]
        snapshot = json.loads((directory / "born-cells.json").read_text(encoding="utf-8"))
        with (directory / "born-cells.csv").open(encoding="utf-8", newline="") as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 90)
        for target in ("P", "D"):
            target_rows = [row for row in rows if row["target"] == target]
            self.assertEqual([int(row["global_bin"]) for row in target_rows], list(range(1, 46)))
            self.assertEqual(sum(row["mc_supported"] == "True" for row in target_rows), 37)
            for row in target_rows:
                cell = snapshot["targets"][target]["cells"][int(row["global_bin"]) - 1]
                self.assertEqual(float(row["value"]), cell["value"])
                self.assertEqual(float(row["stat"]), cell["stat"])
                if row["supported_prediction"] == "False":
                    self.assertEqual((row["mc_a_parallel"], row["mc_stat"]), ("", ""))
        self.assertEqual(panel.call_count, 38)
        for call in panel.call_args_list:
            _, fixed_x_panel, limits = call.args
            self.assertEqual(limits, snapshot["row_y_limits"][fixed_x_panel["row_index"]])
        self.assert_gallery_is_self_contained(directory)
        self.assert_gallery_is_self_contained(self.output / atlas.MEASUREMENT)


if __name__ == "__main__":
    unittest.main()
