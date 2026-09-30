#!/usr/bin/env python3
"""Direct longitudinal asymmetries retain ordinary moments and axis semantics."""
from __future__ import annotations

import copy
import csv
import json
import io
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import run_experimental_campaign as campaign
from test_experimental_target_outputs import FAMILY, FAKE_YODA, mixed_measurement, load_series, tar_fixture
from test_experimental_order_coefficients import magnitude_inputs, independent_ratio_variance


def descriptor():
    result = mixed_measurement()
    result["diagnostics"] = []
    for selection, output in result["outputs"].items():
        output.pop("a1")
        output.update(estimator="a_parallel", observable="A_parallel", axis="q2",
                      projection="x-integrated", reference_status="theory-only")
        result["raw_observables"][selection] = {"ordinary": "ordinary"}
    return result


def calculate(measurement=None):
    selected = descriptor() if measurement is None else measurement
    with mock.patch.object(campaign, "_load_component_series", side_effect=load_series):
        return campaign._load_asymmetry_family_products({}, ROOT, selected, "nominal", FAMILY)[0]


class DirectAparallelTests(unittest.TestCase):
    def test_no_weighted_histograms_are_read_and_no_a1_result_exists(self):
        measurement = descriptor()
        with mock.patch.object(campaign, "_load_component_series", side_effect=load_series) as loader:
            results, _ = campaign._load_asymmetry_family_products({}, ROOT, measurement, "nominal", FAMILY)
        self.assertEqual({call.args[3] for call in loader.call_args_list}, {"ordinary"})
        for result in results.values():
            self.assertNotIn("a1_values", result)
            self.assertNotIn("a1_errors", result)
            self.assertNotIn("sigma_ll_over_d", result)
        self.assertAlmostEqual(results["P"]["apar_values"][0], .2)
        self.assertAlmostEqual(results["D"]["apar_values"][0], .4625/6.25)

    def test_target_ratio_error_matches_independent_poisson_gradients(self):
        measurement = descriptor()
        values = {"P": [14., 8., 8., 10.], "N": [2., 4., 3., 1.]}
        variances = {"P": [4., 1., 2., 3.], "N": [5., 2., 4., 6.]}
        target_uu = {"P": .5, "N": .5}
        target_ll = {"P": .4625, "N": .4625}
        signs = [1., -1., -1., 1.]
        denominator = sum(target_uu[c]*sum(v)/4. for c, v in values.items())
        numerator = sum(target_ll[c]*sum(s*v for s, v in zip(signs, row))/4.
                        for c, row in values.items())
        ratio = numerator/denominator
        variance = sum((target_ll[c]*sign/4.-ratio*target_uu[c]/4.)**2*var/denominator**2
                       for c in values for sign, var in zip(signs, variances[c]))
        deuteron = calculate(measurement)["D"]
        self.assertAlmostEqual(deuteron["apar_values"][0], ratio)
        self.assertAlmostEqual(deuteron["apar_errors"][0]**2, variance)
        # .925 is a numerator factor; the ratio is not the arithmetic P/N mean.
        self.assertNotAlmostEqual(ratio, .925*(.2-.4)/2.)

    def test_order_signs_and_variations_apply_to_direct_estimator(self):
        measurement = copy.deepcopy(campaign.get_measurement("HERMES_2007_I726689"))
        # This order fixture has one synthetic bin; published-cell geometry is
        # covered separately and must not be replaced by this invented bin.
        selections = ("Q2GT1", "Q2GT4", "D_Q2GT1", "D_Q2GT4")
        for block in ("outputs", "raw_observables"):
            measurement[block] = {key: measurement[block][key] for key in selections}
        for selection, output in measurement["outputs"].items():
            output.pop("a1")
            output["estimator"] = "a_parallel"
            measurement["raw_observables"][selection] = {
                "ordinary": measurement["raw_observables"][selection]["ordinary"]
            }
        family = campaign.campaign_family_specs(measurement, False)["nominal"]
        for member in (0, 1):
            groups, reader = magnitude_inputs(member=member)
            with mock.patch.object(campaign, "read_histogram_series", side_effect=reader):
                results, _ = campaign._load_asymmetry_family_products(
                    groups, ROOT, measurement, "nominal", family, member)
            for selection in ("Q2GT1", "D_Q2GT1"):
                ratio, variance, _ = independent_ratio_variance(selection, False)
                self.assertAlmostEqual(results[selection]["apar_values"][0], ratio)
                self.assertAlmostEqual(results[selection]["apar_errors"][0]**2, variance)

    def test_direct_objects_and_annotations_do_not_claim_a1_or_depolarization(self):
        measurement = descriptor()
        results = calculate(measurement)
        objects = campaign._asymmetry_family_objects(FAKE_YODA, measurement, {}, results, {}, {})
        self.assertEqual(len(objects), 10)
        self.assertFalse(any("a1" in item.path.lower() for item in objects))
        apar = next(item for item in objects if item.path == "/TEST/apar_D")
        self.assertEqual(apar.annotations["Observable"], "A_paralleld")
        self.assertEqual(apar.annotations["PlotAxis"], "q2_mean")
        self.assertEqual(apar.annotations["Projection"], "x-integrated")
        self.assertEqual(apar.annotations["ReferenceStatus"], "theory-only")
        self.assertIn("not used", apar.annotations["DepolarizationModel"])
        self.assertIn("no analysis-level A2 neglect", apar.annotations["A2G2Assumption"])

    def test_uncertainties_support_direct_and_mixed_estimators(self):
        direct = calculate()
        with mock.patch.object(campaign, "_load_component_series", side_effect=load_series):
            old, _ = campaign._load_asymmetry_family_products({}, ROOT, mixed_measurement(), "nominal", FAMILY)
        central = {"direct": direct["D"], "legacy": old["P"]}
        variations = {}
        for key, shift in (((0, 0, 1.), 0.), ((1, 0, 1.), .1), ((2, 0, 1.), -.1),
                           ((0, 0, .5), -.02), ((0, 0, 2.), .03)):
            varied = copy.deepcopy(central)
            for result in varied.values():
                result["apar_values"] = [value+shift for value in result["apar_values"]]
                if "a1_values" in result:
                    result["a1_values"] = [value+shift for value in result["a1_values"]]
            variations[key] = varied
        bands = campaign.aggregate_dis_uncertainties(variations)
        self.assertEqual(set(bands["direct"]), {"apar"})
        self.assertEqual(set(bands["legacy"]), {"apar", "a1"})
        self.assertAlmostEqual(bands["direct"]["apar"]["pdf_68"][0], math.sqrt(.02))
        self.assertAlmostEqual(bands["direct"]["apar"]["scale_down"][0], -.02)
        self.assertAlmostEqual(bands["direct"]["apar"]["scale_up"][0], .03)

    def test_unsupported_bins_mask_ratios_cross_sections_and_summary(self):
        measurement = descriptor()
        measurement["outputs"]["D"]["supported_bins"] = []
        results = calculate(measurement)
        self.assertEqual(results["D"]["apar_values"], [None])
        self.assertEqual(results["D"]["parity_pp_mm_values"], [None])
        row = campaign._summary_rows("D", results["D"])[0]
        self.assertFalse(row["supported"])
        self.assertIsNone(row["sigma_uu_pb"])
        self.assertIsNone(row["sigma_ll_stat_pb"])
        with mock.patch.object(FAKE_YODA.BinnedEstimate1D, "maskBin", create=True) as mask:
            objects = campaign._asymmetry_family_objects(FAKE_YODA, measurement, {}, results, {}, {})
            masked = [item for item in objects if item.annotations.get("Target") == "deuteron"]
        self.assertEqual(len(masked), 5)
        self.assertGreaterEqual(mask.call_count, 5)
        measurement["outputs"]["D"]["supported_bins"] = [2]
        with self.assertRaises(campaign.CampaignError):
            calculate(measurement)

    def test_direct00_uses_axis_and_supported_mask(self):
        measurement = descriptor()
        measurement["outputs"]["D"]["supported_bins"] = []
        family = {**FAMILY, "helicities": ["00"], "postprocess": "direct_unpolarized"}
        with mock.patch.object(campaign, "_load_component_series", side_effect=load_series):
            results, diagnostics = campaign._load_direct_unpolarized_products({}, ROOT, measurement, "direct", family)
        with mock.patch.object(FAKE_YODA.BinnedEstimate1D, "maskBin", create=True):
            _, closure = campaign._direct_unpolarized_objects(
                FAKE_YODA, measurement, results, diagnostics, calculate(measurement), {})
        rows = campaign._direct_unpolarized_rows(results, closure, measurement)
        self.assertEqual(closure["D"]["values"], [None])
        self.assertIn("q2_low", rows[0])
        self.assertNotIn("x_low", rows[0])
        self.assertIsNone(rows[1]["sigma_00_pb"])

    def test_axis_aware_mixed_summary_csv_has_no_invented_a1(self):
        results = calculate()
        direct_row = campaign._summary_rows("D", results["D"])[0]
        self.assertNotIn("a1", direct_row)
        self.assertNotIn("x_low", direct_row)
        self.assertEqual(direct_row["q2_low"], 0.)
        old = copy.deepcopy(results["P"])
        old["axis"] = "x"
        old.update(a1_values=[.4], a1_errors=[.03])
        rows = campaign._summary_rows("old", old)+[direct_row]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"summary.csv"
            campaign._write_summary_csv(path, rows)
            with path.open() as stream:
                parsed = list(csv.DictReader(stream))
        self.assertEqual(parsed[1]["a1"], "")
        self.assertEqual(parsed[0]["q2_low"], "")
        self.assertEqual(parsed[1]["x_low"], "")

    def test_estimator_axis_and_supported_bins_are_validated(self):
        campaign._validate_measurement(descriptor(), Path("TEST.json"))
        for field, value in (("estimator", "invalid"), ("axis", "W2"),
                             ("supported_bins", [0]), ("supported_bins", [True]),
                             ("supported_bins", [1, 1])):
            measurement = descriptor()
            measurement["outputs"]["P"][field] = value
            with self.assertRaises(campaign.CampaignError):
                campaign._validate_measurement(measurement, Path("TEST.json"))
        for location, key, value in (("outputs", "a1", "/TEST/fakeA1"),
                                      ("raw_observables", "weighted", "inverseD")):
            measurement = descriptor()
            measurement[location]["P"][key] = value
            with self.assertRaises(campaign.CampaignError):
                campaign._validate_measurement(measurement, Path("TEST.json"))


class SupplementaryReferenceTests(unittest.TestCase):
    def test_live_hermes_references_map_exactly_to_direct_output_paths(self):
        measurement = campaign.get_measurement("HERMES_2007_I726689")
        snapshot = campaign.load_reference_snapshot(measurement)
        self.assertEqual(len(snapshot["datasets"]), 40)
        born = [dataset for dataset in snapshot["datasets"] if dataset.get("selection", "").startswith("Born")]
        self.assertEqual(len(born), 38)
        self.assertEqual(sum(len(dataset["points"]) for dataset in born), 90)
        for dataset in born:
            output = measurement["outputs"][dataset["selection"]]
            self.assertEqual(output["estimator"], "a_parallel")
            self.assertEqual(output["axis"], "q2")
            self.assertEqual(output["apar"], dataset["rivet_path"].removeprefix("/REF"))
            self.assertEqual(output["supported_bins"], [point["bin"] for point in dataset["points"]
                                                       if point["supported_prediction"]])

    def files(self, root):
        primary = {"measurement": "TEST", "datasets": [
            {"id": "P", "rivet_path": "/REF/TEST/A1_P"}], "provenance": {"source": "A1"}}
        supplement = {"measurement": "TEST", "datasets": [
            {"id": "Born_X01_P", "rivet_path": "/REF/TEST/d07-x01-y02", "plot_axis": "q2_mean",
             "points": [{"x_mean": .02, "q2_mean": 1.5}]}], "provenance": {"source": "Born"}}
        for name, data in (("primary.json", primary), ("supplement.json", supplement)):
            (root/name).write_text(json.dumps(data))
        return {"id": "TEST", "reference": {"snapshot": "primary.json", "additional_snapshots": [
            {"path": "supplement.json", "sha256": campaign.sha256_file(root/"supplement.json")}]}}

    def test_supplement_is_checksum_verified_and_exact_path_overlay_is_selected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            measurement = self.files(root)
            with mock.patch.object(campaign, "DISPOL_ROOT", root):
                snapshot = campaign.load_reference_snapshot(measurement)
                self.assertEqual(len(snapshot["datasets"]), 2)
                points = campaign._fixed_reference_overlay_points(measurement, snapshot, "d07-x01-y02")
                self.assertEqual(points[0]["plot_x"], 1.5)
                self.assertEqual(snapshot["provenance"]["additional_references"], [{"source": "Born"}])
                (root/"supplement.json").write_text("{}")
                with self.assertRaisesRegex(campaign.CampaignError, "checksum mismatch"):
                    campaign.load_reference_snapshot(measurement)

    def test_duplicate_paths_ids_and_measurement_mismatch_are_rejected(self):
        for field, value in (("rivet_path", "/REF/TEST/A1_P"), ("id", "P"), ("measurement", "OTHER")):
            with tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                measurement = self.files(root)
                supplement = json.loads((root/"supplement.json").read_text())
                target = supplement if field == "measurement" else supplement["datasets"][0]
                target[field] = value
                (root/"supplement.json").write_text(json.dumps(supplement))
                measurement["reference"]["additional_snapshots"][0]["sha256"] = campaign.sha256_file(root/"supplement.json")
                with mock.patch.object(campaign, "DISPOL_ROOT", root), self.assertRaises(campaign.CampaignError):
                    campaign.load_reference_snapshot(measurement)

    def test_signature_hashes_supplementary_content_and_refuses_stale_checksum(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            measurement = self.files(root)
            measurement.update(analysis={"name": "TEST", "source": "analysis.cc", "info": "analysis.info", "plot": "analysis.plot"},
                               cards={"directory": "cards"})
            (root/"cards").mkdir()
            for path in ("analysis.cc", "analysis.info", "analysis.plot", "cards/input.in"):
                (root/path).write_text(path)
            with mock.patch.object(campaign, "DISPOL_ROOT", root):
                old = campaign.measurement_signature(measurement)
                (root/"supplement.json").write_text((root/"supplement.json").read_text()+"\n")
                with self.assertRaisesRegex(campaign.CampaignError, "checksum mismatch"):
                    campaign.measurement_signature(measurement)
                measurement["reference"]["additional_snapshots"][0]["sha256"] = campaign.sha256_file(root/"supplement.json")
                self.assertNotEqual(campaign.measurement_signature(measurement), old)

    def test_raw_source_files_are_verified_independently_of_snapshot_checksum(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            measurement = self.files(root)
            (root/"raw-data.json").write_text('{"value": 0.2}')
            supplement_path = root/"supplement.json"
            supplement = json.loads(supplement_path.read_text())
            supplement["source_files"] = [{"path": "raw-data.json", "sha256": campaign.sha256_file(root/"raw-data.json")}]
            supplement_path.write_text(json.dumps(supplement))
            measurement["reference"]["additional_snapshots"][0]["sha256"] = campaign.sha256_file(supplement_path)
            with mock.patch.object(campaign, "DISPOL_ROOT", root):
                self.assertEqual(len(campaign.load_reference_snapshot(measurement)["datasets"]), 2)
                (root/"raw-data.json").write_text('{"value": 0.3}')
                with self.assertRaisesRegex(campaign.CampaignError, "source file checksum mismatch"):
                    campaign.load_reference_snapshot(measurement)

    def test_supplementary_yoda_cache_contains_40_panels_without_editing_tracked_reference(self):
        try:
            yoda = campaign._import_yoda()
        except campaign.CampaignError:
            self.skipTest("Native YODA is required to verify the cached reference artifact")
        measurement = campaign.get_measurement("HERMES_2007_I726689")
        tracked = campaign.resolve_dispol_path(measurement["analysis"]["reference_yoda"])
        original = tracked.read_bytes()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            destination = campaign.ensure_reference_yoda(measurement, root)
            self.assertEqual(destination, root/"reference"/tracked.name)
            self.assertEqual(len(yoda.read(str(destination))), 40)
            # A stale local cache is regenerated solely from the pinned snapshots.
            destination.write_bytes(b"stale generated cache")
            campaign.ensure_reference_yoda(measurement, root)
            self.assertEqual(len(yoda.read(str(destination))), 40)
        self.assertEqual(tracked.read_bytes(), original)

    def test_cached_reference_directory_precedes_source_data_and_legacy_is_unchanged(self):
        measurement = campaign.get_measurement("HERMES_2007_I726689")
        root = Path("/private/tmp/hermes-reference-environment-test")
        with mock.patch.dict(campaign.os.environ, {"RIVET_DATA_PATH": "/previous/reference"}):
            environment = campaign.analysis_environment(measurement, root, {})
            self.assertEqual(environment["RIVET_DATA_PATH"].split(campaign.os.pathsep)[0], str(root/"reference"))
            legacy = campaign.get_measurement("COMPASS_2016_I1357198")
            old_environment = campaign.analysis_environment(legacy, root, {})
            self.assertNotIn(str(root/"reference"), old_environment["RIVET_DATA_PATH"].split(campaign.os.pathsep))
        with tempfile.TemporaryDirectory() as directory, mock.patch.object(campaign, "write_reference_yoda") as writer:
            destination = campaign.ensure_reference_yoda(legacy, Path(directory))
        self.assertEqual(destination, campaign.resolve_dispol_path(legacy["analysis"]["reference_yoda"]))
        writer.assert_not_called()

    def test_fetch_supplementary_reference_only_writes_the_ignored_data_cache(self):
        payload, sources, primary = tar_fixture()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            (root/"primary.json").write_text(json.dumps(primary))
            supplement = {"measurement": "TEST", "datasets": [{"id": "Born", "rivet_path": "/REF/TEST/Born"}]}
            (root/"supplement.json").write_text(json.dumps(supplement))
            tracked = root/"tracked.yoda.gz"
            tracked.write_bytes(b"original tracked reference bytes")
            measurement = {"id": "TEST", "analysis": {"reference_yoda": "tracked.yoda.gz"},
                           "reference": {"snapshot": "primary.json", "format": "tar-five-column-multidataset",
                                         "datasets": sources, "source_url": "https://example.test/source",
                                         "source_sha256": campaign.sha256_bytes(payload), "additional_snapshots": [
                                             {"path": "supplement.json", "sha256": campaign.sha256_file(root/"supplement.json")}]}}
            with mock.patch.object(campaign, "DISPOL_ROOT", root), \
                    mock.patch.object(campaign, "CAMPAIGN_ROOT", root/"campaigns"), \
                    mock.patch.object(campaign.urllib.request, "urlopen", return_value=io.BytesIO(payload)), \
                    mock.patch.object(campaign, "write_reference_yoda") as writer:
                campaign.fetch_reference_data(measurement)
            self.assertEqual(writer.call_args.args[1], root/"campaigns/_data_cache/TEST/reference/tracked.yoda.gz")
            self.assertEqual(len(writer.call_args.args[0]["datasets"]), 3)
            self.assertEqual(tracked.read_bytes(), b"original tracked reference bytes")


class PublishedCoordinateTests(unittest.TestCase):
    def render(self, reference_values=(.1, .2), reference_points=None):
        import numpy
        points = reference_points or [
            {"plot_x": 1.3, "value": .1, "stat": .01},
            {"plot_x": 3.1, "value": .2, "stat": .02}]
        data = {"yvals": {"Data": list(reference_values), "MC": [.12, .21]},
                "xpoints": {"Data": [1.5, 3.], "MC": [1.5, 3.]},
                "xedges": {"Data": [1., 2., 4.], "MC": [1., 2., 4.]},
                "xerrs": {"Data": [[.5, 1.], [.5, 1.]], "MC": [[.5, 1.], [.5, 1.]]},
                "ref_xerrs": [[.5, 1.], [.5, 1.]]}
        axis = mock.Mock()
        source = """# curve from input yoda files in main panel
ax.errorbar(dataf['xpoints']['Data'], dataf['yvals']['Data'], label='total')
ax.errorbar(dataf['xpoints']['MC'], dataf['yvals']['MC'], label='MC')
# set plot metadata as defined above
"""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"plot.py"
            path.write_text(source)
            campaign.apply_published_reference_coordinates(path, points)
            campaign.add_experimental_error_overlay(path, points, show_statistical=True)
            exec(path.read_text(), {"dataf": data, "np": numpy, "ax": axis})
        return data, axis

    def test_total_and_statistical_bars_use_means_without_moving_theory(self):
        data, axis = self.render()
        self.assertEqual(data["xpoints"]["Data"], [1.3, 3.1])
        self.assertEqual(data["xpoints"]["MC"], [1.5, 3.])
        self.assertEqual(data["xedges"]["MC"], [1., 2., 4.])
        self.assertEqual(data["ref_xerrs"], data["xerrs"]["Data"])
        for actual, expected in zip(data["xerrs"]["Data"], ([.3, 1.1], [.7, .9])):
            for value, target in zip(actual, expected):
                self.assertAlmostEqual(value, target)
        self.assertEqual(axis.errorbar.call_args_list[0].args[0], [1.3, 3.1])
        self.assertEqual(axis.errorbar.call_args_list[2].args[0], [1.3, 3.1])

    def test_mismatched_reference_curve_is_rejected(self):
        with self.assertRaisesRegex(RuntimeError, "do not match"):
            self.render(reference_values=(.2, .1))
        with self.assertRaisesRegex(RuntimeError, "outside"):
            self.render(reference_points=[{"plot_x": .9, "value": .1, "stat": .01},
                                          {"plot_x": 3.1, "value": .2, "stat": .02}])

    def test_prediction_only_script_remains_unchanged(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"plot.py"
            path.write_text("prediction only\n")
            campaign.apply_published_reference_coordinates(path, None)
            self.assertEqual(path.read_text(), "prediction only\n")

    def test_canvas_layout_is_drawn_before_each_format_save(self):
        events = []
        figure = mock.Mock()
        figure.canvas.draw.side_effect = lambda: events.append("draw")
        pyplot = mock.Mock()
        pyplot.savefig.side_effect = lambda path: events.append(path)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"plot.py"
            path.write_text("plt.savefig('plot.pdf')\nplt.savefig('plot.png')\n")
            campaign.ensure_plot_canvas_draw(path)
            exec(path.read_text(), {"fig": figure, "plt": pyplot})
        self.assertEqual(events, ["draw", "plot.pdf", "draw", "plot.png"])

    def test_mathtext_fallback_uses_matching_builtin_fonts(self):
        for tex_enabled in (False, True):
            pyplot = mock.Mock()
            pyplot.rcParams = {"mathtext.fontset": "custom"}
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory)/"plot.py"
                path.write_text(f"plt.rcParams['text.usetex'] = {tex_enabled!r}\nplt.savefig('plot.png')\n")
                campaign.ensure_plot_canvas_draw(path)
                exec(path.read_text(), {"fig": mock.Mock(), "plt": pyplot})
            self.assertEqual(pyplot.rcParams["mathtext.fontset"], "custom" if tex_enabled else "dejavusans")

    def test_statistical_note_avoids_a_lower_legend(self):
        axis = mock.Mock()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"plot.py"
            path.write_text("ax.legend([], loc='lower center')\n# set plot metadata as defined above\n")
            campaign.add_experimental_error_overlay(path, [{"plot_x": 1.3, "value": .1, "stat": .01}],
                                                    show_statistical=True)
            exec(path.read_text(), {"ax": axis})
        self.assertEqual(axis.text.call_args.args[:2], (.98, .97))
        self.assertEqual(axis.text.call_args.kwargs["ha"], "right")
        self.assertEqual(axis.text.call_args.kwargs["va"], "top")


if __name__ == "__main__":
    unittest.main()
