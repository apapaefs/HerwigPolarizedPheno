"""Helium additions must compose with the current shared DIS interfaces."""
from __future__ import annotations

import contextlib
import io
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import run_experimental_campaign as campaign
import run_phenomenology_campaign as primary
import phenomenology_reference_data as references
import rivet_mkhtml_safe as safe
from test_dis_plot_scale_nbaldona import Estimate
import test_dis_hide_scale_plot_nbaldona as plot_fixtures


class SharedIntegrationTests(unittest.TestCase):
    def test_primary_forwards_dis_sign_and_plot_options(self):
        for command in ("plot", "full"):
            for options in (["--hide-scale-uncertainty"], ["--plot-scale", "0.5"]):
                label = r"$\frac{M-D}{\sigma_D}$"
                args = primary.make_parser().parse_args([
                    command, "--measurement", "E142_1996_I424108", "--tag", "test",
                    "--ratiolabel", label, *options,
                ])
                forwarded = campaign.make_parser().parse_args(primary._legacy_arguments(args))
                self.assertEqual(forwarded.ratiolabel, label)
                self.assertEqual(forwarded.plot_scale, args.plot_scale)
                self.assertEqual(forwarded.hide_scale_uncertainty, args.hide_scale_uncertainty)
        for command in ("postprocess", "full", "check-nlo-signs"):
            args = primary.make_parser().parse_args([
                command, "--measurement", "E142_1996_I424108", "--tag", "test",
                "--nlo-sign-convention", "magnitudes",
            ])
            forwarded = campaign.make_parser().parse_args(primary._legacy_arguments(args))
            self.assertEqual(forwarded.nlo_sign_convention, "magnitudes")

    def test_dis_only_options_are_not_silently_ignored_for_pp(self):
        for command, options in (("plot", ["--plot-scale", "2"]),
                                 ("plot", ["--hide-scale-uncertainty"]),
                                 ("postprocess", ["--nlo-sign-convention", "signed"]),
                                 ("check-nlo-signs", [])):
            with mock.patch.object(primary, "discover_all", return_value={
                "PP": {"process_kind": "polarized_pp"}
            }), contextlib.redirect_stderr(io.StringIO()) as output:
                self.assertEqual(primary.main([
                    command, "--measurement", "PP", "--tag", "test", *options,
                ]), 2)
            self.assertIn("only valid for fixed-target DIS", output.getvalue())

    def test_current_e142_descriptor_allows_proven_plot_metadata_refresh(self):
        measurement = campaign.get_measurement("E142_1996_I424108")
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            campaign.atomic_write_json(directory / "resolved-measurement.json", measurement)
            expected = {"mode": "presentation_only_plot_metadata_refresh"}
            with mock.patch.object(campaign, "authorize_e142_output_schema_migration") as migrate, \
                    mock.patch.object(campaign, "authorize_plot_metadata_refresh", return_value=expected) as refresh:
                actual = campaign._assert_manifest_signatures_current(
                    {"configuration": {"measurement_signature": "historical-plot"}},
                    measurement, campaign_dir=directory, allow_plot_metadata_refresh=True,
                )
            self.assertEqual(actual, expected)
            refresh.assert_called_once()
            migrate.assert_not_called()

    def test_published_coordinate_aliases_and_upstream_default(self):
        point = {"x_mean": .2, "q2_mean": 3.5}
        for axis_fields, expected in (({}, 3.5), ({"plot_axis": "x_mean"}, .2),
                                      ({"plot_coordinate": "x_mean"}, .2)):
            snapshot = {"datasets": [{"rivet_path": "/REF/TEST/d01", "rivet_aliases": ["/REF/TEST/d02"],
                                       "points": [point], **axis_fields}]}
            for stem in ("d01", "d02"):
                result = campaign._fixed_reference_overlay_points({}, snapshot, stem)
                self.assertEqual(result[0]["plot_x"], expected)
        with self.assertRaisesRegex(campaign.CampaignError, "Conflicting"):
            campaign._reference_dataset_axis({"plot_axis": "x_mean", "plot_coordinate": "q2_mean"})

    def test_aliased_references_use_writer_that_preserves_aliases(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            generated = root / "generated.yoda"
            generated.write_bytes(b"reference including nested-window aliases")
            snapshot = {"measurement": "E142_1996_I424108", "datasets": [
                {"bin_edges": [0., 1.], "rivet_aliases": ["/REF/TEST/alias"]}
            ]}
            destination = root / "copied.yoda"
            with mock.patch.object(references, "write_reference_yoda", return_value=generated) as write, \
                    mock.patch.object(campaign, "_import_yoda") as native:
                campaign.write_reference_yoda(snapshot, destination)
            write.assert_called_once_with(snapshot["measurement"], snapshot)
            native.assert_not_called()
            self.assertEqual(destination.read_bytes(), generated.read_bytes())

    def test_selected_scale_handles_direct_apar_and_summary_only_outputs(self):
        measurement = {"analysis": {"name": "TEST"}, "outputs": {
            "direct": {"apar": "/TEST/apar"},
            "summary": {"a1": "/TEST/not-written", "summary_only": True},
        }}
        path = "/TEST/VARIATIONS/p000-u000-mu2/apar_direct"
        with mock.patch.object(campaign, "_import_yoda", return_value=types.SimpleNamespace(
            read=lambda _: {path: Estimate(path)}
        )), mock.patch.object(campaign, "_write_yoda_objects") as write:
            bands = campaign._write_selected_scale_plot_copy(measurement, Path("in"), Path("out"), 2.)
        self.assertEqual(bands, {"apar": {"monte_carlo": [.12, None, .12]}})
        self.assertEqual([item.path for item in write.call_args.args[1]], ["/TEST/apar"])

    def test_hermes_galleries_retained_only_for_default_display(self):
        for options in ([], ["--hide-scale-uncertainty"], ["--plot-scale", "1"]):
            with self.subTest(options=options), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                measurement, _ = plot_fixtures.HideScalePlotTests().fixture(root)
                measurement["id"] = "HERMES_2007_I726689"
                measurement["physics"] = {"acceptance": True}
                suffix = "-central" if "--hide-scale-uncertainty" in options else "-mu1" if options else ""
                output = root / f"plots{suffix}"
                args = campaign.make_parser().parse_args([
                    "plot", "--measurement", measurement["id"], "--tag", "test", *options,
                ])
                def generate(*args):
                    output.mkdir(exist_ok=True)
                    (output / "d01.py").write_text("plt.savefig('d01.png')\n")
                    (output / "d01.png").write_bytes(b"plot")
                def index(*args, **kwargs):
                    destination = output / "index.html"
                    destination.write_text("index")
                    return destination
                with contextlib.ExitStack() as stack:
                    stack.enter_context(contextlib.redirect_stdout(io.StringIO()))
                    stack.enter_context(mock.patch.object(campaign, "DISPOL_ROOT", root))
                    for name, replacement in {
                        "campaign_directory": mock.Mock(return_value=root),
                        "_assert_manifest_signatures_current": mock.Mock(return_value=None),
                        "ensure_reference_yoda": mock.Mock(),
                        "analysis_environment": mock.Mock(return_value={}),
                        "write_campaign_monitor_files": mock.Mock(),
                        "build_campaign_monitor_payload": mock.Mock(),
                        "existing_herwig_monitor": mock.Mock(),
                        "_write_no_scale_plot_copy": mock.Mock(),
                        "_write_selected_scale_plot_copy": mock.Mock(return_value={}),
                        "_run_logged": mock.Mock(side_effect=generate),
                        "add_experimental_error_overlay": mock.Mock(),
                        "add_theory_uncertainty_overlay": mock.Mock(),
                        "plot_script_has_finite_y": mock.Mock(return_value=True),
                        "write_plot_indexes": mock.Mock(side_effect=index),
                    }.items():
                        stack.enter_context(mock.patch.object(campaign, name, replacement))
                    galleries = [stack.enter_context(mock.patch.object(campaign, name, return_value={"test": name}))
                                 for name in ("_hermes_reconstructed_plot_gallery", "_hermes_born_cell_plot_gallery",
                                              "_hermes_data_ratio_plot_gallery")]
                    import hermes_acceptance_control
                    galleries.append(stack.enter_context(mock.patch.object(hermes_acceptance_control, "ensure_gallery", return_value={})))
                    stack.enter_context(mock.patch.object(campaign.subprocess, "run", return_value=types.SimpleNamespace(returncode=0)))
                    campaign.plot_campaign(args, measurement)
                for gallery in galleries:
                    self.assertEqual(gallery.call_count, int(not options))
                manifest = campaign.load_json(root / campaign.MANIFEST_NAME)
                if options:
                    self.assertEqual(manifest["plots"], {"index": "plots/index.html"})
                    self.assertIn("supplemental_galleries", manifest["plots" + suffix.replace("-", "_")])
                else:
                    self.assertIn("reconstructed_born_projections", manifest["plots"])

    def test_yoda_adapter_preserves_mode_for_both_api_generations(self):
        for modern in (False, True):
            for call_options, expected in (({}, False), ({"deviation": True}, True),
                                            ({"deviation_modes": [True, False]}, [True, False])):
                module = types.ModuleType("yoda.plotting.fetch_data")
                captured = []
                if modern:
                    def original(aos, panels, *, deviation_modes=False, **kwargs):
                        self.assertNotIn("deviation", kwargs)
                        captured.append(deviation_modes)
                else:
                    def original(aos, panels, *, deviation=False, **kwargs):
                        self.assertNotIn("deviation_modes", kwargs)
                        captured.append(deviation)
                module.mkCurves1D = original
                yoda = types.ModuleType("yoda")
                yoda.__path__ = []
                plotting = types.ModuleType("yoda.plotting")
                plotting.__path__ = []
                yoda.plotting = plotting
                plotting.fetch_data = module
                with mock.patch.dict(sys.modules, {"yoda": yoda, "yoda.plotting": plotting,
                                                    "yoda.plotting.fetch_data": module}):
                    safe._patch_fetch_data_band_annotations()
                    module.mkCurves1D({}, [], **call_options)
                self.assertEqual(captured, [expected])


if __name__ == "__main__":
    unittest.main()
