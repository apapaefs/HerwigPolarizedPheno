"""Native lower-panel errors must reach the rendered Matplotlib artists."""
from __future__ import annotations

import contextlib
import importlib
import io
import os
from pathlib import Path
import runpy
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import run_experimental_campaign as campaign


class NativeResidualRenderingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tool = shutil.which("rivet-mkhtml")
        if cls.tool is None:
            raise unittest.SkipTest("Native rivet-mkhtml is not on PATH")
        cls.temporary = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.root = Path(cls.temporary.name).resolve()
        environment = mock.patch.dict(os.environ, {
            "MPLCONFIGDIR": str(cls.root / "matplotlib-cache"),
            "MPLBACKEND": "Agg",
            "DISPOL_FORCE_NO_ERROR_BANDS": "1",
        })
        environment.start()
        cls.addClassCleanup(environment.stop)
        try:
            cls.yoda = importlib.import_module("yoda")
            importlib.import_module("rivet")
            cls.np = importlib.import_module("numpy")
            matplotlib = importlib.import_module("matplotlib")
            matplotlib.use("Agg")
            cls.plt = importlib.import_module("matplotlib.pyplot")
            cls.errorbar_type = importlib.import_module("matplotlib.container").ErrorbarContainer
        except ImportError as exc:
            raise unittest.SkipTest(f"Native plotting runtime is unavailable: {exc}") from exc
        if not hasattr(cls.yoda, "Estimate1D"):
            raise unittest.SkipTest("The native plotting check requires YODA Estimate1D")

    def render(self, mode, *, residual):
        directory = self.root / (mode + ("-residual" if residual else "-ratio"))
        directory.mkdir()
        path = "/MC_RESIDUAL_TEST/d01-x01-y01"
        edges = [0., 1., 2.]
        reference = campaign._estimate_from_values(
            self.yoda, edges, "/REF" + path, [.2, .4], [.05, .1], {},
        )
        prediction = campaign._estimate_from_values(
            self.yoda, edges, path, [.3, .6], [.03, .06], {},
        )
        for index, factor in enumerate((1., 2.), start=1):
            prediction.bin(index).setErr(-.04 * factor, .04 * factor, "pdf68")
            prediction.bin(index).setErr(-.12 * factor, .12 * factor, "hard_scale")
        source = directory / "source.yoda"
        campaign._write_yoda_objects(self.yoda, [prediction], source)
        if mode == "no-hard-scale":
            selected = directory / "selected.yoda"
            campaign._write_no_scale_plot_copy(source, selected)
        elif mode == "stat-only":
            selected = directory / "selected.yoda"
            campaign._write_selected_scale_plot_copy(
                {"analysis": {"name": "MC_RESIDUAL_TEST"},
                 "outputs": {"S": {"a1": path}}}, source, selected, 1.,
            )
        else:
            selected = source
        prediction = self.yoda.read(str(selected))[path]
        inputs = directory / "input.yoda"
        campaign._write_yoda_objects(self.yoda, [reference, prediction], inputs)
        config = directory / "panels.plot"
        config.write_text(
            "BEGIN PLOT /MC_RESIDUAL_TEST/.*\n"
            "RatioPlot=1\nRatioPlotYLabel=MC/Data\nLogX=0\nLogY=0\n"
            "Title=Native uncertainty rendering check\nEND PLOT\n",
            encoding="utf-8",
        )
        output = directory / "plots"
        command = [sys.executable, str(ROOT / "scripts/rivet_mkhtml_safe.py"),
                   self.tool, "--dry-run", "--no-rivet-refs", "--config", str(config),
                   "-o", str(output), str(inputs) + ":Title=Synthetic MC"]
        if residual:
            command.append("--deviation")
        result = subprocess.run(command, cwd=directory, env=dict(os.environ),
                                text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                timeout=60)
        self.assertEqual(result.returncode, 0, result.stdout)
        scripts = [p for p in output.rglob("*.py") if not p.name.endswith("__data.py")]
        self.assertEqual(len(scripts), 1, result.stdout)
        script = scripts[0]
        data = runpy.run_path(str(script.with_name(script.stem + "__data.py")))
        mc_labels = [label for label in data["yvals"] if label != "Data"]
        self.assertEqual(len(mc_labels), 1)
        bands = {"monte_carlo": [.03, .06]}
        if mode != "stat-only":
            bands["pdf_68"] = [.04, .08]
        if mode == "all-sources":
            bands.update(scale_down=[-.12, -.24], scale_up=[.12, .24])
        # Exercise the same main-panel overlay as plot_campaign: it disables
        # ordinary main error bars and must leave lower-panel errors visible.
        campaign.add_theory_uncertainty_overlay(script, bands, mc_labels[0])
        campaign.ensure_plot_canvas_draw(script)
        previous = Path.cwd()
        try:
            os.chdir(script.parent)
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                state = runpy.run_path(str(script), run_name="__main__")
        finally:
            os.chdir(previous)
        self.addCleanup(self.plt.close, state["fig"])
        self.assertTrue(script.with_suffix(".png").is_file())
        self.assertTrue(script.with_suffix(".pdf").is_file())
        return state

    def assert_mc_error_segments(self, state, centres, half_widths):
        """Inspect displayed vertical segments, not only precomputed YODA arrays."""
        axis = state["ratio0_ax"]
        containers = [item for item in axis.containers if isinstance(item, self.errorbar_type)]
        labels = list(state["dataf"]["ratio0_yvals"])
        self.assertEqual(len(containers), len(labels))
        mc_labels = [label for label in labels if label != "Data"]
        self.assertEqual(len(mc_labels), 1)
        container = containers[labels.index(mc_labels[0])]
        self.assertTrue(container.has_yerr)
        # ErrorbarContainer stores horizontal collections before vertical ones.
        vertical = container.lines[2][-1]
        segments = self.np.asarray(vertical.get_segments())
        self.assertEqual(segments.shape, (2, 2, 2))
        expected = self.np.column_stack((self.np.asarray(centres) - half_widths,
                                        self.np.asarray(centres) + half_widths))
        self.np.testing.assert_allclose(segments[:, :, 1], expected, rtol=3e-5, atol=3e-5)
        self.np.testing.assert_allclose(segments[:, :, 0], [[.5, .5], [1.5, 1.5]])

    def test_residual_errors_follow_the_selected_uncertainty_sources(self):
        # Total MC errors are .13/.26, .05/.10, or .03/.06. Reference errors
        # are .05/.10, hence the displayed residual half-widths below.
        for mode, error in (("all-sources", 2.6), ("no-hard-scale", 1.), ("stat-only", .6)):
            with self.subTest(mode=mode):
                state = self.render(mode, residual=True)
                self.assert_mc_error_segments(state, [2., 2.], self.np.array([error, error]))
                containers = [item for item in state["ratio0_ax"].containers
                              if isinstance(item, self.errorbar_type)]
                labels = list(state["dataf"]["ratio0_yvals"])
                reference = containers[labels.index("Data")]
                segments = self.np.asarray(reference.lines[2][-1].get_segments())
                self.np.testing.assert_allclose(segments[:, :, 1], 0.)

    def test_true_ratio_keeps_its_existing_error_normalization(self):
        state = self.render("all-sources", residual=False)
        # True ratios divide by the reference values .2/.4, not its errors.
        self.assert_mc_error_segments(state, [1.5, 1.5], self.np.array([.65, .65]))


if __name__ == "__main__":
    unittest.main()
