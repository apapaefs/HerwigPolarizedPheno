#!/usr/bin/env python3
"""Regression tests for the local Rivet HTML compatibility wrapper."""

from __future__ import annotations

import contextlib
import io
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

import numpy as np


DISPOL_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(DISPOL_ROOT / "scripts"))

import rivet_mkhtml_safe as safe  # noqa: E402


class RivetMkhtmlSafeTests(unittest.TestCase):
    def test_large_numpy_lists_are_not_serialized_with_ellipsis(self) -> None:
        yoda_module = types.ModuleType("yoda")
        yoda_module.__path__ = []  # type: ignore[attr-defined]
        plotting_module = types.ModuleType("yoda.plotting")
        plotting_module.__path__ = []  # type: ignore[attr-defined]
        fetch_module = types.ModuleType("yoda.plotting.fetch_data")
        rendered: list[str] = []

        def write_lists(_cmd, _value_type, values):
            rendered.append(str(values["Data"]))

        fetch_module.writeLists = write_lists  # type: ignore[attr-defined]
        yoda_module.plotting = plotting_module  # type: ignore[attr-defined]
        plotting_module.fetch_data = fetch_module  # type: ignore[attr-defined]
        modules = {
            "yoda": yoda_module,
            "yoda.plotting": plotting_module,
            "yoda.plotting.fetch_data": fetch_module,
        }
        with mock.patch.dict(sys.modules, modules):
            safe._patch_yoda_write_lists()
            fetch_module.writeLists(  # type: ignore[attr-defined]
                None, "xpoints", {"Data": np.arange(2332, dtype=float)}
            )
        self.assertEqual(len(rendered), 1)
        self.assertNotIn("...", rendered[0])
        self.assertIn("2.331e+03", rendered[0])

    def test_nested_yoda_path_parent_is_created(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            parent = safe._ensure_plot_output_parent(
                output,
                "/HERMES_2019_I1698889/DIAGNOSTICS/UnpolarizedCosPhi_test",
            )
            self.assertEqual(
                Path(parent),
                output / "HERMES_2019_I1698889" / "DIAGNOSTICS",
            )
            self.assertTrue(Path(parent).is_dir())
            self.assertEqual(
                safe._nested_plot_style_path(
                    "/HERMES_2019_I1698889/DIAGNOSTICS/"
                    "UnpolarizedCosPhi_test"
                ),
                "../../",
            )

    def test_plot_path_cannot_escape_output_directory(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(ValueError, "escapes"):
                safe._ensure_plot_output_parent(temporary, "/../outside")

    def test_script_generator_failures_are_retained_for_nonzero_exit(self) -> None:
        yoda_module = types.ModuleType("yoda")
        yoda_module.__path__ = []  # type: ignore[attr-defined]
        plotting_module = types.ModuleType("yoda.plotting")
        plotting_module.__path__ = []  # type: ignore[attr-defined]
        generator_module = types.ModuleType("yoda.plotting.script_generator")

        def fail_process(yaml_file, _plot_name, outdir, *_args, **_kwargs):
            expected = (
                Path(outdir)
                / "HERMES_2019_I1698889"
                / "DIAGNOSTICS"
            )
            self.assertTrue(expected.is_dir())
            self.assertEqual(yaml_file["stylepath"], "../../")
            raise RuntimeError("synthetic plotting failure")

        generator_module.process = fail_process  # type: ignore[attr-defined]
        yoda_module.plotting = plotting_module  # type: ignore[attr-defined]
        plotting_module.script_generator = generator_module  # type: ignore[attr-defined]

        modules = {
            "yoda": yoda_module,
            "yoda.plotting": plotting_module,
            "yoda.plotting.script_generator": generator_module,
        }
        with tempfile.TemporaryDirectory() as temporary:
            with mock.patch.dict(sys.modules, modules):
                failures = safe._patch_script_generator_traceback()
                with contextlib.redirect_stderr(io.StringIO()):
                    with self.assertRaisesRegex(
                        RuntimeError, "synthetic plotting failure"
                    ):
                        generator_module.process(
                            {"histograms": {}},
                            "/HERMES_2019_I1698889/DIAGNOSTICS/test",
                            temporary,
                        )
                self.assertEqual(len(failures), 1)
                self.assertIn(
                    "/HERMES_2019_I1698889/DIAGNOSTICS/test",
                    failures[0],
                )
                self.assertIn("synthetic plotting failure", failures[0])

    def test_swallowed_generator_failure_makes_wrapper_exit_nonzero(self) -> None:
        yoda_module = types.ModuleType("yoda")
        yoda_module.__path__ = []  # type: ignore[attr-defined]
        plotting_module = types.ModuleType("yoda.plotting")
        plotting_module.__path__ = []  # type: ignore[attr-defined]
        generator_module = types.ModuleType("yoda.plotting.script_generator")

        def fail_process(*_args, **_kwargs):
            raise RuntimeError("failure swallowed by fake rivet-mkhtml")

        generator_module.process = fail_process  # type: ignore[attr-defined]
        yoda_module.plotting = plotting_module  # type: ignore[attr-defined]
        plotting_module.script_generator = generator_module  # type: ignore[attr-defined]
        modules = {
            "yoda": yoda_module,
            "yoda.plotting": plotting_module,
            "yoda.plotting.script_generator": generator_module,
        }

        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "plots"
            tool = Path(temporary) / "fake-rivet-mkhtml.py"
            tool.write_text(
                "from yoda.plotting import script_generator\n"
                "try:\n"
                "    script_generator.process(\n"
                "        {'histograms': {}},\n"
                "        '/ANALYSIS/DIAGNOSTICS/test',\n"
                f"        {str(output)!r},\n"
                "    )\n"
                "except RuntimeError:\n"
                "    pass\n",
                encoding="utf-8",
            )
            with (
                mock.patch.dict(sys.modules, modules),
                mock.patch.object(safe, "_patch_yoda_parseyaml"),
                mock.patch.object(safe, "_patch_yoda_write_lists"),
                mock.patch.object(safe, "_patch_fetch_data_band_annotations"),
                mock.patch.object(sys, "argv", list(sys.argv)),
                contextlib.redirect_stdout(io.StringIO()),
                contextlib.redirect_stderr(io.StringIO()) as stderr,
            ):
                self.assertEqual(safe.main([str(tool)]), 1)
            self.assertIn(
                "rivet-mkhtml suppressed one or more plot-generation failures",
                stderr.getvalue(),
            )


if __name__ == "__main__":
    unittest.main()
