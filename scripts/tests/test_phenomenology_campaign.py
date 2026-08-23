#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import runpy
import sys
import tempfile
import unittest
from pathlib import Path


DISPOL_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(DISPOL_ROOT / "scripts"))
import run_experimental_campaign as experimental  # noqa: E402
import run_phenomenology_campaign as campaign  # noqa: E402


def args(profile: str = "central", **updates: object) -> argparse.Namespace:
    values = {
        "profile": profile, "jobs": 4, "shards": 1, "seed_base": None,
        "lo_events": 100, "smoke": False, "polarized_pdf_members": None,
        "unpolarized_pdf_members": None, "scales": None, "families": None,
        "posnlo_events": None, "negnlo_events": None,
    }
    values.update(updates)
    return argparse.Namespace(**values)


def series(values: list[float], variance: float = 0.0) -> experimental.BinSeries:
    return experimental.BinSeries(
        [float(index) for index in range(len(values) + 1)], values,
        [variance for _ in values],
    )


class RegistryAndMatrixTests(unittest.TestCase):
    def test_combined_registry_and_scientific_labels(self) -> None:
        registry = campaign.discover_all()
        for identifier in (
            "HERMES_2007_I726689", "HERMES_2007_I726689_LEGACY",
            "COMPASS_2010_I843494", "COMPASS_2016_I1357198",
            "COMPASS_2017_I1501480", "STAR_2019_I1708793",
            "PHENIX_2023_I2033856", "STAR_2021_I1850855",
            "STAR_2022_I1949588", "HERMES_2019_I1698889",
        ):
            self.assertIn(identifier, registry)
        for identifier in ("STAR_2019_I1708793", "PHENIX_2023_I2033856"):
            families = registry[identifier]["families"]
            self.assertEqual(families["nominal"]["label"], "LO+PS (polarized; MPI on)")
            self.assertEqual(families["mpi_off"]["label"], "LO+PS (polarized; MPI off)")
            self.assertEqual(families["unpolarized_closure"]["label"],
                             "LO+PS (unpolarized closure; MPI on)")

    def test_central_and_paper_job_matrices(self) -> None:
        registry = campaign.discover_pp_registry()
        star = registry["STAR_2019_I1708793"]
        phenix = registry["PHENIX_2023_I2033856"]
        star_central = campaign.build_job_matrix(star, campaign._resolved_options(args(), star))
        phenix_central = campaign.build_job_matrix(phenix, campaign._resolved_options(args(), phenix))
        self.assertEqual(len(star_central), 12)
        self.assertEqual(len(phenix_central), 4)
        star_paper = campaign.build_job_matrix(star, campaign._resolved_options(args("paper"), star))
        phenix_paper = campaign.build_job_matrix(phenix, campaign._resolved_options(args("paper"), phenix))
        self.assertEqual(len(star_paper), 2436)
        self.assertEqual(len(phenix_paper), 812)
        for jobs in (star_paper, phenix_paper):
            self.assertEqual(len({job["id"] for job in jobs}), len(jobs))
            self.assertEqual(len({job["seed"] for job in jobs}), len(jobs))
        self.assertTrue(all(job["order"] == "LO" for job in star_paper + phenix_paper))

    def test_variation_selectors_do_not_form_replica_cartesian_product(self) -> None:
        selected = args("paper", polarized_pdf_members="0,7",
                        unpolarized_pdf_members="0,9", scales="0.5,1,2")
        self.assertEqual(campaign._variation_points(selected),
                         [(0,0,1.0),(0,9,1.0),(7,0,1.0),(0,0,0.5),(0,0,2.0)])
        measurement = campaign.discover_pp_registry()["STAR_2019_I1708793"]
        selected.families = "all"
        options = campaign._resolved_options(selected, measurement)
        jobs = campaign.build_job_matrix(measurement, options)
        # nominal: 3 channels x 5 points x 4 helicities; MPI-off and closure:
        # 3 x 4 plus 3 x 1 central-only jobs.
        self.assertEqual(len(jobs), 75)

    def test_generated_scale_mpi_pdf_and_unpolarized_overrides(self) -> None:
        measurement = campaign.discover_pp_registry()["STAR_2019_I1708793"]
        job = campaign.build_job_matrix(
            measurement,
            campaign._resolved_options(args(
                "paper", polarized_pdf_members="7",
                unpolarized_pdf_members="central", scales="2",
                families="all",
            ), measurement),
        )
        scale_job = next(item for item in job if item["family"] == "nominal" and
                         item["channel"] == "Wplus" and item["scale"] == 2.0)
        scale_card = campaign._card_text(measurement, scale_job)
        self.assertIn("set /Herwig/MatrixElements/MEqq2W2ff:ScalePreFactor 4", scale_card)
        replica_job = next(item for item in job if item["polarized_pdf_member"] == 7)
        self.assertIn("set /Herwig/Partons/STARPolarizedPDF:Member 7",
                      campaign._card_text(measurement, replica_job))
        mpi_job = next(item for item in job if item["family"] == "mpi_off")
        self.assertIn("PowhegShowerHandler:MPIHandler NULL",
                      campaign._card_text(measurement, mpi_job))
        closure_job = next(item for item in job if item["family"] == "unpolarized_closure")
        closure = campaign._card_text(measurement, closure_job)
        self.assertIn("FirstLongitudinalPolarization 0", closure)
        self.assertIn("SecondLongitudinalPolarization 0", closure)

    def test_data_total_default_and_optional_statistical_overlay(self) -> None:
        measurement = campaign.discover_pp_registry()["STAR_2019_I1708793"]
        snapshot = campaign.validate_vendored("STAR_2019_I1708793")
        points = campaign._pp_reference_overlay_points(
            measurement, snapshot, "Wplus_AL"
        )
        self.assertEqual([point["plot_x"] for point in points],
                         [-1.22, -0.71, -0.24, 0.25, 0.72, 1.22])
        with tempfile.TemporaryDirectory() as temporary:
            source = "fig, ax = None, None\n# set plot metadata as defined above\n"
            default_script = Path(temporary) / "default.py"
            default_script.write_text(source, encoding="utf-8")
            experimental.add_experimental_error_overlay(default_script, points)
            generated_default = default_script.read_text(encoding="utf-8")
            self.assertNotIn("_data_stat", generated_default)
            self.assertNotIn("inner bars: statistical", generated_default)

            component_script = Path(temporary) / "components.py"
            component_script.write_text(source, encoding="utf-8")
            experimental.add_experimental_error_overlay(
                component_script, points, show_statistical=True
            )
            generated = component_script.read_text(encoding="utf-8")
            self.assertIn("inner bars: statistical; outer bars: total", generated)
            self.assertIn("_data_x = [-1.22, -0.71", generated)
            self.assertIn("elinewidth=1.7", generated)

            low_q2_script = Path(temporary) / "low-q2.py"
            low_q2_script.write_text(
                source,
                encoding="utf-8",
            )
            low_q2_point = dict(points[0], nonperturbative_extrapolation=True)
            experimental.add_experimental_error_overlay(
                low_q2_script, [low_q2_point]
            )
            generated_low_q2 = low_q2_script.read_text(encoding="utf-8")
            self.assertIn("ax.scatter(_lowq_x", generated_low_q2)
            self.assertNotIn("_data_stat", generated_low_q2)

    def test_named_theory_uncertainty_overlay(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            script = Path(temporary) / "plot.py"
            script.write_text(
                "import numpy as np\n"
                "class FakeAxis:\n"
                "    def __init__(self):\n"
                "        self.bands = []\n"
                "    def fill_between(self, x, low, high, **kwargs):\n"
                "        self.bands.append((list(x), list(low), list(high), kwargs))\n"
                "        return object()\n"
                "ax = FakeAxis()\n"
                "dataf = {\n"
                "    'yvals': {'Data': [1.0, 2.0], 'prediction.yoda': [0.4, 0.5]},\n"
                "    'xedges': {'Data': [0.0, 1.0, 2.0], 'prediction.yoda': [0.0, 1.0, 2.0]},\n"
                "    'xpoints': {'Data': [0.5, 1.5], 'prediction.yoda': [0.5, 1.5]},\n"
                "}\n"
                "styles = {\n"
                "    'Data': {'yerrorbars': 1, 'color': 'black'},\n"
                "    'prediction.yoda': {'yerrorbars': 1, 'color': 'red'},\n"
                "}\n"
                "legend_handles = {}\n"
                "labels = {'legend': [[], '']}\n"
                "# curve from input yoda files in main panel\n"
                "assert np.allclose(ax.bands[-1][0], [0.0, 1.0, 2.0])\n"
                "assert np.allclose(ax.bands[-1][1], [0.39, 0.39, 0.48])\n"
                "assert np.allclose(ax.bands[-1][2], [0.41, 0.41, 0.52])\n"
                "assert ax.bands[-1][3]['step'] == 'pre'\n"
                "assert ax.bands[-1][3]['color'] == 'red'\n"
                "assert styles['Data']['yerrorbars'] == 1\n"
                "assert styles['prediction.yoda']['yerrorbars'] == 0\n"
                "legend_items = list(legend_handles.values())\n",
                encoding="utf-8",
            )
            experimental.add_theory_uncertainty_overlay(script, {
                "monte_carlo": [0.01, 0.02],
                "pdf_68": [0.03, 0.04],
                "scale_down": [-0.05, -0.06],
                "scale_up": [0.07, 0.08],
            }, "prediction.yoda")
            generated = script.read_text(encoding="utf-8")
            self.assertIn("_theory_label = 'prediction.yoda'", generated)
            self.assertNotIn("next(iter(dataf['yvals']))", generated)
            self.assertIn("styles[_theory_label]['yerrorbars'] = 0", generated)
            self.assertIn("_theory_mc_handle = ax.fill_between", generated)
            self.assertNotIn("_theory_mc_handle = ax.errorbar", generated)
            self.assertIn("PDF 68%", generated)
            self.assertIn("hard scale", generated)
            self.assertIn("_theory_mc = [0.01, 0.02]", generated)
            runpy.run_path(str(script))

    def test_fixed_target_comparison_plotting_is_explicitly_opt_in(self) -> None:
        default = campaign.make_parser().parse_args([
            "plot", "--measurement", "HERMES_2007_I726689",
            "--tag", "plot-default",
        ])
        self.assertNotIn("--plot-comparisons", campaign._legacy_arguments(default))

        comparisons = campaign.make_parser().parse_args([
            "plot", "--measurement", "HERMES_2007_I726689",
            "--tag", "plot-comparisons", "--plot-comparisons",
        ])
        self.assertIn(
            "--plot-comparisons", campaign._legacy_arguments(comparisons)
        )

        full = campaign.make_parser().parse_args([
            "full", "--measurement", "HERMES_2007_I726689",
            "--tag", "full-comparisons", "--plot-comparisons",
            "--plot-data-components",
        ])
        self.assertIn("--plot-comparisons", campaign._legacy_arguments(full))
        self.assertIn("--plot-data-components", campaign._legacy_arguments(full))

        data_components = campaign.make_parser().parse_args([
            "plot", "--measurement", "HERMES_2007_I726689",
            "--tag", "plot-data-components", "--plot-data-components",
        ])
        self.assertIn(
            "--plot-data-components", campaign._legacy_arguments(data_components)
        )
        experimental_components = experimental.make_parser().parse_args([
            "plot", "--measurement", "HERMES_2007_I726689",
            "--tag", "plot-data-components", "--plot-data-components",
        ])
        self.assertTrue(experimental_components.plot_data_components)

        refresh = campaign.make_parser().parse_args([
            "plot", "--measurement", "HERMES_2007_I726689",
            "--tag", "plot-label-refresh",
            "--allow-plot-metadata-refresh",
        ])
        self.assertIn(
            "--allow-plot-metadata-refresh",
            campaign._legacy_arguments(refresh),
        )
        experimental_refresh = experimental.make_parser().parse_args([
            "plot", "--measurement", "HERMES_2007_I726689",
            "--tag", "plot-label-refresh",
            "--allow-plot-metadata-refresh",
        ])
        self.assertTrue(experimental_refresh.allow_plot_metadata_refresh)


class ObservableArithmeticTests(unittest.TestCase):
    def test_star_beam_symmetrization_eta_folding_and_exchange_closure(self) -> None:
        samples = {
            "PP": series([12.0] * 6, 0.4), "PM": series([8.0] * 6, 0.3),
            "MP": series([9.0] * 6, 0.2), "MM": series([11.0] * 6, 0.5),
        }
        prediction = campaign._star_prediction("Wplus", samples)
        self.assertEqual(prediction["Wplus_AL"]["values"], [0.025] * 6)
        self.assertEqual(prediction["Wplus_ALL"]["values"], [0.15] * 3)
        self.assertTrue(all(math.isfinite(value)
                            for value in prediction["Wplus_AL"]["errors"]))
        closures = campaign._beam_exchange_closures(samples)
        self.assertEqual(closures["PP"]["values"], [0.0] * 6)
        self.assertEqual(closures["MM"]["values"], [0.0] * 6)

    def test_phenix_cross_section_average_all_and_parity_closures(self) -> None:
        physical = {
            "PP": series([12.0], 0.4), "PM": series([8.0], 0.3),
            "MP": series([8.0], 0.2), "MM": series([12.0], 0.5),
        }
        output = campaign._phenix_prediction({
            "inclusive_cross_section": physical,
            "isolated_cross_section": physical,
            "isolated_all": physical,
        })
        self.assertAlmostEqual(output["inclusive_cross_section"]["values"][0], 10.0)
        self.assertAlmostEqual(output["isolated_cross_section"]["values"][0], 10.0)
        self.assertAlmostEqual(output["isolated_all"]["values"][0], 0.2)
        self.assertAlmostEqual(output["Parity_PP_MM_isolated_all"]["values"][0], 0.0)
        self.assertAlmostEqual(output["Parity_PM_MP_isolated_all"]["values"][0], 0.0)
        self.assertAlmostEqual(output["SingleSpinA_isolated_all"]["values"][0], 0.0)
        self.assertAlmostEqual(output["SingleSpinB_isolated_all"]["values"][0], 0.0)

    def test_covariance_zero_denominator_and_pulls(self) -> None:
        self.assertEqual(experimental.ratio_with_covariance(1,1,0,1,0), (None, None))
        value, error = experimental.ratio_with_covariance(2,0.4,4,0.6,0.1)
        self.assertAlmostEqual(value, 0.5)
        self.assertTrue(math.isfinite(error))
        pulls, chi2, count = campaign._pulls(
            {"values": [0.2, None], "errors": [0.1, None]},
            [{"value": 0.1, "stat": 0.2, "systematic_combined": 0.1},
             {"value": 0.0, "stat": 0.1, "systematic_combined": 0.0}],
        )
        self.assertIsNotNone(pulls[0]); self.assertIsNone(pulls[1])
        self.assertEqual(count, 1); self.assertGreater(chi2, 0.0)

    def test_pdf_replica_quadrature_scale_envelope_and_mpi_shift(self) -> None:
        measurement = campaign.discover_pp_registry()["PHENIX_2023_I2033856"]
        central_key = ("nominal", "photon", 0, 0, 1.0, "on")
        predictions = {
            central_key: {"isolated_all": {"values": [1.0], "errors": [0.1]}},
            ("nominal", "photon", 1, 0, 1.0, "on"):
                {"isolated_all": {"values": [1.2], "errors": [0.1]}},
            ("nominal", "photon", 2, 0, 1.0, "on"):
                {"isolated_all": {"values": [0.8], "errors": [0.1]}},
            ("nominal", "photon", 0, 1, 1.0, "on"):
                {"isolated_all": {"values": [1.1], "errors": [0.1]}},
            ("nominal", "photon", 0, 2, 1.0, "on"):
                {"isolated_all": {"values": [0.9], "errors": [0.1]}},
            ("nominal", "photon", 0, 0, 0.5, "on"):
                {"isolated_all": {"values": [0.7], "errors": [0.1]}},
            ("nominal", "photon", 0, 0, 2.0, "on"):
                {"isolated_all": {"values": [1.4], "errors": [0.1]}},
            ("mpi_off", "photon", 0, 0, 1.0, "off"):
                {"isolated_all": {"values": [0.95], "errors": [0.1]}},
        }
        # The aggregator expects full replica index lists, but deliberately uses
        # whichever members are available and requires at least two per axis.
        bands = campaign.aggregate_uncertainties(predictions, measurement)["isolated_all"]
        self.assertAlmostEqual(bands["polarized_pdf_68"][0], math.sqrt(0.08))
        self.assertAlmostEqual(bands["unpolarized_pdf_68"][0], math.sqrt(0.02))
        self.assertAlmostEqual(bands["pdf_68"][0], math.sqrt(0.10))
        self.assertAlmostEqual(bands["scale_down"][0], -0.3)
        self.assertAlmostEqual(bands["scale_up"][0], 0.4)
        self.assertAlmostEqual(bands["mpi_shift"][0], -0.05)


class SourceAndAnalysisContractTests(unittest.TestCase):
    def test_phenix_exact_binning_and_invariant_cross_section_formula(self) -> None:
        source = (DISPOL_ROOT / "analyses/rivet/pp/PHENIX_2023_I2033856.cc").read_text()
        self.assertIn("10.0,12.0,14.0,16.0", source)
        self.assertIn("22.0,24.0,26.0,28.0,30.0", source)
        self.assertIn("_inclusive->fill(pt, 1.0/pt)", source)
        self.assertIn("sf/(2.0*M_PI*0.5)", source)
        self.assertIn("particle.pT()/GeV < 0.2", source)
        self.assertIn("particle.E()/GeV < 0.3", source)

    def test_fully_masked_plot_is_suppressed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "plot.py"
            path.write_text("# generated plotting driver\n")
            path.with_name("plot__data.py").write_text(
                "yvals = {'prediction': [float('nan'), float('nan')]}\n"
            )
            self.assertFalse(experimental.plot_script_has_finite_y(path))
            path.with_name("plot__data.py").write_text(
                "yvals = {'prediction': [0.0, 1.0]}\n"
            )
            self.assertTrue(experimental.plot_script_has_finite_y(path))
            path.write_text(
                "# generated plotting driver\n"
                "xLims = (0.0, 2.0)\n"
                "yLims = (nan, nan)\n"
            )
            self.assertFalse(experimental.plot_script_has_finite_y(path))
            path.write_text(
                "# generated plotting driver\n"
                "xLims = (0.0, 2.0)\n"
                "yLims = (-1.0, 1.0)\n"
            )
            self.assertTrue(experimental.plot_script_has_finite_y(path))


if __name__ == "__main__":
    unittest.main()
