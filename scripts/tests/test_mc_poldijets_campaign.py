#!/usr/bin/env python3
"""Focused contracts for the internal MC_POLDIJETS measurement."""

from __future__ import annotations

import argparse
import contextlib
import io
import math
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import run_experimental_campaign as experimental  # noqa: E402
import run_phenomenology_campaign as campaign  # noqa: E402


def arguments(**updates: object) -> argparse.Namespace:
    values: dict[str, object] = {
        "profile": "central",
        "jobs": 8,
        "shards": 1,
        "seed_base": None,
        "lo_events": None,
        "posnlo_events": None,
        "negnlo_events": None,
        "smoke": False,
        "polarized_pdf_members": None,
        "unpolarized_pdf_members": None,
        "scales": None,
        "families": None,
        "jet_kt_min_gev": None,
        "nominal_prediction": None,
        "plot_comparisons": False,
    }
    values.update(updates)
    return argparse.Namespace(**values)


def series(
    values: list[float], variances: list[float]
) -> experimental.BinSeries:
    return experimental.BinSeries(
        [float(index) for index in range(len(values) + 1)],
        values,
        variances,
    )


class MCPOLDIJETSCampaignTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.measurement = campaign.discover_pp_registry()["MC_POLDIJETS"]

    def test_internal_definition_matches_descriptor(self) -> None:
        snapshot = campaign._measurement_snapshot(self.measurement)
        self.assertEqual(snapshot["kind"], "internal_observable_definition")
        self.assertEqual(snapshot["sqrt_s_gev"], 510.0)
        configured = self.measurement["channels"]["dijets"]["raw_objects"]
        self.assertEqual(set(snapshot["observables"]), set(configured))
        self.assertEqual(len(configured), 15)
        for observable, raw_object in configured.items():
            self.assertEqual(
                snapshot["observables"][observable]["raw_object"], raw_object
            )
        self.assertNotIn("reference_yoda", self.measurement["analysis"])

    def test_spin_on_off_matrix_and_card_delta(self) -> None:
        options = campaign._resolved_options(
            arguments(families="nominal,shower_spin_off"), self.measurement
        )
        jobs = campaign.build_job_matrix(self.measurement, options)
        self.assertEqual(len(jobs), 8)
        self.assertEqual(len({job["seed"] for job in jobs}), 8)
        self.assertEqual(options["jet_kt_min_gev"], 3.0)
        self.assertEqual(
            {job["helicity"] for job in jobs}, {"PP", "PM", "MP", "MM"}
        )
        expected_instance = (
            "MC_POLDIJETS:ETAMAX=1.5:LEVEL=HARDPARTON:PTJ1MIN=5.0:"
            "PTJ2MIN=4.0:PTJ3MIN=2.0:R=0.5"
        )
        self.assertEqual(
            {job["analysis_instance"] for job in jobs}, {expected_instance}
        )

        nominal = next(
            job for job in jobs
            if job["family"] == "nominal" and job["helicity"] == "PP"
        )
        control = next(
            job for job in jobs
            if job["family"] == "shower_spin_off" and job["helicity"] == "PP"
        )
        nominal_card = campaign._card_text(self.measurement, nominal)
        control_card = campaign._card_text(self.measurement, control)
        self.assertIn("ShowerHandler:SpinCorrelations Yes", nominal_card)
        self.assertIn("ShowerHandler:SpinCorrelations No", control_card)
        self.assertIn("JetKtCut:MinKT 3.0*GeV", nominal_card)
        common = (
            ROOT
            / self.measurement["cards"]["directory"]
            / self.measurement["cards"]["common"]
        ).read_text(encoding="utf-8")
        self.assertIn("Luminosity:Energy 510.0*GeV", common)
        normalized_nominal = nominal_card.replace(
            "set /Herwig/Shower/ShowerHandler:SpinCorrelations Yes\n", ""
        )
        normalized_control = control_card.replace(
            "set /Herwig/Shower/ShowerHandler:SpinCorrelations No\n", ""
        ).replace(control["stem"], nominal["stem"])
        self.assertEqual(normalized_nominal, normalized_control)

    def test_cross_section_delta_sigma_and_asymmetry_arithmetic(self) -> None:
        samples = {
            "PP": series([12.0], [0.4]),
            "PM": series([8.0], [0.3]),
            "MP": series([8.0], [0.2]),
            "MM": series([12.0], [0.5]),
        }
        prediction = campaign._mc_poldijets_prediction({"jet1_pt": samples})
        self.assertAlmostEqual(
            prediction["SigmaUU_jet1_pt"]["values"][0], 10.0
        )
        self.assertAlmostEqual(
            prediction["DeltaSigmaLL_jet1_pt"]["values"][0], 2.0
        )
        self.assertAlmostEqual(
            prediction["ALL_jet1_pt"]["values"][0], 0.2
        )
        expected_error = math.sqrt((0.4 + 0.3 + 0.2 + 0.5) / 16.0)
        self.assertAlmostEqual(
            prediction["SigmaUU_jet1_pt"]["errors"][0], expected_error
        )
        self.assertAlmostEqual(
            prediction["DeltaSigmaLL_jet1_pt"]["errors"][0], expected_error
        )
        self.assertAlmostEqual(
            prediction["SingleSpinA_jet1_pt"]["values"][0], 0.0
        )
        self.assertAlmostEqual(
            prediction["SingleSpinB_jet1_pt"]["values"][0], 0.0
        )
        self.assertAlmostEqual(
            prediction["Parity_PP_MM_jet1_pt"]["values"][0], 0.0
        )

    def test_primary_paths_have_no_external_points(self) -> None:
        snapshot = campaign._measurement_snapshot(self.measurement)
        for prefix in ("ALL_", "DeltaSigmaLL_", "SigmaUU_"):
            observable = prefix + "jet1_pt"
            self.assertEqual(
                campaign._reference_path(
                    self.measurement, observable, snapshot
                ),
                f"/MC_POLDIJETS/{observable}",
            )
            self.assertTrue(
                campaign._primary_reference_observable(
                    self.measurement, snapshot, observable
                )
            )
            self.assertIsNone(
                campaign._reference_points(
                    self.measurement, snapshot, observable
                )
            )
        self.assertIsNone(
            campaign._reference_path(
                self.measurement, "SingleSpinA_jet1_pt", snapshot
            )
        )

    def test_loose_source_and_blue_control_contracts(self) -> None:
        source = (ROOT / self.measurement["analysis"]["source"]).read_text(
            encoding="utf-8"
        )
        for raw_object in self.measurement["channels"]["dijets"][
            "raw_objects"
        ].values():
            self.assertIn(f'"{raw_object}"', source)
        self.assertIn("_leadingPtMin = getOption<double>(\"PTJ1MIN\", 5.0)", source)
        self.assertIn("_subleadingPtMin = getOption<double>(\"PTJ2MIN\", 4.0)", source)
        self.assertNotIn("dphi <", source)
        self.assertNotIn("deta <", source)
        control = self.measurement["families"]["shower_spin_off"]
        self.assertEqual(control["plot_options"]["LineColor"], "#0077BB")
        self.assertEqual(control["shower_spin_correlations"], "off")

    def test_fetch_data_validates_without_external_download(self) -> None:
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            status = campaign.main(
                ["fetch-data", "--measurement", "MC_POLDIJETS"]
            )
        self.assertEqual(status, 0)
        self.assertIn("no external data or reference YODA", output.getvalue())

    def test_plotting_falls_back_to_mathtext_without_dvipng(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            style = output / "default.mplstyle"
            style.write_text(
                "font.size: 11\ntext.usetex: True\n",
                encoding="utf-8",
            )
            rendering = campaign._configure_plot_text_rendering(
                output, {"PATH": str(output / "missing-tools")}
            )
            self.assertEqual(
                rendering,
                {
                    "mode": "mathtext",
                    "missing_tools": ["latex", "dvipng"],
                },
            )
            self.assertIn(
                "text.usetex: False", style.read_text(encoding="utf-8")
            )


if __name__ == "__main__":
    unittest.main()
