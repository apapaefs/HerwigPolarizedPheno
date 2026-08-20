#!/usr/bin/env python3
from __future__ import annotations

import contextlib
import hashlib
import io
import json
import math
import re
import sys
import tempfile
import unittest
from pathlib import Path


DISPOL_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(DISPOL_ROOT / "scripts"))
import run_experimental_campaign as campaign  # noqa: E402


PROTON_ID = "COMPASS_2016_I1357198"
DEUTERON_ID = "COMPASS_2017_I1501480"


def snapshot(identifier: str) -> dict:
    return json.loads(
        (DISPOL_ROOT / "data" / "experimental" / identifier / "reference.json")
        .read_text(encoding="utf-8")
    )


def descriptor(identifier: str) -> dict:
    return campaign.get_measurement(identifier)


class CompassReferenceTests(unittest.TestCase):
    def test_compass_2010_cpp_and_reference_use_the_exact_published_edges(self) -> None:
        expected = [
            0.004, 0.005, 0.006, 0.008, 0.010, 0.020, 0.030, 0.040,
            0.060, 0.100, 0.150, 0.200, 0.250, 0.350, 0.500, 0.700,
        ]
        data = snapshot("COMPASS_2010_I843494")
        self.assertEqual(data["bin_edges"], expected)
        source = (
            DISPOL_ROOT / "analyses" / "rivet" / "dis" /
            "COMPASS_2010_I843494.cc"
        ).read_text(encoding="utf-8")
        match = re.search(r"_xEdges\s*=\s*\{([^}]*)\}", source, re.DOTALL)
        self.assertIsNotNone(match)
        cpp_edges = [
            float(token.strip())
            for token in match.group(1).split(",")
            if token.strip()
        ]
        self.assertEqual(cpp_edges, expected)
        self.assertNotIn(0.014, cpp_edges)

    def test_proton_values_errors_edges_and_means(self) -> None:
        data = snapshot(PROTON_ID)
        self.assertEqual(
            data["bin_edges"],
            [0.0025, 0.004, 0.005, 0.006, 0.008, 0.010, 0.014, 0.020,
             0.030, 0.040, 0.060, 0.100, 0.150, 0.200, 0.250, 0.350,
             0.500, 0.700],
        )
        points = data["points"]
        self.assertEqual(
            [point["x_mean"] for point in points],
            [0.0036, 0.0045, 0.0055, 0.007, 0.009, 0.0119, 0.0167,
             0.0244, 0.0346, 0.0488, 0.0768, 0.122, 0.172, 0.223,
             0.292, 0.407, 0.57],
        )
        self.assertEqual(
            [point["q2_mean"] for point in points],
            [1.1, 1.23, 1.39, 1.61, 1.91, 2.33, 3.03, 4.11, 5.6,
             7.64, 11.7, 18, 24.8, 31.3, 39.5, 52, 67.4],
        )
        self.assertEqual(
            [point["value"] for point in points],
            [0.02, 0.017, 0.02, 0.0244, 0.019, 0.0431, 0.0719,
             0.0788, 0.088, 0.114, 0.166, 0.264, 0.318, 0.337,
             0.389, 0.484, 0.73],
        )
        self.assertEqual(
            [point["stat"] for point in points],
            [0.017, 0.012, 0.012, 0.0093, 0.01, 0.0086, 0.0091,
             0.0097, 0.013, 0.013, 0.014, 0.019, 0.027, 0.036,
             0.037, 0.055, 0.11],
        )
        self.assertEqual(
            [point["systematic_combined"] for point in points],
            [0.007, 0.005, 0.005, 0.0041, 0.006, 0.0045, 0.006,
             0.0065, 0.01, 0.009, 0.013, 0.019, 0.024, 0.03,
             0.029, 0.051, 0.09],
        )

    def test_deuteron_values_errors_edges_and_means(self) -> None:
        data = snapshot(DEUTERON_ID)
        self.assertEqual(
            data["bin_edges"],
            [0.004, 0.005, 0.006, 0.008, 0.010, 0.020, 0.030, 0.040,
             0.060, 0.100, 0.150, 0.200, 0.250, 0.350, 0.500, 0.700],
        )
        points = data["points"]
        self.assertEqual(
            [point["x_mean"] for point in points],
            [0.0046, 0.0055, 0.00699, 0.00899, 0.0141, 0.0244, 0.0346,
             0.0487, 0.0766, 0.121, 0.171, 0.222, 0.29, 0.405, 0.567],
        )
        self.assertEqual(
            [point["q2_mean"] for point in points],
            [1.1, 1.22, 1.39, 1.62, 2.19, 3.29, 4.43, 6.06, 9, 13.5,
             18.6, 23.8, 31.1, 43.9, 60.8],
        )
        self.assertEqual(
            [point["value"] for point in points],
            [-0.0054, 0.0003, -0.0011, -0.0087, -0.0011, 0.0075,
             0.0095, 0.0159, 0.0527, 0.095, 0.121, 0.16, 0.19,
             0.317, 0.494],
        )
        self.assertEqual(
            [point["stat"] for point in points],
            [0.0074, 0.0058, 0.0042, 0.0049, 0.0032, 0.0048, 0.0064,
             0.0063, 0.007, 0.01, 0.015, 0.021, 0.023, 0.037, 0.082],
        )
        self.assertEqual(
            [point["systematic_combined"] for point in points],
            [0.0048, 0.0043, 0.0023, 0.0031, 0.0024, 0.0034, 0.0042,
             0.0044, 0.0072, 0.011, 0.016, 0.02, 0.022, 0.036, 0.084],
        )

    def test_raw_hepdata_payloads_and_explicit_proton_override(self) -> None:
        for identifier, expected_hash in (
            (PROTON_ID, "5bf7e5ae4895321dd22ec535fa4dc41df9fb98e0c5bf6b6b05e313d8b2a8a404"),
            (DEUTERON_ID, "6bfd1254974854b6ddc6fafb14ffe2698b7b948fe7a07b3829a62d16b62e3d8d"),
        ):
            measurement = descriptor(identifier)
            raw_path = DISPOL_ROOT / measurement["reference"]["raw_snapshot"]
            payload = raw_path.read_bytes()
            self.assertEqual(hashlib.sha256(payload).hexdigest(), expected_hash)
            rows = campaign.parse_hepdata_table_json(payload, measurement["reference"])
            campaign.validate_reference_snapshot(snapshot(identifier), rows)
            self.assertEqual(len(rows), measurement["reference"]["expected_rows"])
        proton_raw = json.loads(
            (DISPOL_ROOT / descriptor(PROTON_ID)["reference"]["raw_snapshot"])
            .read_text(encoding="utf-8")
        )
        self.assertEqual(float(proton_raw["values"][4]["x"][0]["high"]), 0.001)
        proton_rows = campaign.parse_hepdata_table_json(
            (DISPOL_ROOT / descriptor(PROTON_ID)["reference"]["raw_snapshot"]).read_bytes(),
            descriptor(PROTON_ID)["reference"],
        )
        self.assertEqual(proton_rows[4]["x_high"], 0.010)
        self.assertEqual(proton_rows[1]["x_mean_source"], "bin midpoint")

    def test_g1_is_preserved_but_not_the_reference_observable(self) -> None:
        self.assertEqual(
            [point["g1_value"] for point in snapshot(PROTON_ID)["points"]],
            [0.6, 0.43, 0.44, 0.43, 0.27, 0.51, 0.642, 0.514, 0.424,
             0.401, 0.376, 0.372, 0.298, 0.224, 0.166, 0.095, 0.0396],
        )
        self.assertEqual(
            [point["g1_value"] for point in snapshot(DEUTERON_ID)["points"]],
            [-0.13, 0, -0.016, -0.121, -0.01, 0.043, 0.043, 0.051,
             0.111, 0.123, 0.101, 0.093, 0.0675, 0.05, 0.0203],
        )
        self.assertEqual(snapshot(PROTON_ID)["observable"], "A1p")
        self.assertEqual(snapshot(DEUTERON_ID)["observable"], "A1d")

    def test_r1998_massive_muon_depolarization_and_eta_regression(self) -> None:
        self.assertAlmostEqual(campaign.r1998(0.0036, 1.1), 0.4016700739043282, places=14)
        self.assertAlmostEqual(
            campaign.compass_depolarization(0.0036, 0.8, 1.1),
            0.7995584671341194,
            places=14,
        )
        self.assertAlmostEqual(
            campaign.compass_eta(0.0036, 0.8, 1.1),
            0.002101403250116935,
            places=14,
        )
        self.assertAlmostEqual(campaign.r1998(0.122, 18.0), 0.08798078107921133, places=14)
        self.assertAlmostEqual(
            campaign.compass_depolarization(0.122, 0.4, 18.0),
            0.4368845056116429,
            places=14,
        )
        with self.assertRaises(ValueError):
            campaign.compass_depolarization(0.1, 0.0, 4.0)


class CompassCampaignTests(unittest.TestCase):
    def test_registry_and_nominal_comparison_matrices(self) -> None:
        proton = descriptor(PROTON_ID)
        deuteron = descriptor(DEUTERON_ID)
        self.assertEqual(len(campaign.build_job_matrix(proton, 100, 100, 1, 1)), 8)
        self.assertEqual(len(campaign.build_job_matrix(deuteron, 100, 100, 1, 1)), 16)
        proton_all = campaign.build_job_matrix(
            proton, 100, 100, 1, 1, include_comparisons=True, lo_events=100
        )
        deuteron_all = campaign.build_job_matrix(
            deuteron, 100, 100, 1, 1, include_comparisons=True, lo_events=100
        )
        self.assertEqual(len(proton_all), 30)
        self.assertEqual(len(deuteron_all), 60)
        self.assertEqual({job["component"] for job in deuteron_all}, {"P", "N"})
        self.assertEqual(len({job["id"] for job in deuteron_all}), 60)
        self.assertEqual(len({job["seed"] for job in deuteron_all}), 60)

    def test_deuteron_target_coefficients_and_arithmetic(self) -> None:
        measurement = descriptor(DEUTERON_ID)
        family = campaign.campaign_family_specs(measurement, False)["nominal"]
        self.assertEqual(
            campaign._target_component_coefficients(measurement, family, "unpolarized"),
            {"P": 0.5, "N": 0.5},
        )
        self.assertEqual(
            campaign._target_component_coefficients(measurement, family, "longitudinal"),
            {"P": 0.4625, "N": 0.4625},
        )
        edges = [0.0, 1.0]
        values = {
            "P:PP": 4.0, "P:PM": 2.0, "P:MP": 2.0, "P:MM": 4.0,
            "N:PP": 2.0, "N:PM": 1.0, "N:MP": 1.0, "N:MM": 2.0,
        }
        ordinary = {
            label: campaign.BinSeries(edges, [value], [0.0])
            for label, value in values.items()
        }
        weighted = {
            label: campaign.BinSeries(edges, [2.0 * value], [0.0])
            for label, value in values.items()
        }
        covariance = {
            label: campaign.BinSeries(edges, [math.sqrt(2.0) * value], [0.0])
            for label, value in values.items()
        }
        uu = {
            f"{component}:{helicity}": 0.5 * 0.25
            for component in ("P", "N")
            for helicity in ("PP", "PM", "MP", "MM")
        }
        signs = {"PP": 0.25, "PM": -0.25, "MP": -0.25, "MM": 0.25}
        ll = {
            f"{component}:{helicity}": 0.4625 * signs[helicity]
            for component in ("P", "N")
            for helicity in signs
        }
        parity = {
            helicity: campaign.linear_combine_series(
                {
                    component: ordinary[f"{component}:{helicity}"]
                    for component in ("P", "N")
                },
                {"P": 0.5, "N": 0.5},
            )
            for helicity in ("PP", "PM", "MP", "MM")
        }
        result = campaign._asymmetry_outputs(
            ordinary, weighted, covariance, uu, ll, parity
        )
        self.assertAlmostEqual(result["sigma_uu"].values[0], 2.25)
        self.assertAlmostEqual(result["sigma_ll"].values[0], 0.69375)
        self.assertAlmostEqual(result["apar_values"][0], 0.69375 / 2.25)
        self.assertAlmostEqual(result["a1_values"][0], 2.0 * 0.69375 / 2.25)
        self.assertAlmostEqual(result["parity_pp_mm_values"][0], 0.0)
        self.assertAlmostEqual(result["parity_pm_mp_values"][0], 0.0)

    def test_card_families_and_target_invariants(self) -> None:
        proton_dir = DISPOL_ROOT / "cards" / "experimental" / PROTON_ID
        deuteron_dir = DISPOL_ROOT / "cards" / "experimental" / DEUTERON_ID
        self.assertEqual(len(list(proton_dir.glob(f"{PROTON_ID}_[PM][PM]-*.in"))), 8)
        self.assertEqual(len(list((proton_dir / "comparisons").glob("*.in"))), 22)
        self.assertEqual(len(list(deuteron_dir.glob(f"{DEUTERON_ID}_[PN]_[PM][PM]-*.in"))), 16)
        self.assertEqual(len(list((deuteron_dir / "comparisons").glob("*.in"))), 44)
        for identifier, energy, q2max, w2min in (
            (PROTON_ID, "200.0*GeV", "190.0*GeV2", "12.0*GeV2"),
            (DEUTERON_ID, "160.0*GeV", "100.0*GeV2", "16.0*GeV2"),
        ):
            common = (
                DISPOL_ROOT / "cards" / "experimental" / identifier /
                f"{identifier}-Common.in"
            ).read_text(encoding="utf-8")
            for required in (
                "EventHandler:BeamA /Herwig/Particles/mu+",
                f"FixedTargetLuminosity:BeamEMaxA {energy}",
                f"NeutralCurrentCut:MaxQ2 {q2max}",
                f"NeutralCurrentCut:MinW2 {w2min}",
                "NNPDF40_nlo_pch_as_01180",
                "NNPDFpol20_nlo_as_01180",
                "PowhegMEDISNCPol:GammaZ Gamma",
                "UsePOWHEGRealSpinVertex Yes",
                "ShowerHandler:Interactions QCD",
                "ShowerHandler:MPIHandler NULL",
            ):
                self.assertIn(required, common)
            self.assertNotIn("HadronizationHandler  NULL", common)
            self.assertNotIn("CascadeHandler NULL", common)
            self.assertNotIn("DecayHandler NULL", common)
        neutron = (
            deuteron_dir / f"{DEUTERON_ID}_N_PP-POSNLO.in"
        ).read_text(encoding="utf-8")
        self.assertIn("EventHandler:BeamB /Herwig/Particles/n0", neutron)
        self.assertIn("FixedTargetLuminosity:TargetParticle /Herwig/Particles/n0", neutron)

    def test_analysis_sources_lock_fiducial_and_model_choices(self) -> None:
        helper = (
            DISPOL_ROOT / "analyses" / "rivet" / "dis" / "COMPASSInclusiveDIS.hh"
        ).read_text(encoding="utf-8")
        self.assertIn("lepton.pid() != -13", helper)
        self.assertIn("hadron.pid() != 2212 && hadron.pid() != 2112", helper)
        self.assertIn("Rivet/Projections/PromptFinalState.hh", helper)
        self.assertIn("0.0485, 0.5470, 2.0621", helper)
        self.assertIn("muonMass = 0.1056583755", helper)
        proton = (
            DISPOL_ROOT / "analyses" / "rivet" / "dis" / f"{PROTON_ID}.cc"
        ).read_text(encoding="utf-8")
        deuteron = (
            DISPOL_ROOT / "analyses" / "rivet" / "dis" / f"{DEUTERON_ID}.cc"
        ).read_text(encoding="utf-8")
        proton_info = (
            DISPOL_ROOT / "analyses" / "rivet" / "dis" / f"{PROTON_ID}.info"
        ).read_text(encoding="utf-8")
        deuteron_info = (
            DISPOL_ROOT / "analyses" / "rivet" / "dis" / f"{DEUTERON_ID}.info"
        ).read_text(encoding="utf-8")
        legacy_proton = (
            DISPOL_ROOT / "analyses" / "rivet" / "dis" /
            "COMPASS_2010_I843494.cc"
        ).read_text(encoding="utf-8")
        for identifier in (PROTON_ID, DEUTERON_ID):
            plot = (
                DISPOL_ROOT / "analyses" / "rivet" / "dis" /
                f"{identifier}.plot"
            ).read_text(encoding="utf-8")
            self.assertIn(r"\mathrm{GeV}", plot)
            self.assertIn(r"\theta_\mu", plot)
            self.assertNotIn("\t", plot)
        self.assertIn("dis.Q2 > 190.0", proton)
        self.assertIn("dis.W2 <= 12.0", proton)
        self.assertIn("dis.Q2 > 100.0", deuteron)
        self.assertIn("dis.W2 <= 16.0", deuteron)
        for source in (proton, deuteron, legacy_proton):
            self.assertIn(
                'PromptFinalState(Cuts::pid == -13), "PromptMuons"',
                source,
            )
            self.assertIn(
                'apply<PromptFinalState>(event, "PromptMuons")',
                source,
            )
        self.assertIn("Beams: [[ANTIMUON, PROTON]]", proton_info)
        self.assertIn(
            "Beams: [[ANTIMUON, PROTON], [ANTIMUON, NEUTRON]]",
            deuteron_info,
        )

        self.assertIn(
            "_xEdges = {0.004, 0.005, 0.006, 0.008, 0.010, 0.020, 0.030,",
            legacy_proton,
        )
        self.assertIn(
            "0.040, 0.060, 0.100, 0.150, 0.200, 0.250, 0.350,",
            legacy_proton,
        )
        self.assertIn("0.500, 0.700};", legacy_proton)
        self.assertNotIn("0.010, 0.014, 0.020", legacy_proton)

    def test_dry_run_counts_are_resolved_and_non_mutating(self) -> None:
        for identifier, expected_nominal, expected_all in (
            (PROTON_ID, 8, 30),
            (DEUTERON_ID, 16, 60),
        ):
            for comparisons, expected in ((False, expected_nominal), (True, expected_all)):
                argv = [
                    "prepare", "--measurement", identifier,
                    "--tag", "compass-unit-dry-run", "--smoke", "--dry-run",
                ]
                if comparisons:
                    argv.append("--comparisons")
                output = io.StringIO()
                with contextlib.redirect_stdout(output):
                    self.assertEqual(campaign.main(argv), 0)
                plan = json.loads(output.getvalue())
                self.assertEqual(plan["logical_jobs"], expected)
                self.assertEqual(plan["shard_jobs"], expected)
                self.assertEqual(sum(job["events"] for job in plan["jobs"]), 100 * expected)

    def test_incomplete_deuteron_component_is_rejected(self) -> None:
        measurement = descriptor(DEUTERON_ID)
        jobs = campaign.build_job_matrix(measurement, 100, 100, 1, 3000)
        manifest = {"configuration": {"shards": 1}, "jobs": jobs[:-1]}
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for job in jobs[:-1]:
                output = root / job["output_yoda"]
                output.parent.mkdir(parents=True, exist_ok=True)
                output.write_text("non-empty", encoding="utf-8")
            with self.assertRaisesRegex(
                campaign.CampaignError, "missing family/helicity/order components"
            ):
                campaign._require_complete_matrix(manifest, measurement, root)


if __name__ == "__main__":
    unittest.main()
