#!/usr/bin/env python3
from __future__ import annotations

import contextlib
import hashlib
import io
import json
import math
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest import mock


DISPOL_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(DISPOL_ROOT / "scripts"))
import run_experimental_campaign as campaign  # noqa: E402


MEASUREMENT_ID = "HERMES_2007_I726689"
SNAPSHOT_PATH = DISPOL_ROOT / "data" / "experimental" / MEASUREMENT_ID / "reference.json"
CARD_DIR = DISPOL_ROOT / "cards" / "experimental" / MEASUREMENT_ID
COMPARISON_CARD_DIR = CARD_DIR / "comparisons"
NOMINAL_COMPATIBILITY_SIGNATURE = "33d36d153b40601a6bb3992c199a3d31618c6ef2dc78ea08598cacbb11f72a86"


class HermesReferenceTests(unittest.TestCase):
    def setUp(self) -> None:
        self.reference = json.loads(SNAPSHOT_PATH.read_text(encoding="utf-8"))
        self.snapshot = self.reference["datasets"][0]

    def test_all_published_values_and_errors(self) -> None:
        self.assertEqual(
            self.snapshot["bin_edges"],
            [0.0212, 0.0295, 0.0362, 0.0444, 0.0568, 0.0727, 0.0929, 0.119,
             0.152, 0.194, 0.249, 0.318, 0.406, 0.52, 0.665, 0.9],
        )
        points = self.snapshot["points"]
        self.assertEqual(
            [point["x_mean"] for point in points],
            [0.0264, 0.0329, 0.0403, 0.0506, 0.0648, 0.0829, 0.1059, 0.1354,
             0.173, 0.2209, 0.282, 0.3598, 0.4583, 0.5819, 0.7248],
        )
        self.assertEqual(
            [point["q2_mean"] for point in points],
            [1.12, 1.25, 1.38, 1.54, 2.01, 2.45, 2.97, 3.59, 4.31, 5.13,
             6.11, 7.24, 8.53, 10.16, 12.21],
        )
        self.assertEqual(
            [point["value"] for point in points],
            [0.1113, 0.0975, 0.0897, 0.1061, 0.109, 0.1997, 0.2011, 0.2681,
             0.3246, 0.3361, 0.4094, 0.5169, 0.6573, 0.6647, 1.1976],
        )
        self.assertEqual(
            [point["stat"] for point in points],
            [0.0351, 0.0268, 0.0236, 0.0201, 0.0202, 0.0206, 0.0218, 0.0236,
             0.0263, 0.0305, 0.037, 0.0486, 0.0714, 0.1302, 0.2763],
        )
        self.assertEqual(
            [point["systematic_combined"] for point in points],
            [0.0095, 0.0087, 0.0082, 0.01, 0.0093, 0.0147, 0.015, 0.019,
             0.0223, 0.0218, 0.0263, 0.0333, 0.0416, 0.0469, 0.0827],
        )
        for point in points:
            quadrature = math.sqrt(sum(value * value for value in point["systematics"].values()))
            self.assertAlmostEqual(quadrature, point["systematic_combined"], delta=1.0e-4)

    def test_official_table_parser_and_archive_checksum(self) -> None:
        lines = [
            "* HERMES test table",
            *(
                f" {point['x_mean']:.4f} {point['q2_mean']:.4f} {point['value']:.4f} "
                f"{point['stat']:.4f} {point['systematic_combined']:.4f}"
                for point in self.snapshot["points"]
            ),
        ]
        member = ("\n".join(lines) + "\n").encode("ascii")
        buffer = io.BytesIO()
        with tarfile.open(fileobj=buffer, mode="w") as archive:
            info = tarfile.TarInfo("A1p_a15")
            info.size = len(member)
            archive.addfile(info, io.BytesIO(member))
        rows = campaign.extract_tar_reference(
            buffer.getvalue(), "A1p_a15", hashlib.sha256(member).hexdigest()
        )
        campaign.validate_reference_snapshot(self.snapshot, rows)
        self.assertEqual(len(rows), 15)
        with self.assertRaises(campaign.CampaignError):
            campaign.extract_tar_reference(buffer.getvalue(), "A1p_a15", "0" * 64)

    def test_r1990_and_depolarization_regression(self) -> None:
        self.assertAlmostEqual(campaign.r1990(0.0264, 1.12), 0.34166533177946749213, places=14)
        self.assertAlmostEqual(campaign.depolarization(0.0264, 0.7, 1.12), 0.7030502026628281, places=14)
        self.assertAlmostEqual(campaign.r1990(0.173, 4.31), 0.18625742296198792093, places=14)
        self.assertAlmostEqual(campaign.depolarization(0.173, 0.4, 4.31), 0.405734444688617, places=14)
        with self.assertRaises(ValueError):
            campaign.depolarization(0.1, 0.0, 4.0)

    def test_generated_reference_yoda(self) -> None:
        try:
            import yoda
        except ImportError:
            self.skipTest("YODA Python bindings are not active")
        reference = DISPOL_ROOT / "analyses" / "rivet" / "dis" / f"{MEASUREMENT_ID}.yoda.gz"
        objects = yoda.read(str(reference))
        obj = objects[f"/REF/{MEASUREMENT_ID}/d14-x01-y01"]
        self.assertEqual(list(obj.xEdges()), self.snapshot["bin_edges"])
        self.assertEqual(obj.numBins(), 15)
        self.assertEqual(
            obj.annotation("PublishedXMeans"),
            [point["x_mean"] for point in self.snapshot["points"]],
        )
        self.assertEqual(
            obj.annotation("PublishedQ2MeansGeV2"),
            [point["q2_mean"] for point in self.snapshot["points"]],
        )
        for index, point in enumerate(self.snapshot["points"], start=1):
            bin_object = obj.bin(index)
            self.assertAlmostEqual(bin_object.val(), point["value"])
            self.assertAlmostEqual(bin_object.errAvg("stat"), point["stat"])
            for source, error in point["systematics"].items():
                self.assertAlmostEqual(bin_object.errAvg(source), error)


class CampaignArithmeticTests(unittest.TestCase):
    def test_normalized_shards_and_signed_orders(self) -> None:
        first = campaign.BinSeries([0.0, 1.0], [10.0], [4.0])
        second = campaign.BinSeries([0.0, 1.0], [14.0], [9.0])
        combined = campaign.combine_shard_series([first, second], [100, 300])
        self.assertAlmostEqual(combined.values[0], 13.0)
        self.assertAlmostEqual(combined.variances[0], (0.25**2) * 4.0 + (0.75**2) * 9.0)
        signed = campaign.add_independent_series(
            [combined, campaign.BinSeries([0.0, 1.0], [-2.0], [1.0])]
        )
        self.assertAlmostEqual(signed.values[0], 11.0)
        self.assertAlmostEqual(signed.variances[0], combined.variances[0] + 1.0)

    def test_helicity_signs_and_covariance_propagation(self) -> None:
        ordinary = {
            "PP": campaign.BinSeries([0.0, 1.0], [12.0], [0.4]),
            "PM": campaign.BinSeries([0.0, 1.0], [8.0], [0.3]),
            "MP": campaign.BinSeries([0.0, 1.0], [6.0], [0.2]),
            "MM": campaign.BinSeries([0.0, 1.0], [10.0], [0.5]),
        }
        weighted = {
            key: campaign.BinSeries([0.0, 1.0], [1.5 * value.values[0]], [2.25 * value.variances[0]])
            for key, value in ordinary.items()
        }
        covariance = {
            key: campaign.BinSeries([0.0, 1.0], [0.0], [1.5 * value.variances[0]])
            for key, value in ordinary.items()
        }
        uu = {key: 0.25 for key in ordinary}
        ll = {"PP": 0.25, "PM": -0.25, "MP": -0.25, "MM": 0.25}
        result = campaign._asymmetry_outputs(ordinary, weighted, covariance, uu, ll)
        self.assertAlmostEqual(result["sigma_uu"].values[0], 9.0)
        self.assertAlmostEqual(result["sigma_ll"].values[0], 2.0)
        self.assertAlmostEqual(result["apar_values"][0], 2.0 / 9.0)
        self.assertAlmostEqual(result["a1_values"][0], 1.5 * 2.0 / 9.0)
        self.assertIsNotNone(result["a1_errors"][0])

    def test_ratio_empty_bin_and_parity_residual(self) -> None:
        self.assertEqual(campaign.ratio_with_covariance(1.0, 1.0, 0.0, 1.0, 0.0), (None, None))
        value, error = campaign.ratio_with_covariance(2.0, 0.4, 4.0, 0.6, 0.1)
        expected_variance = 0.4 / 16.0 + 4.0 * 0.6 / 256.0 - 4.0 * 0.1 / 64.0
        self.assertAlmostEqual(value, 0.5)
        self.assertAlmostEqual(error, math.sqrt(expected_variance))
        parity, parity_error = campaign.parity_residual(12.0, 0.4, 10.0, 0.5)
        self.assertAlmostEqual(parity, 2.0 / 22.0)
        expected = (20.0 / 22.0**2) ** 2 * 0.4 + (-24.0 / 22.0**2) ** 2 * 0.5
        self.assertAlmostEqual(parity_error, math.sqrt(expected))


class CampaignRegistryAndCardTests(unittest.TestCase):
    def setUp(self) -> None:
        self.measurement = campaign.get_measurement(MEASUREMENT_ID)

    def test_registry_and_sixteen_job_matrix(self) -> None:
        registry = campaign.discover_registry()
        self.assertIn(MEASUREMENT_ID, registry)
        jobs = campaign.build_job_matrix(self.measurement, 100000, 10000, 1, 726689)
        self.assertEqual(len(jobs), 16)
        self.assertEqual({job["helicity"] for job in jobs}, {"PP", "PM", "MP", "MM"})
        self.assertEqual({job["order"] for job in jobs}, {"POSNLO", "NEGNLO"})
        self.assertEqual(sum(job["events"] for job in jobs), 880000)
        self.assertEqual(len({job["seed"] for job in jobs}), 16)
        sharded = campaign.build_job_matrix(self.measurement, 101, 11, 2, 9000)
        self.assertEqual(len(sharded), 32)
        for component in ("P", "N"):
            for helicity in ("PP", "PM", "MP", "MM"):
                self.assertEqual(sum(job["events"] for job in sharded if job["component"] == component and job["helicity"] == helicity and job["order"] == "POSNLO"), 101)
                self.assertEqual(sum(job["events"] for job in sharded if job["component"] == component and job["helicity"] == helicity and job["order"] == "NEGNLO"), 11)

    def test_expanded_comparison_matrix(self) -> None:
        jobs = campaign.build_job_matrix(
            self.measurement,
            300000,
            30000,
            10,
            1726689,
            include_comparisons=True,
            lo_events=300000,
        )
        logical = {
            (job["family"], job["component"], job["helicity"], job["order"])
            for job in jobs
        }
        self.assertEqual(len(logical), 60)
        self.assertEqual(len(jobs), 600)
        self.assertEqual(sum(job["events"] for job in jobs), 10_980_000)
        self.assertEqual(len({job["seed"] for job in jobs}), 600)
        by_family = {
            family: [job for job in jobs if job["family"] == family]
            for family in {job["family"] for job in jobs}
        }
        self.assertEqual(len(by_family["nominal"]), 160)
        self.assertEqual(len(by_family["unpolarized_nlo"]), 40)
        self.assertEqual(len(by_family["polarized_lo"]), 80)
        self.assertEqual(len(by_family["no_real_spin_nlo"]), 160)
        self.assertEqual(len(by_family["no_shower_spin_nlo"]), 160)
        self.assertEqual(
            {(job["helicity"], job["order"]) for job in by_family["unpolarized_nlo"]},
            {("00", "POSNLO"), ("00", "NEGNLO")},
        )
        self.assertEqual({job["order"] for job in by_family["polarized_lo"]}, {"LO"})
        self.assertTrue(
            all(job["card_input"].startswith("comparisons/") for job in jobs if job["family"] != "nominal")
        )

    def test_paper_profile_pdf_scale_matrix_and_generated_card(self) -> None:
        points = campaign._variation_points(
            "paper", "0,7", "0,9", "0.5,1,2"
        )
        self.assertEqual(
            points,
            [(0, 0, 1.0), (0, 9, 1.0), (7, 0, 1.0),
             (0, 0, 0.5), (0, 0, 2.0)],
        )
        jobs = campaign.build_job_matrix(
            self.measurement,
            100,
            10,
            1,
            1726689,
            include_comparisons=True,
            lo_events=100,
            variation_points=points,
        )
        # Five nominal points x sixteen target/helicity/sign jobs, plus
        # 44 central-only comparison jobs.
        self.assertEqual(len(jobs), 124)
        self.assertEqual(len({job["id"] for job in jobs}), 124)
        self.assertEqual(len({job["seed"] for job in jobs}), 124)
        scale_job = next(
            job for job in jobs
            if job["family"] == "nominal" and job["scale"] == 2.0
            and job["helicity"] == "PP" and job["order"] == "POSNLO"
        )
        source = CARD_DIR / f"{scale_job['base_stem']}.in"
        card = campaign._dis_variation_card_text(
            source, CARD_DIR / self.measurement["cards"]["common"], scale_job
        )
        self.assertIn("HERMESPDF:Member 0", card)
        self.assertIn("HERMESDiffPDF:Member 0", card)
        self.assertIn("PowhegMEDISNCPol:ScaleFactor 2", card)
        self.assertIn(f"saverun {scale_job['stem']} EventGenerator", card)
        replica_job = next(job for job in jobs if job["polarized_pdf_member"] == 7)
        replica_card = campaign._dis_variation_card_text(
            CARD_DIR / f"{replica_job['base_stem']}.in",
            CARD_DIR / self.measurement["cards"]["common"],
            replica_job,
        )
        self.assertIn("HERMESDiffPDF:Member 7", replica_card)
        manifest = {
            "configuration": {
                "shards": 1, "comparisons": True,
                "variation_points": [list(point) for point in points],
            },
            "jobs": jobs,
        }
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for job in jobs:
                job["status"] = "success"
                output = root / job["output_yoda"]
                output.parent.mkdir(parents=True, exist_ok=True)
                output.write_text("non-empty", encoding="utf-8")
            groups = campaign._require_complete_matrix(
                manifest, self.measurement, root
            )
            self.assertEqual(len(groups), 124)

    def test_dis_pdf_replica_quadrature_and_scale_envelope(self) -> None:
        def result(a1: float, apar: float) -> dict[str, dict[str, object]]:
            return {
                "Q2GT1": {
                    "a1_values": [a1], "a1_errors": [0.03],
                    "apar_values": [apar], "apar_errors": [0.02],
                }
            }

        predictions = {
            (0, 0, 1.0): result(1.0, 0.5),
            (1, 0, 1.0): result(1.1, 0.55),
            (2, 0, 1.0): result(0.9, 0.45),
            (0, 1, 1.0): result(1.2, 0.60),
            (0, 2, 1.0): result(0.8, 0.40),
            (0, 0, 0.5): result(0.7, 0.35),
            (0, 0, 2.0): result(1.4, 0.70),
        }
        bands = campaign.aggregate_dis_uncertainties(predictions)
        a1 = bands["Q2GT1"]["a1"]
        self.assertAlmostEqual(a1["polarized_pdf_68"][0], math.sqrt(0.02))
        self.assertAlmostEqual(a1["unpolarized_pdf_68"][0], math.sqrt(0.08))
        self.assertAlmostEqual(a1["pdf_68"][0], math.sqrt(0.10))
        self.assertAlmostEqual(a1["scale_down"][0], -0.3)
        self.assertAlmostEqual(a1["scale_up"][0], 0.4)
        self.assertEqual(a1["monte_carlo"], [0.03])

    def test_dis_pdf_replica_variance_is_about_replica_mean(self) -> None:
        def result(a1: float) -> dict[str, dict[str, object]]:
            return {
                "Q2GT1": {
                    "a1_values": [a1], "a1_errors": [0.03],
                    "apar_values": [0.5], "apar_errors": [0.02],
                }
            }

        predictions = {
            (0, 0, 1.0): result(10.0),
            (1, 0, 1.0): result(1.0),
            (2, 0, 1.0): result(2.0),
            (3, 0, 1.0): result(4.0),
        }
        bands = campaign.aggregate_dis_uncertainties(predictions)
        expected = math.sqrt(
            ((1.0 - 7.0 / 3.0) ** 2 +
             (2.0 - 7.0 / 3.0) ** 2 +
             (4.0 - 7.0 / 3.0) ** 2) / 2.0
        )
        self.assertAlmostEqual(
            bands["Q2GT1"]["a1"]["polarized_pdf_68"][0], expected
        )

    def test_dual_target_requires_fresh_signature_and_comparisons_are_separate(self) -> None:
        self.assertNotEqual(
            campaign.measurement_signature(self.measurement),
            NOMINAL_COMPATIBILITY_SIGNATURE,
        )
        comparison = campaign.comparison_signature(self.measurement)
        self.assertEqual(len(comparison), 64)
        nominal_options = {
            "jobs": 4,
            "shards": 10,
            "seed_base": 1726689,
            "posnlo_events": 300000,
            "negnlo_events": 30000,
            "lo_events": 300000,
            "smoke": False,
            "comparisons": False,
        }
        nominal_configuration = campaign._manifest_configuration(
            self.measurement,
            "hermes_prod_300k_20260721",
            nominal_options,
            campaign.measurement_signature(self.measurement),
        )
        self.assertNotIn("comparisons", nominal_configuration)
        self.assertNotIn("lo_events", nominal_configuration)
        comparison_options = {**nominal_options, "comparisons": True}
        comparison_configuration = campaign._manifest_configuration(
            self.measurement,
            "hermes_comparisons",
            comparison_options,
            campaign.measurement_signature(self.measurement),
        )
        self.assertTrue(comparison_configuration["comparisons"])
        self.assertEqual(comparison_configuration["lo_events"], 300000)
        self.assertEqual(comparison_configuration["comparison_signature"], comparison)

    def test_postprocessing_refuses_stale_measurement_signatures(self) -> None:
        current = {
            "configuration": {
                "measurement_signature": campaign.measurement_signature(
                    self.measurement
                )
            }
        }
        campaign._assert_manifest_signatures_current(current, self.measurement)
        stale = {
            "configuration": {"measurement_signature": "stale-analysis"}
        }
        with self.assertRaisesRegex(
            campaign.CampaignError, "intentionally refused"
        ):
            campaign._assert_manifest_signatures_current(
                stale, self.measurement
            )

    def test_plot_metadata_refresh_proves_only_the_plot_file_changed(self) -> None:
        historical_plot = b"historical plot metadata\n"
        recorded_signature = "generation-signature"
        manifest = {
            "status": "complete",
            "configuration": {
                "measurement_signature": recorded_signature,
            },
            "jobs": [{"status": "success"}],
            "runtime": {
                "provenance": {
                    "source_control": {"commit": "a" * 40}
                }
            },
        }
        signature = lambda payload: (
            recorded_signature if payload == historical_plot else "unexpected"
        )
        with mock.patch.object(
            campaign,
            "_git_file_at_commit",
            return_value=historical_plot,
        ):
            provenance = campaign.authorize_plot_metadata_refresh(
                manifest,
                self.measurement,
                current_signature="current-signature",
                signature_with_plot_bytes=signature,
            )
        self.assertEqual(
            provenance["mode"], "presentation_only_plot_metadata_refresh"
        )
        self.assertEqual(
            provenance["generation_measurement_signature"],
            recorded_signature,
        )
        self.assertEqual(
            provenance["current_measurement_signature"], "current-signature"
        )

        with mock.patch.object(
            campaign,
            "_git_file_at_commit",
            return_value=historical_plot,
        ):
            with self.assertRaisesRegex(
                campaign.CampaignError, "more than the Rivet .plot metadata"
            ):
                campaign.authorize_plot_metadata_refresh(
                    manifest,
                    self.measurement,
                    current_signature="current-signature",
                    signature_with_plot_bytes=lambda _payload: "other-change",
                )

    def test_card_family_invariants(self) -> None:
        common = (CARD_DIR / f"{MEASUREMENT_ID}-Common.in").read_text(encoding="utf-8")
        for required in (
            "FixedTargetLuminosity:BeamEMaxA 27.6*GeV",
            "NNPDF40_nlo_pch_as_01180",
            "NNPDFpol20_nlo_as_01180",
            "PowhegMEDISNCPol:GammaZ Gamma",
            "UsePOWHEGRealSpinVertex Yes",
            "ShowerHandler:SpinCorrelations Yes",
            "ShowerHandler:Interactions QCD",
            "ShowerHandler:MPIHandler NULL",
            "UseNativeDISWindowGeneration Yes",
            f"Rivet:Analyses 0 {MEASUREMENT_ID}",
        ):
            self.assertIn(required, common)
        self.assertNotIn("HadronizationHandler  NULL", common)
        self.assertNotIn("CascadeHandler NULL", common)
        self.assertNotIn("DecayHandler NULL", common)
        self.assertEqual(len(list(CARD_DIR.glob(f"{MEASUREMENT_ID}_[PN]_[PM][PM]-*.in"))), 16)
        self.assertIn("/Herwig/Particles/n0:PDF /Herwig/Partons/HERMESPDF", common)
        for component in ("P", "N"):
            target = "p+" if component == "P" else "n0"
            for helicity, signs in self.measurement["cards"]["helicities"].items():
                for order, contribution in self.measurement["cards"]["orders"].items():
                    text = (CARD_DIR / f"{MEASUREMENT_ID}_{component}_{helicity}-{order}.in").read_text(encoding="utf-8")
                    self.assertIn(f"FirstLongitudinalPolarization {signs[0]}", text)
                    self.assertIn(f"SecondLongitudinalPolarization {signs[1]}", text)
                    self.assertIn(f"Contribution {contribution}", text)
                    self.assertIn(f"BeamB /Herwig/Particles/{target}", text)
                    self.assertIn(f"TargetParticle /Herwig/Particles/{target}", text)

    def test_analysis_uses_prompt_scattered_positron(self) -> None:
        source = (
            DISPOL_ROOT / "analyses" / "rivet" / "dis" /
            "HERMES_2007_I726689.cc"
        ).read_text(encoding="utf-8")
        self.assertIn(
            'PromptFinalState(Cuts::pid == -11), "PromptPositrons"',
            source,
        )
        self.assertIn(
            'apply<PromptFinalState>(event, "PromptPositrons")',
            source,
        )

    def test_comparison_card_family_invariants(self) -> None:
        cards = sorted(COMPARISON_CARD_DIR.glob("*.in"))
        self.assertEqual(len(cards), 44)
        families = campaign.campaign_family_specs(self.measurement, include_comparisons=True)
        self.assertEqual(
            {family_id: family["label"] for family_id, family in families.items()},
            {
                "nominal": "NLO+PS (polarized; full spin)",
                "unpolarized_nlo": "NLO+PS (unpolarized beams)",
                "polarized_lo": "LO+PS (polarized)",
                "no_real_spin_nlo": "NLO+PS (polarized; Born spin only)",
                "no_shower_spin_nlo": "NLO+PS (polarized; shower spin off)",
            },
        )
        for family_id, family in families.items():
            annotations = campaign._family_annotations(
                self.measurement, family_id, family
            )
            self.assertNotEqual(annotations["ShowerSpinCorrelations"], "unspecified")
            if family_id == "nominal":
                continue
            for helicity, signs in family["helicities"].items():
                for order, contribution in family["orders"].items():
                    stem = family["stem_pattern"].format(component="P", helicity=helicity, order=order)
                    text = (COMPARISON_CARD_DIR / f"{stem}.in").read_text(encoding="utf-8")
                    self.assertIn("read ../HERMES_2007_I726689-Common.in", text)
                    self.assertIn(f"FirstLongitudinalPolarization {signs[0]}", text)
                    self.assertIn(f"SecondLongitudinalPolarization {signs[1]}", text)
                    self.assertIn(f"saverun {stem} EventGenerator", text)
                    if order in {"POSNLO", "NEGNLO"}:
                        self.assertIn(f"Contribution {contribution}", text)
        for path in COMPARISON_CARD_DIR.glob("*NOREALSPIN*.in"):
            text = path.read_text(encoding="utf-8")
            self.assertIn("UsePOWHEGRealSpinVertex No", text)
            self.assertIn("ShowerHandler:SpinCorrelations Yes", text)
        for path in COMPARISON_CARD_DIR.glob("*NOSHOWERSPIN*.in"):
            text = path.read_text(encoding="utf-8")
            self.assertIn("UsePOWHEGRealSpinVertex No", text)
            self.assertIn("ShowerHandler:SpinCorrelations No", text)
        for path in COMPARISON_CARD_DIR.glob("*_LO_*-LO.in"):
            text = path.read_text(encoding="utf-8")
            self.assertIn("erase /Herwig/MatrixElements/SubProcess:MatrixElements 0", text)
            self.assertIn("MatrixElements[0] /Herwig/MatrixElements/MEDISNCPol", text)
            self.assertIn("MEDISNCPol:GammaZ Gamma", text)
            self.assertIn("ShowerHandler:HardEmission None", text)
            self.assertNotIn(":Contribution", text)
        for path in COMPARISON_CARD_DIR.glob("*UNPOL*.in"):
            text = path.read_text(encoding="utf-8")
            self.assertIn("FirstLongitudinalPolarization 0", text)
            self.assertIn("SecondLongitudinalPolarization 0", text)

    def test_dry_run_is_resolved_and_non_mutating(self) -> None:
        tag = "unit-dry-run-do-not-create"
        destination = campaign.campaign_directory(MEASUREMENT_ID, tag)
        before = destination.exists()
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            status = campaign.main([
                "prepare", "--measurement", MEASUREMENT_ID, "--tag", tag, "--smoke", "--dry-run"
            ])
        self.assertEqual(status, 0)
        plan = json.loads(output.getvalue())
        self.assertEqual(plan["logical_jobs"], 16)
        self.assertEqual(plan["shard_jobs"], 16)
        self.assertEqual(sum(job["events"] for job in plan["jobs"]), 1600)
        self.assertEqual(destination.exists(), before)

    def test_comparison_dry_run_is_resolved_and_non_mutating(self) -> None:
        tag = "unit-comparison-dry-run-do-not-create"
        destination = campaign.campaign_directory(MEASUREMENT_ID, tag)
        before = destination.exists()
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            status = campaign.main(
                [
                    "prepare",
                    "--measurement",
                    MEASUREMENT_ID,
                    "--tag",
                    tag,
                    "--smoke",
                    "--comparisons",
                    "--dry-run",
                ]
            )
        self.assertEqual(status, 0)
        plan = json.loads(output.getvalue())
        self.assertEqual(plan["logical_jobs"], 60)
        self.assertEqual(plan["shard_jobs"], 60)
        self.assertEqual(sum(job["events"] for job in plan["jobs"]), 6000)
        self.assertEqual(destination.exists(), before)

    def test_progress_cli_matches_validation_runner_defaults(self) -> None:
        parser = campaign.make_parser()
        args = parser.parse_args(
            ["campaign", "--measurement", MEASUREMENT_ID, "--tag", "tracker-options"]
        )
        self.assertEqual(args.progress_interval, 5.0)
        self.assertEqual(args.max_listed, 12)
        custom = parser.parse_args(
            [
                "full",
                "--measurement",
                MEASUREMENT_ID,
                "--tag",
                "tracker-options-custom",
                "--progress-interval",
                "2.5",
                "--max-listed",
                "3",
                "--comparisons",
                "--lo-events",
                "321",
            ]
        )
        self.assertEqual(custom.progress_interval, 2.5)
        self.assertEqual(custom.max_listed, 3)
        self.assertTrue(custom.comparisons)
        self.assertEqual(custom.lo_events, 321)

    def test_live_progress_markers_payload_and_monitor_files(self) -> None:
        jobs = campaign.build_job_matrix(self.measurement, 100, 100, 1, 3000)
        jobs[0]["status"] = "success"
        jobs[1]["status"] = "failed"
        jobs[2]["status"] = "queued"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            log = root / "active.log"
            log.write_text(
                "event> init 100\nevent> 25/100\revent>     50    100\n",
                encoding="utf-8",
            )
            self.assertEqual(campaign.latest_progress_marker(log), "50/100 (50.0%)")
            activity = campaign.JobActivity(
                job=jobs[2],
                tag="tracker-active-a00",
                log_path=log,
                started_at=0.0,
            )
            herwig = campaign.build_herwig_progress_payload(
                jobs,
                jobs[3:],
                [activity],
                max_listed=12,
                started_at=0.0,
            )
            self.assertEqual(
                {key: herwig[key] for key in ("completed", "running", "pending", "failed", "total")},
                {"completed": 1, "running": 1, "pending": 13, "failed": 1, "total": 16},
            )
            self.assertEqual(herwig["active_rows"][0][2], "50/100 (50.0%)")
            payload = campaign.build_campaign_monitor_payload(
                measurement=self.measurement,
                tag="tracker-test",
                phase="running-herwig",
                started_at=0.0,
                herwig=herwig,
            )
            campaign.write_campaign_monitor_files(root, payload)
            status_json = json.loads(campaign.campaign_status_json_path(root).read_text())
            status_text = campaign.campaign_status_txt_path(root).read_text()
            self.assertEqual(status_json["phase"], "running-herwig")
            self.assertEqual(status_json["herwig"]["running"], 1)
            self.assertIn("Shards: completed 1/16 | running 1 | pending 13 | failed 1", status_text)
            self.assertIn("Active shards:", status_text)
            self.assertIn("Progress", status_text)

    def test_campaign_writes_final_tracker_state(self) -> None:
        options = {"jobs": 2, "shards": 1, "seed_base": 5000, "posnlo_events": 100,
                   "negnlo_events": 100, "smoke": True}
        jobs = campaign.build_job_matrix(self.measurement, 100, 100, 1, 5000)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            root.mkdir(parents=True, exist_ok=True)
            manifest = {
                "manifest_version": 1,
                "measurement": MEASUREMENT_ID,
                "tag": "tracker-integration",
                "created_at": campaign.utc_now(),
                "configuration": campaign._manifest_configuration(
                    self.measurement,
                    "tracker-integration",
                    options,
                    campaign.measurement_signature(self.measurement),
                ),
                "runtime": {"test_runtime": True},
                "jobs": jobs,
                "history": [],
                "status": "prepared",
            }
            campaign.atomic_write_json(root / campaign.MANIFEST_NAME, manifest)

            def fake_job(job, measurement, campaign_dir, runtime, campaign_tag, active, active_lock):
                output = campaign_dir / job["output_yoda"]
                output.parent.mkdir(parents=True, exist_ok=True)
                output.write_text("fake-yoda", encoding="utf-8")
                return campaign.JobResult(str(job["id"]), True, 0, output.stat().st_size, "ok")

            args = campaign.make_parser().parse_args(
                [
                    "campaign",
                    "--measurement",
                    MEASUREMENT_ID,
                    "--tag",
                    "tracker-integration",
                    "--smoke",
                    "--jobs",
                    "2",
                    "--seed-base",
                    "5000",
                    "--progress-interval",
                    "-1",
                ]
            )
            args._tracker_started_at = 0.0
            with mock.patch.object(campaign, "campaign_directory", return_value=root), mock.patch.object(
                campaign, "_run_tracked_campaign_job", side_effect=fake_job
            ), contextlib.redirect_stdout(io.StringIO()):
                campaign.run_campaign(args, self.measurement)
            final_status = json.loads(campaign.campaign_status_json_path(root).read_text())
            final_manifest = json.loads((root / campaign.MANIFEST_NAME).read_text())
            self.assertEqual(final_status["phase"], "campaign-complete")
            self.assertEqual(final_status["herwig"]["completed"], 16)
            self.assertEqual(final_status["herwig"]["pending"], 0)
            self.assertEqual(final_manifest["status"], "complete")

    def test_resume_and_fresh_seed_recovery(self) -> None:
        jobs = campaign.build_job_matrix(self.measurement, 100, 100, 1, 1000)
        jobs[0]["status"] = "success"
        jobs[1]["status"] = "failed"
        manifest = {"jobs": jobs}
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            successful = root / jobs[0]["output_yoda"]
            successful.parent.mkdir(parents=True)
            successful.write_text("non-empty", encoding="utf-8")
            scheduled, blocked = campaign.pending_jobs(manifest, root, recover_failed=False)
            self.assertEqual(len(blocked), 1)
            self.assertEqual(len(scheduled), 14)
            failed_seed = jobs[1]["seed"]
            scheduled, blocked = campaign.pending_jobs(manifest, root, recover_failed=True)
            self.assertFalse(blocked)
            recovered = next(job for job in scheduled if job["id"] == jobs[1]["id"])
            self.assertNotEqual(recovered["seed"], failed_seed)
            self.assertEqual(recovered["attempt"], 1)

    def test_postprocessing_refuses_an_incomplete_component(self) -> None:
        jobs = campaign.build_job_matrix(self.measurement, 100, 100, 1, 2000)
        manifest = {"configuration": {"shards": 1}, "jobs": jobs}
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for job in jobs[:-1]:
                output = root / job["output_yoda"]
                output.parent.mkdir(parents=True, exist_ok=True)
                output.write_text("non-empty", encoding="utf-8")
            with self.assertRaisesRegex(campaign.CampaignError, "incomplete components"):
                campaign._require_complete_matrix(manifest, self.measurement, root)

    def test_expanded_postprocessing_requires_every_comparison_component(self) -> None:
        jobs = campaign.build_job_matrix(
            self.measurement,
            100,
            100,
            1,
            2100,
            include_comparisons=True,
            lo_events=100,
        )
        manifest = {
            "configuration": {"shards": 1, "comparisons": True},
            "jobs": jobs[:-1],
        }
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for job in manifest["jobs"]:
                output = root / job["output_yoda"]
                output.parent.mkdir(parents=True, exist_ok=True)
                output.write_text("non-empty", encoding="utf-8")
            with self.assertRaisesRegex(
                campaign.CampaignError,
                "missing family/helicity/order components",
            ):
                campaign._require_complete_matrix(manifest, self.measurement, root)

    def test_comparison_plot_plan_is_nominal_by_default_and_opt_in(self) -> None:
        entries = [
            ("nominal", "NLO+PS (polarized; full spin)"),
            ("unpolarized_nlo", "NLO+PS (unpolarized beams)"),
            ("polarized_lo", "LO+PS (polarized)"),
            ("no_real_spin_nlo", "NLO+PS (polarized; Born spin only)"),
            ("no_shower_spin_nlo", "NLO+PS (polarized; shower spin off)"),
        ]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            predictions = []
            for family, label in entries:
                relative = Path("postprocess") / f"{family}.yoda"
                output = root / relative
                output.parent.mkdir(parents=True, exist_ok=True)
                output.write_text("non-empty", encoding="utf-8")
                predictions.append(
                    {"family": family, "label": label, "yoda": str(relative)}
                )
            manifest = {
                "configuration": {
                    "comparisons": True,
                    "measurement_signature": campaign.measurement_signature(
                        self.measurement
                    ),
                    "comparison_signature": campaign.comparison_signature(
                        self.measurement
                    ),
                },
                "postprocess": {"predictions": predictions},
                "runtime": {"tools": {"rivet-mkhtml": "/test/rivet-mkhtml"}},
            }
            campaign.atomic_write_json(root / campaign.MANIFEST_NAME, manifest)
            args = campaign.make_parser().parse_args(
                [
                    "plot",
                    "--measurement",
                    MEASUREMENT_ID,
                    "--tag",
                    "comparison-plot-plan",
                    "--dry-run",
                ]
            )
            output = io.StringIO()
            with mock.patch.object(campaign, "campaign_directory", return_value=root):
                with contextlib.redirect_stdout(output):
                    campaign.plot_campaign(args, self.measurement)
            command = output.getvalue()
            family, label = entries[0]
            self.assertIn(f"{family}.yoda:Title={label}", command)
            for family, label in entries[1:]:
                self.assertNotIn(f"{family}.yoda:Title={label}", command)
            self.assertIn("-M /HERMES_2007_I726689/Sigma(LL|UU)_", command)

            comparison_args = campaign.make_parser().parse_args(
                [
                    "plot",
                    "--measurement",
                    MEASUREMENT_ID,
                    "--tag",
                    "comparison-plot-plan",
                    "--dry-run",
                    "--plot-comparisons",
                ]
            )
            comparison_output = io.StringIO()
            with mock.patch.object(campaign, "campaign_directory", return_value=root):
                with contextlib.redirect_stdout(comparison_output):
                    campaign.plot_campaign(comparison_args, self.measurement)
            comparison_command = comparison_output.getvalue()
            for family, label in entries:
                self.assertIn(f"{family}.yoda:Title={label}", comparison_command)
            self.assertIn("-M /HERMES_2007_I726689/SigmaLL_", comparison_command)

    def test_sequential_rivet_plot_indexes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "plots"
            analysis_dir = output / MEASUREMENT_ID
            analysis_dir.mkdir(parents=True)
            scripts = []
            for name in ("A1_Q2GT1", "Aparallel_Q2GT1"):
                script = analysis_dir / f"{name}.py"
                script.write_text("# generated by rivet-mkhtml\n", encoding="utf-8")
                (analysis_dir / f"{name}.png").write_bytes(b"png")
                (analysis_dir / f"{name}.pdf").write_bytes(b"pdf")
                scripts.append(script)
            index = campaign.write_plot_indexes(output, self.measurement, scripts)
            self.assertTrue(index.is_file())
            self.assertIn(f"{MEASUREMENT_ID}/index.html", index.read_text(encoding="utf-8"))
            child = analysis_dir / "index.html"
            child_text = child.read_text(encoding="utf-8")
            self.assertIn("A1_Q2GT1.png", child_text)
            self.assertIn("Aparallel_Q2GT1.pdf", child_text)


if __name__ == "__main__":
    unittest.main()
