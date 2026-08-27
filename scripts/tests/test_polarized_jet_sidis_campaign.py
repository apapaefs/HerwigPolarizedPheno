#!/usr/bin/env python3
"""Focused validation for STAR jets, HERMES SIDIS, and polarized QCD 2->2."""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import io
import json
import math
import re
import sys
import unittest
from pathlib import Path

import numpy as np


DISPOL_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(DISPOL_ROOT / "scripts"))

import polarized_jet_sidis_reference_data as reference  # noqa: E402
import polarized_sidis_postprocess as sidis  # noqa: E402
import run_experimental_campaign as experimental  # noqa: E402
import run_phenomenology_campaign as campaign  # noqa: E402


def plot_metadata(path: Path, object_path: str) -> dict[str, str]:
    metadata: dict[str, str] = {}
    text = path.read_text(encoding="utf-8")
    blocks = re.finditer(
        r"(?ms)^\s*#?\s*BEGIN PLOT\s+(.+?)\s*$\n(.*?)"
        r"^\s*#?\s*END PLOT\s*$",
        text,
    )
    for block in blocks:
        if re.match(block.group(1), object_path) is None:
            continue
        for line in block.group(2).splitlines():
            if "=" in line and not line.lstrip().startswith("#"):
                key, value = line.split("=", 1)
                metadata[key.strip()] = value.strip()
    return metadata


def args(profile: str = "central", **updates: object) -> argparse.Namespace:
    values: dict[str, object] = {
        "profile": profile,
        "jobs": 4,
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
    }
    values.update(updates)
    return argparse.Namespace(**values)


def series(
    values: list[float],
    variances: list[float] | None = None,
) -> experimental.BinSeries:
    return experimental.BinSeries(
        [float(index) for index in range(len(values) + 1)],
        values,
        variances or [0.0 for _ in values],
    )


class ReferenceDataTests(unittest.TestCase):
    def test_hermes_sidis_primary_and_pull_coordinates_are_labeled(self) -> None:
        plot = (
            DISPOL_ROOT / "analyses" / "rivet" / "dis" /
            "HERMES_2019_I1698889.plot"
        )
        expected = {
            "x": "$x$",
            "xz": "Flattened $(x,z)$ bin",
            "xpt": "Flattened $(x,P_{hT})$ bin",
            "xzpt": "Flattened $(x,z,P_{hT})$ bin",
        }
        for suffix, label in expected.items():
            for prefix in ("", "Pull_"):
                stem = f"{prefix}proton_piplus_{suffix}"
                with self.subTest(stem=stem):
                    metadata = plot_metadata(
                        plot, f"/HERMES_2019_I1698889/{stem}"
                    )
                    self.assertEqual(metadata.get("XLabel"), label)
                    self.assertTrue(metadata.get("YLabel"))

    def test_star_primary_and_pull_coordinates_are_labeled(self) -> None:
        plot = (
            DISPOL_ROOT / "analyses" / "rivet" / "pp" /
            "STAR_2022_I1949588.plot"
        )
        expected = {
            "inclusive_ALL": "Parton-jet $p_T$ [GeV]",
            "dijet_A_ALL": "Parton-dijet mass [GeV]",
            "Pull_inclusive": "Parton-jet $p_T$ [GeV]",
            "Pull_dijet_A": "Parton-dijet mass [GeV]",
        }
        for stem, label in expected.items():
            with self.subTest(stem=stem):
                metadata = plot_metadata(
                    plot, f"/STAR_2022_I1949588/{stem}"
                )
                self.assertEqual(metadata.get("XLabel"), label)
                self.assertTrue(metadata.get("YLabel"))

    def test_every_vendored_value_reconstructs_from_pinned_sources(self) -> None:
        for measurement in sorted(reference.MEASUREMENTS):
            generated = reference.normalized_from_raw(measurement)
            vendored = json.loads(
                (
                    DISPOL_ROOT / reference.REFERENCE_PATHS[measurement]
                ).read_text(encoding="utf-8")
            )
            reference.assert_normalized_snapshot_matches(generated, vendored)
            self.assertEqual(reference.validate_vendored(measurement), vendored)

    def test_snapshot_comparison_allows_only_last_bit_float_roundoff(self) -> None:
        reference.assert_normalized_snapshot_matches(
            {"value": [1.0]}, {"value": [1.0 + 8.0e-16]}
        )
        with self.assertRaisesRegex(
            reference.NewReferenceDataError, "float mismatch"
        ):
            reference.assert_normalized_snapshot_matches(
                {"value": [1.0]}, {"value": [1.0 + 1.0e-12]}
            )

    def test_star_complete_table_inventories_and_checksums(self) -> None:
        for measurement, expected_count in (
            ("STAR_2021_I1850855", 21),
            ("STAR_2022_I1949588", 20),
        ):
            manifest = json.loads(
                (
                    DISPOL_ROOT
                    / reference.SOURCE_MANIFEST_PATHS[measurement]
                ).read_text(encoding="utf-8")
            )
            self.assertEqual(manifest["record_version"], 1)
            self.assertEqual(
                [entry["number"] for entry in manifest["tables"]],
                list(range(1, expected_count + 1)),
            )
            record = DISPOL_ROOT / manifest["record_path"]
            self.assertEqual(
                hashlib.sha256(record.read_bytes()).hexdigest(),
                manifest["record_sha256"],
            )
            for entry in manifest["tables"]:
                payload = (DISPOL_ROOT / entry["path"]).read_bytes()
                self.assertEqual(
                    hashlib.sha256(payload).hexdigest(), entry["sha256"]
                )
            if measurement == "STAR_2022_I1949588":
                binning = manifest["analysis_binning"]
                payload = (DISPOL_ROOT / binning["path"]).read_bytes()
                self.assertEqual(
                    hashlib.sha256(payload).hexdigest(), binning["sha256"]
                )
                self.assertEqual(
                    json.loads(payload)["measurement"], measurement
                )

    def test_star_point_counts_bins_covariance_and_globals(self) -> None:
        star200 = reference.validate_vendored("STAR_2021_I1850855")
        self.assertEqual(
            {key: len(value["points"]) for key, value in star200["datasets"].items()},
            {
                "inclusive_forward": 11,
                "inclusive_central": 11,
                "inclusive_combined": 11,
                "dijet_same_sign": 7,
                "dijet_opposite_sign": 7,
            },
        )
        self.assertEqual(
            star200["datasets"]["inclusive_forward"]["bin_edges"],
            [
                6.0,
                7.1,
                8.4,
                9.9,
                11.7,
                13.8,
                16.3,
                19.2,
                22.7,
                26.8,
                31.6,
                37.3,
            ],
        )
        self.assertEqual(
            star200["datasets"]["dijet_same_sign"]["bin_edges"],
            [17.0, 19.0, 23.0, 28.0, 34.0, 41.0, 58.0, 82.0],
        )
        self.assertEqual(star200["primary_covariance"]["dimension"], 36)
        self.assertAlmostEqual(
            star200["global_uncertainties"]["relative_luminosity_absolute"],
            7.0e-4,
        )
        self.assertAlmostEqual(
            star200["global_uncertainties"]["polarization_relative"], 0.061
        )

        star510 = reference.validate_vendored("STAR_2022_I1949588")
        self.assertEqual(
            {key: len(value["points"]) for key, value in star510["datasets"].items()},
            {
                "inclusive": 14,
                "dijet_A": 12,
                "dijet_B": 13,
                "dijet_C": 12,
                "dijet_D": 12,
            },
        )
        self.assertEqual(star510["primary_covariance"]["dimension"], 63)
        self.assertEqual(
            star510["datasets"]["inclusive"]["bin_edges"],
            [
                7.0, 8.2, 9.6, 11.2, 13.1, 15.3, 17.9, 20.9,
                24.5, 28.7, 33.6, 39.3, 46.0, 53.8, 62.8,
            ],
        )
        self.assertEqual(
            star510["datasets"]["dijet_A"]["bin_edges"],
            [
                12.0, 14.0, 17.0, 20.0, 24.0, 29.0, 34.0,
                41.0, 49.0, 59.0, 70.0, 84.0, 101.0,
            ],
        )
        self.assertEqual(
            star510["datasets"]["dijet_D"]["bin_edges"],
            [
                14.0, 17.0, 20.0, 24.0, 29.0, 34.0, 41.0,
                49.0, 59.0, 70.0, 84.0, 101.0, 121.0,
            ],
        )
        first_inclusive = star510["datasets"]["inclusive"]["points"][0]
        self.assertAlmostEqual(
            first_inclusive["coordinate_low"], 0.02719 * 255.0
        )
        self.assertAlmostEqual(
            first_inclusive["coordinate_high"], 0.03391 * 255.0
        )
        self.assertEqual(first_inclusive["low"], 7.0)
        self.assertEqual(first_inclusive["high"], 8.2)
        self.assertEqual(
            len(star510["datasets"]["inclusive"]["display_bin_edges"]),
            15,
        )
        self.assertAlmostEqual(
            star510["global_uncertainties"]["relative_luminosity_absolute"],
            4.7e-4,
        )
        self.assertAlmostEqual(
            star510["global_uncertainties"]["polarization_relative"], 0.064
        )
        reported_systematic = (
            first_inclusive[
                "reported_systematic_including_relative_luminosity"
            ]
        )
        point_to_point_systematic = math.sqrt(
            reported_systematic**2 - (4.7e-4)**2
        )
        self.assertAlmostEqual(reported_systematic, 6.0e-4)
        self.assertAlmostEqual(
            first_inclusive["systematic"], point_to_point_systematic
        )
        self.assertAlmostEqual(
            first_inclusive["systematic_combined"],
            point_to_point_systematic,
        )
        self.assertAlmostEqual(
            first_inclusive["point_to_point_total"],
            math.hypot(first_inclusive["stat"], point_to_point_systematic),
        )
        self.assertAlmostEqual(
            math.hypot(first_inclusive["systematic"], 4.7e-4),
            reported_systematic,
        )

        for snapshot in (star200, star510):
            covariance = np.asarray(
                snapshot["primary_covariance"]["covariance"], dtype=float
            )
            self.assertTrue(np.allclose(covariance, covariance.T, atol=1.0e-14))
            self.assertGreaterEqual(
                float(np.linalg.eigvalsh(covariance)[0]), -1.0e-12
            )
            for index, label in enumerate(
                snapshot["primary_covariance"]["ordering"]
            ):
                dataset, point = label.rsplit(":", 1)
                expected = snapshot["datasets"][dataset]["points"][
                    int(point) - 1
                ]["point_to_point_total"]
                self.assertAlmostEqual(covariance[index, index], expected**2)

    def test_star_510_systematic_only_dijet_block(self) -> None:
        snapshot = reference.validate_vendored("STAR_2022_I1949588")
        manifest = json.loads(
            (
                DISPOL_ROOT
                / reference.SOURCE_MANIFEST_PATHS["STAR_2022_I1949588"]
            ).read_text(encoding="utf-8")
        )
        entry = next(item for item in manifest["tables"] if item["number"] == 12)
        table = json.loads(
            (DISPOL_ROOT / entry["path"]).read_text(encoding="utf-8")
        )
        row = table["values"][0]
        rho = float(row["y"][0]["value"])
        first = snapshot["datasets"]["dijet_A"]["points"][0]
        second = snapshot["datasets"]["dijet_B"]["points"][0]
        expected = rho * first["systematic"] * second["systematic"]
        self.assertAlmostEqual(
            math.hypot(first["systematic"], 4.7e-4),
            first["reported_systematic_including_relative_luminosity"],
        )
        self.assertAlmostEqual(
            math.hypot(second["systematic"], 4.7e-4),
            second["reported_systematic_including_relative_luminosity"],
        )
        ordering = snapshot["primary_covariance"]["ordering"]
        i = ordering.index("dijet_A:1")
        j = ordering.index("dijet_B:1")
        self.assertAlmostEqual(
            snapshot["primary_covariance"]["covariance"][i][j], expected
        )

    def test_hermes_inventory_grids_covariances_and_global_metadata(self) -> None:
        snapshot = reference.validate_vendored("HERMES_2019_I1698889")
        self.assertEqual(
            snapshot["projection_point_counts"],
            {"x": 54, "xz": 126, "xpt": 108, "xzpt": 477, "difference": 36},
        )
        self.assertEqual(len(snapshot["datasets"]), 28)
        self.assertEqual(len(snapshot["diagnostics"]), 6)
        self.assertEqual(
            sum(len(item["points"]) for item in snapshot["diagnostics"].values()),
            180,
        )
        for diagnostic in snapshot["diagnostics"].values():
            self.assertEqual(
                diagnostic["observable"], "A_parallel_cosphi_amplitude"
            )
            self.assertEqual(diagnostic["source_column_label"], "2<cos(phi)>")
            self.assertIn("PUBLISHED_AParallelCosPhi", diagnostic["rivet_path"])
            self.assertIn(
                "aparallel_cosphi_amplitude", diagnostic["points"][0]
            )
            self.assertNotIn("two_cosphi_moment", diagnostic["points"][0])
        self.assertEqual(
            snapshot["binning"]["x"],
            [0.023, 0.04, 0.055, 0.075, 0.1, 0.14, 0.2, 0.3, 0.4, 0.6],
        )
        self.assertEqual(
            snapshot["binning"]["three_dimensional_x"],
            [0.023, 0.04, 0.055, 0.075, 0.1, 0.14, 0.2, 0.3, 0.4, 0.6],
        )
        self.assertEqual(
            snapshot["provenance"]["normalization_overrides"][0][
                "supplemental_header_value"
            ],
            0.45,
        )
        self.assertEqual(
            snapshot["binning"]["azimuthal_x_fine"],
            [0.023, 0.04, 0.055, 0.075, 0.14, 0.6],
        )
        self.assertEqual(
            snapshot["binning"]["azimuthal_z_fine"],
            [0.2, 0.32, 0.44, 0.56, 0.68, 0.8],
        )
        self.assertAlmostEqual(
            snapshot["global_uncertainties"]["proton_polarization_relative"],
            0.066,
        )
        self.assertAlmostEqual(
            snapshot["global_uncertainties"]["deuteron_polarization_relative"],
            0.057,
        )
        missing = []
        for identifier, dataset in snapshot["datasets"].items():
            if dataset["projection"] == "difference":
                self.assertEqual(dataset["observable"], "A1_charge_difference")
                self.assertNotIn("aparallel_statistical_covariance", dataset)
                continue
            covariance = dataset["aparallel_statistical_covariance"]
            if covariance is None:
                missing.append(identifier)
                continue
            size = len(dataset["points"])
            self.assertEqual((len(covariance), len(covariance[0])), (size, size))
            self.assertGreaterEqual(
                float(np.linalg.eigvalsh(np.asarray(covariance))[0]), -1.0e-12
            )
        self.assertEqual(missing, ["deuteron_kminus_xpt"])

    def test_hermes_unpolarized_cosphi_never_overlays_published_spin_fit(
        self,
    ) -> None:
        snapshot = reference.validate_vendored("HERMES_2019_I1698889")
        measurement = campaign.discover_pp_registry()["HERMES_2019_I1698889"]
        identifier, diagnostic = next(iter(snapshot["diagnostics"].items()))
        observable = f"UnpolarizedCosPhi_{identifier}"
        self.assertIsNone(
            campaign._reference_path(measurement, observable, snapshot)
        )
        self.assertIsNone(
            campaign._reference_points(measurement, snapshot, observable)
        )
        self.assertIsNone(
            campaign._pp_reference_overlay_points(
                measurement,
                snapshot,
                Path(diagnostic["rivet_path"]).name,
            )
        )
        plot = (
            DISPOL_ROOT / "analyses/rivet/dis/HERMES_2019_I1698889.plot"
        ).read_text(encoding="utf-8")
        self.assertIn(
            "BEGIN PLOT /HERMES_2019_I1698889/DIAGNOSTICS/"
            "UnpolarizedCosPhi_.*",
            plot,
        )
        self.assertIn(r"YLabel=$2\langle\cos\phi\rangle_{UU}$", plot)
        self.assertIn(
            "BEGIN PLOT /HERMES_2019_I1698889/"
            "PUBLISHED_AParallelCosPhi/.*",
            plot,
        )

    def test_hermes_sparse_3d_rows_map_to_canonical_cells(self) -> None:
        snapshot = reference.validate_vendored("HERMES_2019_I1698889")
        pion_tail = [72, 73, 75, 76, 77, 78, 79, 80]
        expected_tails = {
            "proton_piplus_xzpt": pion_tail,
            "proton_piminus_xzpt": pion_tail,
            "deuteron_piplus_xzpt": pion_tail,
            "deuteron_piminus_xzpt": pion_tail,
            "deuteron_kplus_xzpt": [72, 75, 76, 77, 78, 79, 80],
            "deuteron_kminus_xzpt": [72, 75, 76, 77, 78, 79],
        }
        for identifier, expected_tail in expected_tails.items():
            points = snapshot["datasets"][identifier]["points"]
            self.assertEqual(
                [point["published_row"] for point in points],
                list(range(1, len(points) + 1)),
            )
            self.assertEqual(
                [point["flat_bin"] for point in points[:72]],
                list(range(72)),
            )
            self.assertEqual(
                [point["flat_bin"] for point in points[72:]],
                expected_tail,
            )
            self.assertEqual(
                len({point["flat_bin"] for point in points}),
                len(points),
            )
            for point in points:
                self.assertEqual(
                    point["flat_bin"],
                    reference._hermes_xzpt_flat_bin(
                        point["x"], point["z"], point["pt"]
                    ),
                )

    def test_hermes_archive_checksum_inventory_and_html_rejection(self) -> None:
        archive = DISPOL_ROOT / reference.HERMES_ARCHIVE
        self.assertEqual(
            hashlib.sha256(archive.read_bytes()).hexdigest(),
            reference.HERMES_ARCHIVE_SHA256,
        )
        members, _ = reference._validated_hermes_members(archive.read_bytes())
        self.assertEqual(len(members), 57)
        with self.assertRaisesRegex(
            reference.NewReferenceDataError, "not a ZIP archive"
        ):
            reference._zip_members(b"<html>Cloudflare challenge</html>")


class CampaignMatrixAndCardTests(unittest.TestCase):
    def test_postprocessing_refuses_stale_measurement_signatures(self) -> None:
        measurement = campaign.discover_pp_registry()["HERMES_2019_I1698889"]
        current = {
            "configuration": {
                "measurement_signature": campaign._signature(measurement)
            }
        }
        campaign._assert_manifest_signature_current(current, measurement)
        stale = {
            "configuration": {"measurement_signature": "stale-analysis"}
        }
        with self.assertRaisesRegex(
            campaign.CampaignError, "intentionally refused"
        ):
            campaign._assert_manifest_signature_current(stale, measurement)

    def test_new_job_matrices_and_family_selection(self) -> None:
        registry = campaign.discover_pp_registry()
        for measurement_id in (
            "STAR_2021_I1850855",
            "STAR_2022_I1949588",
        ):
            measurement = registry[measurement_id]
            central = campaign.build_job_matrix(
                measurement, campaign._resolved_options(args(), measurement)
            )
            paper = campaign.build_job_matrix(
                measurement,
                campaign._resolved_options(args("paper"), measurement),
            )
            all_families = campaign.build_job_matrix(
                measurement,
                campaign._resolved_options(
                    args("paper", families="all"), measurement
                ),
            )
            self.assertEqual(len(central), 4)
            self.assertEqual(len(paper), 812)
            self.assertEqual(len(all_families), 821)
            self.assertEqual({job["family"] for job in paper}, {"nominal"})
            self.assertEqual(
                {job["analysis_instance"] for job in central},
                {f"{measurement_id}:LEVEL=HARDPARTON"},
            )
            self.assertEqual(
                {job["family"] for job in all_families},
                {
                    "nominal",
                    "shower_spin_off",
                    "hadron_mpi_on",
                    "unpolarized_closure",
                },
            )
            hadron = next(
                job
                for job in all_families
                if job["family"] == "hadron_mpi_on"
            )
            self.assertEqual(
                hadron["analysis_instance"],
                f"{measurement_id}:LEVEL=HADRON",
            )

        hermes = registry["HERMES_2019_I1698889"]
        central = campaign.build_job_matrix(
            hermes, campaign._resolved_options(args(), hermes)
        )
        paper = campaign.build_job_matrix(
            hermes, campaign._resolved_options(args("paper"), hermes)
        )
        self.assertEqual(len(central), 16)
        self.assertEqual(len(paper), 3248)
        self.assertEqual({job["target_component"] for job in central}, {"P", "N"})
        self.assertEqual(
            {job["contribution"] for job in central}, {"POSNLO", "NEGNLO"}
        )
        for jobs in (central, paper):
            self.assertEqual(len({job["id"] for job in jobs}), len(jobs))
            self.assertEqual(len({job["seed"] for job in jobs}), len(jobs))

    def test_scale_mapping_and_star_card_invariants(self) -> None:
        measurement = campaign.discover_pp_registry()["STAR_2021_I1850855"]
        options = campaign._resolved_options(
            args(
                "paper",
                polarized_pdf_members="7",
                unpolarized_pdf_members="9",
                scales="2",
            ),
            measurement,
        )
        jobs = campaign.build_job_matrix(measurement, options)
        scale_job = next(job for job in jobs if job["scale"] == 2.0)
        text = campaign._card_text(measurement, scale_job)
        self.assertIn("MEQCD2to2:ScalePreFactor 4", text)
        self.assertIn("STAR200JetPolarizedPDF:Member 0", text)
        replica_job = next(
            job for job in jobs if job["polarized_pdf_member"] == 7
        )
        replica_text = campaign._card_text(measurement, replica_job)
        self.assertIn("STAR200JetPolarizedPDF:Member 7", replica_text)
        self.assertIn("HardLOPDF:Member 0", replica_text)
        common = (
            DISPOL_ROOT
            / measurement["cards"]["directory"]
            / measurement["cards"]["common"]
        ).read_text(encoding="utf-8")
        self.assertIn("NNPDF40_nlo_pch_as_01180", common)
        self.assertIn("NNPDFpol20_nlo_as_01180", common)
        self.assertIn("SpinCorrelations Yes", common)
        self.assertIn("Interactions QCD", common)
        self.assertIn("MPIHandler NULL", common)
        base = (
            DISPOL_ROOT
            / measurement["cards"]["directory"]
            / "STAR_2021_I1850855_jets_PP-LO.in"
        ).read_text(encoding="utf-8")
        self.assertIn("/Herwig/MatrixElements/MEQCD2to2", base)
        self.assertNotIn("MEQCD2to2Fast", base)

    def test_star_shower_spin_off_family_changes_only_the_shower_control(
        self,
    ) -> None:
        registry = campaign.discover_pp_registry()
        for measurement_id in (
            "STAR_2021_I1850855",
            "STAR_2022_I1949588",
        ):
            measurement = registry[measurement_id]
            options = campaign._resolved_options(
                args(families="nominal,shower_spin_off"), measurement
            )
            jobs = campaign.build_job_matrix(measurement, options)
            self.assertEqual(len(jobs), 8)
            self.assertEqual(
                {job["family"] for job in jobs},
                {"nominal", "shower_spin_off"},
            )
            control_only = campaign.build_job_matrix(
                measurement,
                campaign._resolved_options(
                    args(families="shower_spin_off"), measurement
                ),
            )
            self.assertEqual(len(control_only), 4)
            self.assertEqual(
                {job["family"] for job in control_only},
                {"shower_spin_off"},
            )
            self.assertEqual(
                {
                    job["helicity"]
                    for job in jobs
                    if job["family"] == "shower_spin_off"
                },
                {"PP", "PM", "MP", "MM"},
            )
            self.assertTrue(
                all(
                    job["shower_spin_correlations"] == "off"
                    for job in jobs
                    if job["family"] == "shower_spin_off"
                )
            )
            spin_off = next(
                job
                for job in jobs
                if job["family"] == "shower_spin_off"
                and job["helicity"] == "PP"
            )
            nominal = next(
                job
                for job in jobs
                if job["family"] == "nominal"
                and job["helicity"] == "PP"
            )
            card = campaign._card_text(measurement, spin_off)
            self.assertIn(
                "set /Herwig/Shower/ShowerHandler:SpinCorrelations No",
                card,
            )
            self.assertIn("FirstLongitudinalPolarization", card)
            self.assertIn("SecondLongitudinalPolarization", card)
            self.assertIn("/Herwig/MatrixElements/MEQCD2to2", card)
            normalized_control = card.replace(
                "set /Herwig/Shower/ShowerHandler:SpinCorrelations No\n",
                "",
            ).replace(spin_off["stem"], nominal["stem"])
            self.assertEqual(
                normalized_control, campaign._card_text(measurement, nominal)
            )

            family = measurement["families"]["shower_spin_off"]
            argument = campaign._prediction_plot_argument(
                Path("prediction-shower_spin_off.yoda"),
                family["label"],
                family,
            )
            self.assertIn(":LineColor=#0077BB", argument)
            self.assertIn("polarized hard; shower spin off", argument)

    def test_star_generator_cut_is_immutable_and_star_only(self) -> None:
        registry = campaign.discover_pp_registry()
        measurement = registry["STAR_2022_I1949588"]
        nominal = campaign._resolved_options(args(), measurement)
        self.assertEqual(nominal["jet_kt_min_gev"], 4.0)
        nominal_jobs = campaign.build_job_matrix(measurement, nominal)
        self.assertTrue(
            all("-kt4gev-" in job["id"] for job in nominal_jobs)
        )
        self.assertTrue(
            all(job["jet_kt_min_gev"] == 4.0 for job in nominal_jobs)
        )
        self.assertIn(
            "set /Herwig/Cuts/JetKtCut:MinKT 4.0*GeV",
            campaign._card_text(measurement, nominal_jobs[0]),
        )

        scan = campaign._resolved_options(
            args(jet_kt_min_gev=3.0), measurement
        )
        scan_jobs = campaign.build_job_matrix(measurement, scan)
        self.assertTrue(all("-kt3gev-" in job["id"] for job in scan_jobs))
        self.assertIn(
            "set /Herwig/Cuts/JetKtCut:MinKT 3.0*GeV",
            campaign._card_text(measurement, scan_jobs[0]),
        )

        with self.assertRaisesRegex(campaign.CampaignError, "3, 4, and 5"):
            campaign._resolved_options(
                args(jet_kt_min_gev=2.0), measurement
            )
        with self.assertRaisesRegex(campaign.CampaignError, "only valid"):
            campaign._resolved_options(
                args(jet_kt_min_gev=4.0),
                registry["HERMES_2019_I1698889"],
            )

        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(
                campaign.main(
                    [
                        "prepare",
                        "--measurement",
                        "COMPASS_2010_I843494",
                        "--tag",
                        "unit-invalid-cut",
                        "--jet-kt-min-gev",
                        "4",
                        "--dry-run",
                    ]
                ),
                2,
            )

    def test_sidis_card_invariants_and_scale_mapping(self) -> None:
        measurement = campaign.discover_pp_registry()["HERMES_2019_I1698889"]
        options = campaign._resolved_options(
            args("paper", scales="0.5"), measurement
        )
        jobs = campaign.build_job_matrix(measurement, options)
        scale_job = next(job for job in jobs if job["scale"] == 0.5)
        text = campaign._card_text(measurement, scale_job)
        for scale_object in measurement["cards"]["scale_objects"]:
            self.assertIn(f"{scale_object}:ScaleFactor 0.5", text)
        common = (
            DISPOL_ROOT
            / measurement["cards"]["directory"]
            / measurement["cards"]["common"]
        ).read_text(encoding="utf-8")
        for token in (
            "BeamEMaxA 27.6*GeV",
            "GammaZ Gamma",
            "UsePOWHEGRealSpinVertex Yes",
            "UseNativeDISWindowGeneration Yes",
            "SpinCorrelations Yes",
            "Interactions QCD",
            "MaxQ2 1.0e10*GeV2",
        ):
            self.assertIn(token, common)
        neutron = next(
            job
            for job in jobs
            if job["target_component"] == "N"
            and job["contribution"] == "NEGNLO"
        )
        neutron_text = campaign._card_text(measurement, neutron)
        self.assertIn("TargetParticle /Herwig/Particles/n0", neutron_text)
        self.assertIn("Contribution NegativeNLO", neutron_text)

    def test_parser_defaults_and_comparison_alias_contract(self) -> None:
        parsed = campaign.make_parser().parse_args(
            [
                "full",
                "--measurement",
                "STAR_2021_I1850855",
                "--tag",
                "unit",
                "--profile",
                "paper",
            ]
        )
        self.assertIsNone(parsed.families)
        self.assertFalse(parsed.comparisons)
        self.assertIsNone(parsed.nominal_prediction)
        explicit = campaign.make_parser().parse_args(
            [
                "full",
                "--measurement",
                "STAR_2021_I1850855",
                "--tag",
                "unit",
                "--families",
                "nominal,hadron_mpi_on",
                "--include-diagnostics",
            ]
        )
        self.assertEqual(explicit.families, "nominal,hadron_mpi_on")
        self.assertTrue(explicit.include_diagnostics)
        reuse = campaign.make_parser().parse_args(
            [
                "full",
                "--measurement",
                "STAR_2022_I1949588",
                "--tag",
                "reuse",
                "--families",
                "shower_spin_off",
                "--seed-base",
                "2949588",
                "--plot-comparisons",
                "--nominal-prediction",
                "/tmp/nominal/prediction.yoda",
            ]
        )
        self.assertEqual(
            reuse.nominal_prediction, Path("/tmp/nominal/prediction.yoda")
        )
        measurement = campaign.discover_pp_registry()["STAR_2022_I1949588"]
        with self.assertRaisesRegex(
            campaign.CampaignError, "comparison-family-only"
        ):
            campaign._resolved_options(
                args(
                    families="nominal,shower_spin_off",
                    nominal_prediction=Path("/tmp/nominal/prediction.yoda"),
                    plot_comparisons=True,
                ),
                measurement,
            )
        with self.assertRaisesRegex(
            campaign.CampaignError, "requires --plot-comparisons"
        ):
            campaign._resolved_options(
                args(
                    families="shower_spin_off",
                    nominal_prediction=Path("/tmp/nominal/prediction.yoda"),
                    seed_base=2949588,
                    plot_comparisons=False,
                ),
                measurement,
            )
        with self.assertRaisesRegex(
            campaign.CampaignError, "explicit --seed-base"
        ):
            campaign._resolved_options(
                args(
                    families="shower_spin_off",
                    nominal_prediction=Path("/tmp/nominal/prediction.yoda"),
                    plot_comparisons=True,
                ),
                measurement,
            )

    def test_star_operational_bins_are_preserved_before_display_mapping(
        self,
    ) -> None:
        snapshot = reference.validate_vendored("STAR_2022_I1949588")
        dataset = snapshot["datasets"]["inclusive"]
        predictions = {
            "inclusive": {
                "edges": list(dataset["bin_edges"]),
                "values": [0.0] * len(dataset["points"]),
                "errors": [0.0] * len(dataset["points"]),
            }
        }
        campaign._apply_star_display_binning(predictions, snapshot)
        self.assertEqual(
            predictions["inclusive"]["analysis_edges"],
            dataset["bin_edges"],
        )
        self.assertEqual(
            predictions["inclusive"]["edges"],
            dataset["display_bin_edges"],
        )

        wrong = {
            "inclusive": {
                "edges": [0.0] + list(dataset["bin_edges"])[1:],
                "values": [0.0] * len(dataset["points"]),
                "errors": [0.0] * len(dataset["points"]),
            }
        }
        with self.assertRaisesRegex(
            campaign.CampaignError, "analysis bins disagree"
        ):
            campaign._apply_star_display_binning(wrong, snapshot)


class StarComparisonPolicyTests(unittest.TestCase):
    def setUp(self) -> None:
        self.measurement = campaign.discover_pp_registry()["STAR_2022_I1949588"]
        self.snapshot = reference.validate_vendored("STAR_2022_I1949588")

    def test_inclusive_mask_starts_at_analysis_bin_five(self) -> None:
        size = len(self.snapshot["datasets"]["inclusive"]["points"])
        mask = campaign._comparison_mask(
            self.measurement, self.snapshot, "inclusive", size
        )
        self.assertEqual(mask[:6], [False, False, False, False, True, True])
        self.assertEqual(sum(mask), 10)
        self.assertEqual(
            self.measurement["spin_density_policy"],
            {
                "name": "radial_bloch_ball_projection",
                "status": "interim_star_production",
                "herwigpol_source_commit":
                "ad0c5486ad137e0e34c6274959c1fbf8323c5937",
                "strict_negative_isr_guard": True,
                "event_or_branching_veto": False,
            },
        )

        invalid = json.loads(json.dumps(self.measurement))
        invalid["comparison_policy"]["primary_bin_masks"]["inclusive"][
            "minimum_analysis_low_gev"
        ] = 13.2
        with self.assertRaisesRegex(campaign.CampaignError, "not the configured"):
            campaign._comparison_mask(
                invalid, self.snapshot, "inclusive", size
            )

    def test_primary_pulls_and_overlays_exclude_low_bins(self) -> None:
        dataset = self.snapshot["datasets"]["inclusive"]
        size = len(dataset["points"])
        prediction = {
            "values": [float(index) for index in range(size)],
            "errors": [0.0] * size,
        }
        mask = campaign._comparison_mask(
            self.measurement, self.snapshot, "inclusive", size
        )
        pulls, _chi2, count = campaign._pulls(
            prediction, dataset["points"], mask
        )
        self.assertEqual(pulls[:4], [None, None, None, None])
        self.assertEqual(count, 10)

        points = campaign._pp_reference_overlay_points(
            self.measurement,
            self.snapshot,
            Path(dataset["rivet_path"]).name,
        )
        self.assertIsNotNone(points)
        self.assertEqual(len(points or []), 10)

        summary = {
            "uncertainties": {
                "inclusive": {
                    "monte_carlo": [1.0] * size,
                    "pdf_68": [2.0] * size,
                }
            }
        }
        bands = campaign._pp_theory_uncertainty_bands(
            self.measurement,
            self.snapshot,
            summary,
            Path(dataset["rivet_path"]).name,
        )
        self.assertEqual(bands["monte_carlo"][:4], [None] * 4)
        self.assertEqual(bands["pdf_68"][4:], [2.0] * 10)

    def test_correlated_fit_selects_exactly_59_points(self) -> None:
        prediction = {
            observable: {
                "values": [0.0] * len(dataset["points"]),
                "errors": [0.0] * len(dataset["points"]),
            }
            for observable, dataset in self.snapshot["datasets"].items()
            if not dataset.get("alternate_projection")
        }
        result = campaign._star_correlated_goodness_of_fit(
            prediction, self.snapshot, self.measurement
        )
        self.assertEqual(result["points"], 59)
        self.assertNotIn("inclusive:1", result["ordering"])
        self.assertNotIn("inclusive:4", result["ordering"])
        self.assertIn("inclusive:5", result["ordering"])

    def test_smoke_fit_reports_partial_coverage_without_weakening_production(
        self,
    ) -> None:
        size = len(self.snapshot["datasets"]["inclusive"]["points"])
        values: list[float | None] = [None] * size
        errors: list[float | None] = [None] * size
        values[4] = 0.0
        errors[4] = 1.0
        partial = {"inclusive": {"values": values, "errors": errors}}
        with self.assertRaisesRegex(campaign.CampaignError, "expected 59"):
            campaign._star_correlated_goodness_of_fit(
                partial, self.snapshot, self.measurement
            )
        result = campaign._star_correlated_goodness_of_fit(
            partial,
            self.snapshot,
            self.measurement,
            require_complete=False,
        )
        self.assertEqual(result["points"], 1)
        self.assertEqual(result["expected_points"], 59)
        self.assertEqual(result["coverage_status"], "partial smoke sample")
        self.assertIn("not evaluated", result["status"])
        self.assertNotIn("chi2_profiled", result)


class ObservableArithmeticTests(unittest.TestCase):
    def test_posneg_charge_difference_and_covariance(self) -> None:
        positive = series([10.0], [4.0])
        negative = series([3.0], [9.0])
        covariance = series([0.0], [2.0])
        difference = sidis.difference_series(
            positive, negative, covariance
        )
        self.assertEqual(difference.values, [7.0])
        self.assertEqual(difference.variances, [9.0])
        combined = sidis.add_nlo_components(
            {
                "POSNLO": series([5.0], [1.0]),
                "NEGNLO": series([-2.0], [0.25]),
            }
        )
        self.assertEqual(combined.values, [3.0])
        self.assertEqual(combined.variances, [1.25])

    def test_proton_deuteron_helicity_arithmetic_and_masking(self) -> None:
        proton_samples = {
            "P:PP": series([12.0], [0.4]),
            "P:PM": series([8.0], [0.3]),
            "P:MP": series([8.0], [0.2]),
            "P:MM": series([12.0], [0.5]),
        }
        proton = sidis.asymmetry(proton_samples, "proton")
        self.assertAlmostEqual(proton["values"][0], 0.2)
        deuteron_samples = {
            **proton_samples,
            "N:PP": series([8.0], [0.2]),
            "N:PM": series([12.0], [0.3]),
            "N:MP": series([12.0], [0.4]),
            "N:MM": series([8.0], [0.5]),
        }
        deuteron = sidis.asymmetry(deuteron_samples, "deuteron")
        self.assertAlmostEqual(deuteron["values"][0], 0.0)
        empty = {
            label: series([0.0], [0.0])
            for label in ("P:PP", "P:PM", "P:MP", "P:MM")
        }
        masked = sidis.asymmetry(empty, "proton")
        self.assertIsNone(masked["values"][0])
        self.assertIsNone(masked["errors"][0])

    def test_published_deuteron_correction_is_not_reapplied_to_theory(self) -> None:
        denominator, longitudinal = sidis.target_coefficients("deuteron")
        self.assertTrue(all(value == 0.125 for value in denominator.values()))
        self.assertAlmostEqual(longitudinal["P:PP"], 0.125)
        self.assertAlmostEqual(longitudinal["N:PM"], -0.125)

    def test_star_all_and_closure_signs(self) -> None:
        samples = {
            "PP": series([12.0], [0.4]),
            "PM": series([8.0], [0.3]),
            "MP": series([8.0], [0.2]),
            "MM": series([12.0], [0.5]),
        }
        prediction = campaign._star_jet_prediction({"inclusive": samples})
        self.assertAlmostEqual(prediction["inclusive"]["values"][0], 0.2)
        self.assertAlmostEqual(
            prediction["SigmaUU_inclusive"]["values"][0], 10.0
        )
        self.assertAlmostEqual(
            prediction["SigmaUU_inclusive"]["errors"][0],
            math.sqrt((0.4 + 0.3 + 0.2 + 0.5) / 16.0),
        )
        self.assertAlmostEqual(
            prediction["SingleSpinA_inclusive"]["values"][0], 0.0
        )
        self.assertAlmostEqual(
            prediction["SingleSpinB_inclusive"]["values"][0], 0.0
        )
        self.assertAlmostEqual(
            prediction["Parity_PP_MM_inclusive"]["values"][0], 0.0
        )

    def test_star_global_nuisances_are_profiled_once(self) -> None:
        result = campaign._star_correlated_goodness_of_fit(
            {"inclusive": {"values": [1.0], "errors": [0.0]}},
            {
                "datasets": {
                    "inclusive": {"points": [{"value": 0.0}]}
                },
                "primary_covariance": {
                    "ordering": ["inclusive:1"],
                    "covariance": [[1.0]],
                },
                "global_uncertainties": {
                    "relative_luminosity_absolute": 1.0,
                    "polarization_relative": 0.0,
                },
            },
        )
        self.assertAlmostEqual(
            result["chi2_correlated_without_global_nuisances"], 1.0
        )
        self.assertAlmostEqual(result["nuisance_pulls"]["relative_luminosity"], 0.5)
        self.assertAlmostEqual(result["nuisance_pulls"]["polarization"], 0.0)
        self.assertAlmostEqual(result["chi2_profiled"], 0.5)
        self.assertEqual(
            result["covariance"],
            "published point-to-point covariance plus diagonal Monte Carlo "
            "statistical variance",
        )

    def test_replica_band_uses_sample_variance_and_independent_quadrature(
        self,
    ) -> None:
        central_key = ("nominal", "jets", 0, 0, 1.0, "off")
        predictions = {
            central_key: {
                "inclusive": {
                    "edges": [0.0, 1.0],
                    "values": [10.0],
                    "errors": [0.5],
                }
            },
            ("nominal", "jets", 1, 0, 1.0, "off"): {
                "inclusive": {
                    "edges": [0.0, 1.0],
                    "values": [9.0],
                    "errors": [0.5],
                }
            },
            ("nominal", "jets", 2, 0, 1.0, "off"): {
                "inclusive": {
                    "edges": [0.0, 1.0],
                    "values": [13.0],
                    "errors": [0.5],
                }
            },
        }
        sigma = campaign._replica_sigma(
            predictions,
            [
                ("nominal", "jets", 1, 0, 1.0, "off"),
                ("nominal", "jets", 2, 0, 1.0, "off"),
            ],
            "inclusive",
            [10.0],
        )
        self.assertAlmostEqual(sigma[0], math.sqrt(8.0))


class DensityMatrixAndSourceContractTests(unittest.TestCase):
    @staticmethod
    def contract(
        amplitudes: np.ndarray, rho_first: np.ndarray, rho_second: np.ndarray
    ) -> float:
        value = 4.0 * np.einsum(
            "ai,bj,ijo,abo->",
            rho_first,
            rho_second,
            amplitudes,
            amplitudes.conjugate(),
        )
        return float(np.real(value))

    @staticmethod
    def interference(
        first: np.ndarray,
        second: np.ndarray,
        rho_first: np.ndarray,
        rho_second: np.ndarray,
    ) -> float:
        value = 4.0 * np.einsum(
            "ai,bj,ijo,abo->",
            rho_first,
            rho_second,
            first,
            second.conjugate(),
        )
        return float(np.real(value))

    @staticmethod
    def physical_spin_one(rho: np.ndarray) -> np.ndarray:
        physical = np.zeros((3, 3), dtype=complex)
        block = np.ix_([0, 2], [0, 2])
        trace = float(np.real(np.trace(rho[block])))
        if trace > 0.0:
            physical[block] = rho[block] / trace
        return physical

    def test_generic_spin_one_identity_projects_to_two_transverse_states(
        self,
    ) -> None:
        generic = np.eye(3, dtype=complex) / 3.0
        physical = self.physical_spin_one(generic)
        self.assertTrue(
            np.allclose(
                physical,
                np.diag([0.5, 0.0, 0.5]).astype(complex),
            )
        )
        polarized = np.asarray(
            [
                [0.7, 0.0, 0.15j],
                [0.0, 0.0, 0.0],
                [-0.15j, 0.0, 0.3],
            ],
            dtype=complex,
        )
        self.assertTrue(
            np.allclose(self.physical_spin_one(polarized), polarized)
        )

    def test_general_density_contraction_unpolarized_limit_and_positivity(
        self,
    ) -> None:
        random = np.random.default_rng(1850855)
        for _ in range(50):
            amplitudes = (
                random.normal(size=(2, 2, 5))
                + 1j * random.normal(size=(2, 2, 5))
            )
            identity = np.eye(2, dtype=complex) / 2.0
            unpolarized = self.contract(amplitudes, identity, identity)
            self.assertAlmostEqual(
                unpolarized, float(np.sum(np.abs(amplitudes) ** 2)), places=11
            )
            pure_values = []
            for first in range(2):
                for second in range(2):
                    rho_first = np.zeros((2, 2), dtype=complex)
                    rho_second = np.zeros((2, 2), dtype=complex)
                    rho_first[first, first] = 1.0
                    rho_second[second, second] = 1.0
                    pure_values.append(
                        self.contract(amplitudes, rho_first, rho_second)
                    )
            self.assertAlmostEqual(
                sum(pure_values) / 4.0, unpolarized, places=11
            )
            vectors = []
            for _beam in range(2):
                vector = (
                    random.normal(size=2) + 1j * random.normal(size=2)
                )
                vector /= np.linalg.norm(vector)
                vectors.append(vector)
            rho_first = np.outer(vectors[0], vectors[0].conjugate())
            rho_second = np.outer(vectors[1], vectors[1].conjugate())
            self.assertTrue(np.allclose(rho_first, rho_first.conjugate().T))
            self.assertTrue(np.allclose(rho_second, rho_second.conjugate().T))
            self.assertGreaterEqual(
                self.contract(amplitudes, rho_first, rho_second), -1.0e-11
            )

    def test_beam_exchange_and_normalized_flow_probabilities(self) -> None:
        random = np.random.default_rng(1949588)
        amplitude = (
            random.normal(size=(2, 2, 4))
            + 1j * random.normal(size=(2, 2, 4))
        )
        symmetric = amplitude + np.swapaxes(amplitude, 0, 1)
        rho_a = np.asarray([[0.7, 0.2j], [-0.2j, 0.3]], dtype=complex)
        rho_b = np.asarray([[0.4, 0.1], [0.1, 0.6]], dtype=complex)
        self.assertAlmostEqual(
            self.contract(symmetric, rho_a, rho_b),
            self.contract(symmetric, rho_b, rho_a),
            places=11,
        )
        flow_weights = np.asarray(
            [
                self.contract(symmetric, rho_a, rho_b),
                self.contract(0.7 * symmetric, rho_a, rho_b),
                self.contract(0.2j * symmetric, rho_a, rho_b),
            ]
        )
        probabilities = flow_weights / np.sum(flow_weights)
        self.assertTrue(np.all(probabilities >= 0.0))
        self.assertAlmostEqual(float(np.sum(probabilities)), 1.0)

    def test_every_qcd_channel_color_contraction_has_legacy_unpolarized_limit(
        self,
    ) -> None:
        """Check every production color topology with an independent sum.

        The random tensors stand for the fixed-point complex helicity
        amplitudes already evaluated by Herwig's QCD vertices.  The reference
        expressions below are the pre-change raw helicity sums and do not call
        the production contraction helper.
        """

        random = np.random.default_rng(20260722)
        identity = np.eye(2, dtype=complex) / 2.0
        pure = []
        for first in range(2):
            for second in range(2):
                rho_first = np.zeros((2, 2), dtype=complex)
                rho_second = np.zeros((2, 2), dtype=complex)
                rho_first[first, first] = 1.0
                rho_second[second, second] = 1.0
                pure.append((first, second, rho_first, rho_second))

        def tensor() -> np.ndarray:
            return (
                random.normal(size=(2, 2, 4))
                + 1j * random.normal(size=(2, 2, 4))
            )

        channel_topologies = {
            "gg2qqbar": "two_flow",
            "qqbar2gg": "two_flow",
            "qg2qg": "two_flow",
            "qbarg2qbarg": "two_flow",
            "gg2gg": "three_flow",
            "qq2qq": "two_diagram",
            "qbarqbar2qbarqbar": "two_diagram",
            "qqbar2qqbar": "two_diagram",
        }
        c1 = 4.0 * (9.0**2 / 8.0 - 3.0 * 9.0 / 8.0 + 1.0 - 0.75 / 9.0)
        c2 = 4.0 * (-0.25 * 9.0 + 1.0 - 0.75 / 9.0)
        for channel, topology in channel_topologies.items():
            count = 3 if topology == "three_flow" else 2
            amplitudes = [tensor() for _ in range(count)]

            def contracted(rho_first: np.ndarray, rho_second: np.ndarray) -> float:
                norms = [
                    self.contract(item, rho_first, rho_second)
                    for item in amplitudes
                ]
                if topology == "two_flow":
                    return (
                        norms[0]
                        + norms[1]
                        - 0.25
                        * self.interference(
                            amplitudes[0],
                            amplitudes[1],
                            rho_first,
                            rho_second,
                        )
                    )
                if topology == "two_diagram":
                    return (
                        norms[0]
                        + norms[1]
                        + 2.0
                        / 3.0
                        * self.interference(
                            amplitudes[0],
                            amplitudes[1],
                            rho_first,
                            rho_second,
                        )
                    )
                return c1 * sum(norms) + 2.0 * c2 * sum(
                    self.interference(
                        amplitudes[first],
                        amplitudes[second],
                        rho_first,
                        rho_second,
                    )
                    for first, second in ((0, 1), (0, 2), (1, 2))
                )

            legacy_norms = [
                float(np.sum(np.abs(item) ** 2)) for item in amplitudes
            ]
            if topology == "two_flow":
                legacy = (
                    legacy_norms[0]
                    + legacy_norms[1]
                    - 0.25
                    * float(
                        np.real(
                            np.sum(amplitudes[0] * amplitudes[1].conjugate())
                        )
                    )
                )
            elif topology == "two_diagram":
                legacy = (
                    legacy_norms[0]
                    + legacy_norms[1]
                    + 2.0
                    / 3.0
                    * float(
                        np.real(
                            np.sum(amplitudes[0] * amplitudes[1].conjugate())
                        )
                    )
                )
            else:
                legacy = c1 * sum(legacy_norms) + 2.0 * c2 * sum(
                    float(
                        np.real(
                            np.sum(
                                amplitudes[first]
                                * amplitudes[second].conjugate()
                            )
                        )
                    )
                    for first, second in ((0, 1), (0, 2), (1, 2))
                )
            self.assertAlmostEqual(
                contracted(identity, identity),
                legacy,
                places=10,
                msg=channel,
            )
            self.assertAlmostEqual(
                sum(
                    contracted(rho_first, rho_second)
                    for _, _, rho_first, rho_second in pure
                )
                / 4.0,
                legacy,
                places=10,
                msg=channel,
            )

class RivetSourceContractTests(unittest.TestCase):
    def test_every_fixed_target_analysis_uses_a_prompt_scattered_lepton(self) -> None:
        contracts = {
            "COMPASS_2010_I843494.cc": (-13, "PromptMuons"),
            "COMPASS_2016_I1357198.cc": (-13, "PromptMuons"),
            "COMPASS_2017_I1501480.cc": (-13, "PromptMuons"),
            "HERMES_2007_I726689.cc": (-11, "PromptPositrons"),
            "HERMES_2007_I726689_LEGACY.cc": (-11, "PromptPositrons"),
            "HERMES_2019_I1698889.cc": (-11, "PromptPositrons"),
        }
        for filename, (pid, projection) in contracts.items():
            source = (
                DISPOL_ROOT / "analyses" / "rivet" / "dis" / filename
            ).read_text(encoding="utf-8")
            self.assertIn("PromptFinalState", source, filename)
            self.assertIn(
                f'PromptFinalState(Cuts::pid == {pid}), "{projection}"',
                source,
                filename,
            )
            self.assertIn(
                f'apply<PromptFinalState>(event, "{projection}")',
                source,
                filename,
            )

    def test_star_algorithms_provenance_bins_and_topologies(self) -> None:
        helper = (
            DISPOL_ROOT / "analyses/rivet/pp/STARPolarizedJets.hh"
        ).read_text(encoding="utf-8")
        star200 = (
            DISPOL_ROOT / "analyses/rivet/pp/STAR_2021_I1850855.cc"
        ).read_text(encoding="utf-8")
        star510 = (
            DISPOL_ROOT / "analyses/rivet/pp/STAR_2022_I1949588.cc"
        ).read_text(encoding="utf-8")
        self.assertIn("fastjet::antikt_algorithm", helper)
        self.assertIn("signal_process_vertex", helper)
        self.assertIn("hasRemnantAncestor", helper)
        self.assertIn("No all-event fallback", helper)
        self.assertIn("Only outgoing hard-process lines", helper)
        self.assertNotIn("hardLineAncestors", helper)
        self.assertNotIn(
            "terminalPartons(incoming", helper
        )
        self.assertIn("hardPartonJets(event, 0.6", star200)
        self.assertIn("first.pT()/GeV <= 8.0", star200)
        self.assertIn("second.pT()/GeV <= 6.0", star200)
        self.assertIn("std::abs(first.eta()) >= 0.8", star200)
        self.assertIn("etaProduct > 0.0", star200)
        self.assertIn("hardPartonJets(event, 0.5", star510)
        for exact_bins in (
            "{7.0,8.2,9.6,11.2,13.1,15.3,17.9,20.9,24.5,28.7,",
            "{12.0,14.0,17.0,20.0,24.0,29.0,34.0,41.0,49.0,",
            "{14.0,17.0,20.0,24.0,29.0,34.0,41.0,49.0,59.0,",
        ):
            self.assertIn(exact_bins, star510)
        self.assertIn("first.pT()/GeV <= 7.0", star510)
        self.assertIn("second.pT()/GeV <= 5.0", star510)
        self.assertIn("std::abs(first.eta()-second.eta()) >= 1.6", star510)
        for name in ("_dijetA", "_dijetB", "_dijetC", "_dijetD"):
            self.assertIn(name, star510)

    def test_hermes_kinematics_pid_cuts_and_every_grid(self) -> None:
        source = (
            DISPOL_ROOT / "analyses/rivet/dis/HERMES_2019_I1698889.cc"
        ).read_text(encoding="utf-8")
        for token in (
            "out.Q2 = -out.q.mass2()/GeV2",
            "out.x = -out.q.mass2()/(2.0*pdotq)",
            "out.y = pdotq/pdotk",
            "out.W2 = (out.target+out.q).mass2()/GeV2",
            "out.theta = out.k.angle(out.kp)",
            "out.z = (dis.target*momentum)/pdotq",
            "out.xF = 2.0*hcm.p3().dot",
            "dis.Q2 <= 1.0",
            "dis.W2 <= 10.0",
            "dis.y >= 0.85",
            "dis.theta <= 0.04",
            "dis.theta >= 0.22",
            "hadron.xF <= 0.10",
            "hadron.momentum > 4.0",
            "hadron.momentum < 13.8",
            "hadron.momentum > 2.0",
            "hadron.momentum < 15.0",
            "icx*7+iz",
            "icx*6+ipt",
            "ix3*9+iz3*3+ipt3",
        ):
            self.assertIn(token, source)
        self.assertIn(
            "{0.023,0.040,0.055,0.075,0.100,0.140,0.200,0.300,0.400,0.600}",
            source,
        )
        self.assertIn(
            "{0.023,0.040,0.055,0.075,0.140,0.600}", source
        )
        self.assertIn(
            "{0.20,0.32,0.44,0.56,0.68,0.80}", source
        )
        self.assertIn("std::sqrt(positive*negative)", source)
        self.assertIn("counts.moment[index]*counts.momentDenominator[index]", source)
        self.assertNotIn("0.926", source)


if __name__ == "__main__":
    unittest.main()
