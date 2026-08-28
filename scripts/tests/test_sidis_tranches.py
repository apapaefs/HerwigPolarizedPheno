#!/usr/bin/env python3
from __future__ import annotations

import argparse
import contextlib
import copy
import csv
import gzip
import hashlib
import importlib.util
import io
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import compass_sidis_postprocess as postprocess  # noqa: E402
import phenomenology_reference_data as reference  # noqa: E402
import run_experimental_campaign as experimental  # noqa: E402
import run_phenomenology_campaign as campaign  # noqa: E402
import sidis_tranche_reference_data as tranche_reference  # noqa: E402


IDS = (
    "COMPASS_2026_I3096394",
    "COMPASS_2025_I2840545",
    "COMPASS_2010_I862410",
    "HERMES_2013_I1208547",
    "COMPASS_2018_I1624692",
    "COMPASS_2020_I1788430",
)


def arguments(profile: str = "central", **updates: object) -> argparse.Namespace:
    values = {
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
        "nominal_prediction": None,
        "plot_comparisons": False,
    }
    values.update(updates)
    return argparse.Namespace(**values)


def series(values: list[float], variances: list[float] | None = None,
           edges: list[float] | None = None) -> experimental.BinSeries:
    return experimental.BinSeries(
        edges or [float(index) for index in range(len(values)+1)],
        values,
        variances or [0.0 for _ in values],
    )


def load_controller():
    path = ROOT / "campaigns/control/sidis-tranches-20260828/controller.py"
    spec = importlib.util.spec_from_file_location("sidis_tranche_controller", path)
    if spec is None or spec.loader is None:  # pragma: no cover
        raise RuntimeError("Could not load SIDIS tranche controller")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TrancheReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.snapshots = {
            identifier: tranche_reference.validate_vendored(
                identifier,
                full_covariance=identifier == "HERMES_2013_I1208547",
            )
            for identifier in IDS
        }

    def test_complete_sources_and_checksums(self) -> None:
        expected = {
            "COMPASS_2026_I3096394": ("10.17182/hepdata.169859.v1", 3),
            "COMPASS_2025_I2840545": ("10.17182/hepdata.159544.v1", 3),
            "HERMES_2013_I1208547": ("10.17182/hepdata.62097.v1", 64),
            "COMPASS_2018_I1624692": ("10.17182/hepdata.83542.v1", 162),
        }
        for identifier, (doi, count) in expected.items():
            manifest = json.loads((
                ROOT / f"data/phenomenology/{identifier}/source-manifest.json"
            ).read_text(encoding="utf-8"))
            self.assertEqual(manifest["record_doi"], doi)
            self.assertEqual(manifest["record_version"], 1)
            self.assertEqual(len(manifest["tables"]), count)
            for item in [
                {"path": manifest["record_path"], "sha256": manifest["record_sha256"]},
                *manifest["tables"],
            ]:
                source = ROOT / item["path"]
                self.assertEqual(
                    hashlib.sha256(source.read_bytes()).hexdigest(),
                    item["sha256"],
                )
        hermes = json.loads((
            ROOT / "data/phenomenology/HERMES_2013_I1208547/source-manifest.json"
        ).read_text(encoding="utf-8"))["full_archive"]
        archive = ROOT / hermes["path"]
        self.assertEqual(archive.stat().st_size, 25_937_678)
        self.assertEqual(
            hashlib.sha256(archive.read_bytes()).hexdigest(),
            "e54ae9e96fb417f43c7845e11319977c21b4c0f7349f00ca987658e00b181e1b",
        )
        compass = json.loads((
            ROOT / "data/phenomenology/COMPASS_2010_I862410/source-manifest.json"
        ).read_text(encoding="utf-8"))
        for key in ("source_archive", "paper_pdf"):
            source = ROOT / compass[key]["path"]
            self.assertEqual(source.stat().st_size, compass[key]["bytes"])
            self.assertEqual(hashlib.sha256(source.read_bytes()).hexdigest(),
                             compass[key]["sha256"])
        ratios = json.loads((
            ROOT / "data/phenomenology/COMPASS_2020_I1788430/source-manifest.json"
        ).read_text(encoding="utf-8"))
        self.assertEqual(ratios["hepdata_audit"]["record_http_status"], 404)
        self.assertEqual(
            ratios["hepdata_audit"]["result"],
            "no official HEPData submission found",
        )
        for key in (
            "source_archive", "paper_pdf", "inherited_kaon_selection_source",
        ):
            source = ROOT / ratios[key]["path"]
            self.assertEqual(source.stat().st_size, ratios[key]["bytes"])
            self.assertEqual(
                hashlib.sha256(source.read_bytes()).hexdigest(),
                ratios[key]["sha256"],
            )

    def test_exact_sparse_cell_maps_and_correction_provenance(self) -> None:
        expected = {
            "COMPASS_2026_I3096394": {
                "hplus": (306, "682ec9d8103269a7c0a8c6abbd376a0a970e87a9d1085f19207fe97515e03401"),
                "kplus": (305, "4f6d3318e9f17400fe07b2e5fe64f0555bbc89d92948f231385e22c7f9cd9833"),
            },
            "COMPASS_2025_I2840545": {
                "hplus": (302, "62ab8829c504fb6b510db48766a5aedaf1fbba4e3b8ac1694320d106b7e4ee69"),
                "kplus": (298, "6aadf44efdc977e472ca900c89b2c0fd76bcf1831b890289a4a5f8b40abc0b73"),
            },
            "COMPASS_2018_I1624692": {
                "hminus": (2332, "2c23642ed449b8959fe2a51add91c0660045a33d44f4873d604b2c429651dc81"),
                "hplus": (2332, "118d5e295d59350f76056c1258af4b153bfd5eeb5723b56385058a8dc94bfe45"),
            },
            "COMPASS_2020_I1788430": {
                "pbar_over_p_xz": (18, "8d26d1728ebd0b5da6d3003db51610d7462a9fb323b2a149f0ff51a14b84efd9"),
                "pbar_over_p_lowx_zp": (34, "0cbf7c6795d6b22438eca42d61b93ec737f589707ce2452c7c664c0544cb2e49"),
                "kminus_over_kplus_lowx_zp": (15, "594bfd350080835b915e1fc708008f1c507fc4540579c042290a6cbc4770e562"),
            },
        }
        for identifier, datasets in expected.items():
            for name, (count, digest) in datasets.items():
                dataset = self.snapshots[identifier]["datasets"][name]
                self.assertEqual(len(dataset["points"]), count)
                payload = json.dumps(
                    dataset["valid_cell_map"], separators=(",", ":")
                ).encode()
                self.assertEqual(hashlib.sha256(payload).hexdigest(), digest)
                for point in dataset["points"]:
                    self.assertGreater(point["z_high"]-point["z_low"], 0.0)
                    if identifier != "COMPASS_2020_I1788430":
                        self.assertTrue(point["corrections"])
        pt = self.snapshots["COMPASS_2018_I1624692"]
        self.assertEqual(pt["source_count_audit"]["hepdata_v1_rows"], 4664)
        self.assertEqual(pt["source_count_audit"]["assessment_claim"], 4918)
        self.assertIn("never invent", pt["source_count_audit"]["disposition"])

    def test_systematic_splits_and_compass_2010_paper_audit(self) -> None:
        for identifier in ("COMPASS_2026_I3096394", "COMPASS_2025_I2840545"):
            for dataset in self.snapshots[identifier]["datasets"].values():
                for point in dataset["points"]:
                    self.assertAlmostEqual(
                        point["systematic_correlated_80pct"],
                        .8*point["systematic"],
                    )
                    self.assertAlmostEqual(
                        point["systematic_uncorrelated_60pct"],
                        .6*point["systematic"],
                    )
        snapshot = self.snapshots["COMPASS_2010_I862410"]
        self.assertEqual(len(snapshot["statistical_correlation_blocks"]), 12)
        self.assertEqual(len(snapshot["statistical_covariance"]), 48)
        for dataset in snapshot["datasets"].values():
            for point in dataset["points"]:
                correlated = .06*float(point["a1"])
                self.assertAlmostEqual(
                    point["systematic_correlated_6pct"], correlated
                )
                self.assertAlmostEqual(
                    point["systematic"]**2,
                    correlated**2 + point["systematic_uncorrelated"]**2,
                )
        audit = json.loads((
            ROOT / "data/phenomenology/COMPASS_2010_I862410/paper-extraction-audit.json"
        ).read_text(encoding="utf-8"))
        self.assertEqual(len(audit["asymmetry_rows"]), 12)
        self.assertEqual(len(audit["correlation_rows"]), 6)
        self.assertTrue(all(len(row) == 12 for row in audit["correlation_rows"]))
        self.assertEqual(len(audit["deuteron_correction_rows"]), 12)
        self.assertEqual(audit["missing_rows"], [])

        ratios = self.snapshots["COMPASS_2020_I1788430"]
        self.assertEqual(sum(
            len(dataset["points"]) for dataset in ratios["datasets"].values()
        ), 67)
        for dataset in ratios["datasets"].values():
            for point in dataset["points"]:
                self.assertAlmostEqual(
                    point["systematic"]**2,
                    point["systematic_correlated_sqrt75"]**2
                    + point["systematic_uncorrelated_half"]**2,
                )
        ratio_audit = json.loads((
            ROOT / "data/phenomenology/COMPASS_2020_I1788430/"
            "paper-extraction-audit.json"
        ).read_text(encoding="utf-8"))
        self.assertEqual(
            [ratio_audit["tables"][key]["row_count"] for key in (
                "pbar_over_p_xz", "pbar_over_p_lowx_zp",
                "kminus_over_kplus_lowx_zp",
            )],
            [18,34,15],
        )
        self.assertEqual(ratio_audit["missing_rows"], [])

    def test_hermes_full_archive_covariance_and_projection_contract(self) -> None:
        snapshot = reference.validate_vendored(
            "HERMES_2013_I1208547", full_covariance=True
        )
        self.assertEqual(len(snapshot["datasets"]), 40)
        self.assertEqual(sum(len(item["points"]) for item in snapshot["datasets"].values()), 6592)
        expected_sizes = {
            "z-3D": 80, "zpt-3D": 90, "zx-3D": 180,
            "zQ2-3D": 180, "zxpt-3D": 294,
        }
        for dataset in snapshot["datasets"].values():
            size = expected_sizes[dataset["binning"]]
            covariance = dataset["statistical_covariance"]
            self.assertEqual(len(covariance), size)
            self.assertTrue(all(len(row) == size for row in covariance))
            for index, point in enumerate(dataset["points"]):
                self.assertAlmostEqual(covariance[index][index], point["stat"]**2,
                                       delta=2.e-12+2.e-5*point["stat"]**2)
                self.assertGreaterEqual(point["z_low"], .2)
        self.assertEqual(len(snapshot["readable_projections"]), 104)
        for projection in snapshot["readable_projections"].values():
            self.assertTrue(projection["integrated_projection"])
            self.assertIn("ProjectionHadronNumerator_",
                          projection["raw_objects"]["numerator"])
            self.assertIn("ProjectionDISDenominator_",
                          projection["raw_objects"]["denominator"])
            self.assertEqual(len(projection["points"])+1, len(projection["edges"]))
        audit = json.loads((
            ROOT / "data/phenomenology/HERMES_2013_I1208547/"
            "archive-hepdata-projection-audit.json"
        ).read_text(encoding="utf-8"))
        self.assertEqual(
            [item["table"] for item in audit["hepdata_tables"]],
            list(range(1, 65)),
        )
        self.assertEqual(audit["numeric_comparison_count"], 4544)
        self.assertEqual(audit["discrepancy_count"], 0)
        self.assertTrue(any(
            item["raw_asymmetric_pair_count"]
            for item in snapshot["covariance_normalization"]["datasets"].values()
        ))
        self.assertEqual(sum(
            item["cauchy_schwarz_replacement_count"]
            for item in snapshot["covariance_normalization"]["datasets"].values()
        ), 9)

    def test_reference_yoda_inventory_includes_integrated_projections(self) -> None:
        identifier = "HERMES_2013_I1208547"
        snapshot = self.snapshots[identifier]
        with gzip.open(
            ROOT / f"analyses/rivet/dis/{identifier}.yoda.gz",
            "rt", encoding="utf-8",
        ) as stream:
            paths = {
                line.split(maxsplit=2)[2].strip()
                for line in stream
                if line.startswith("BEGIN YODA_ESTIMATE1D_V3 ")
            }
        expected = {
            "/REF" + str(dataset["flat_rivet_path"])
            for dataset in snapshot["datasets"].values()
        }
        expected.update(
            "/REF" + str(item["rivet_path"])
            for dataset in snapshot["datasets"].values()
            for item in dataset["slices"]
        )
        expected.update(
            "/REF" + str(item["rivet_path"])
            for item in snapshot["readable_projections"].values()
        )
        self.assertEqual(paths, expected)

    def test_compass_2020_reference_yoda_inventory_is_complete(self) -> None:
        identifier = "COMPASS_2020_I1788430"
        snapshot = self.snapshots[identifier]
        expected = {
            "/REF" + str(dataset["flat_rivet_path"])
            for dataset in snapshot["datasets"].values()
        }
        expected.update(
            "/REF" + str(item["rivet_path"])
            for dataset in snapshot["datasets"].values()
            for item in dataset["slices"]
        )
        self.assertEqual(len(expected), 19)
        # Reference YODA is a generated product under the workspace contract.
        # When present (after fetch-data/prepare), it must exactly mirror the
        # authoritative normalized JSON; a clean clone need not contain it.
        reference_yoda = ROOT / f"analyses/rivet/dis/{identifier}.yoda.gz"
        if reference_yoda.is_file():
            with gzip.open(reference_yoda, "rt", encoding="utf-8") as stream:
                paths = {
                    line.split(maxsplit=2)[2].strip()
                    for line in stream
                    if line.startswith("BEGIN YODA_ESTIMATE1D_V3 ")
                }
            self.assertEqual(paths, expected)


class TrancheRunnerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.registry = campaign.discover_pp_registry()

    def test_schema6_targets_active_axes_and_exact_job_counts(self) -> None:
        central = dict(zip(IDS, (4, 2, 8, 4, 4, 4)))
        paper = dict(zip(IDS, (412, 206, 1624, 412, 412, 412)))
        for identifier in IDS:
            measurement = self.registry[identifier]
            self.assertEqual(measurement["schema_version"], 6)
            self.assertTrue(measurement["postprocess_config"]["target_outputs"])
            central_jobs = campaign.build_job_matrix(
                measurement, campaign._resolved_options(arguments(), measurement)
            )
            paper_jobs = campaign.build_job_matrix(
                measurement,
                campaign._resolved_options(arguments("paper"), measurement),
            )
            self.assertEqual(len(central_jobs), central[identifier])
            self.assertEqual(len(paper_jobs), paper[identifier])
            if measurement["process_kind"] == "unpolarized_sidis":
                self.assertFalse(measurement["pdf_ensembles"]["polarized"]["active"])
                self.assertTrue(measurement["pdf_ensembles"]["unpolarized"]["active"])
                self.assertEqual(set(measurement["cards"]["helicities"]), {"00"})
                with self.assertRaisesRegex(campaign.CampaignError, "central polarized-PDF"):
                    campaign._resolved_options(
                        arguments(polarized_pdf_members="1"), measurement
                    )
            else:
                self.assertEqual(
                    set(measurement["cards"]["helicities"]),
                    set(postprocess.HELICITIES),
                )
                self.assertTrue(measurement["pdf_ensembles"]["polarized"]["active"])
                self.assertEqual(
                    measurement["postprocess_config"]["longitudinal_target_scale"],
                    1.0,
                )

    def test_schema6_rejects_invalid_target_components(self) -> None:
        measurement = copy.deepcopy(self.registry["COMPASS_2026_I3096394"])
        measurement["cards"]["target_components"] = {}
        with self.assertRaises(campaign.CampaignError):
            campaign._validate_descriptor(measurement, Path("synthetic.json"))
        measurement = copy.deepcopy(self.registry["COMPASS_2026_I3096394"])
        measurement["postprocess_config"]["target_outputs"] = {
            "D": {"P": .5, "X": .5}
        }
        with self.assertRaises(campaign.CampaignError):
            campaign._validate_descriptor(measurement, Path("synthetic.json"))

    def test_proton_and_isoscalar_target_order_density_and_masking(self) -> None:
        numerator = {"P": series([4.0], [1.0]), "N": series([2.0], [1.0])}
        denominator = {"P": series([10.0], [1.0]), "N": series([6.0], [1.0])}
        covariance = {"P": series([0.0], [.2]), "N": series([0.0], [.1])}
        proton = postprocess.multiplicity_target_combination(
            {"P": numerator["P"]}, {"P": denominator["P"]},
            {"P": covariance["P"]}, {"P": 1.0}, [.1],
        )
        deuteron = postprocess.multiplicity_target_combination(
            numerator, denominator, covariance, {"P": .5, "N": .5}, [.1],
        )
        self.assertAlmostEqual(proton["values"][0], 4.0)
        self.assertAlmostEqual(deuteron["values"][0], 3.75)
        swapped = postprocess.multiplicity_target_combination(
            {"N": numerator["N"], "P": numerator["P"]},
            {"N": denominator["N"], "P": denominator["P"]},
            {"N": covariance["N"], "P": covariance["P"]},
            {"N": .5, "P": .5}, [.1],
        )
        self.assertEqual(swapped["values"], deuteron["values"])
        masked = postprocess.multiplicity_target_combination(
            {"P": series([1.0])}, {"P": series([0.0])},
            {"P": series([0.0])}, {"P": 1.0}, [.1],
        )
        self.assertEqual(masked["values"], [None])

    def test_charge_ratio_target_order_covariance_and_tablewise_fit(self) -> None:
        negative = {"P": series([4.0], [1.0]), "N": series([2.0], [1.0])}
        positive = {"P": series([10.0], [4.0]), "N": series([6.0], [1.0])}
        covariance = {"P": series([0.0], [.5]), "N": series([0.0], [.25])}
        ratio = postprocess.charge_ratio_target_combination(
            negative, positive, covariance, {"P": .5, "N": .5}
        )
        self.assertAlmostEqual(ratio["values"][0], 3.0/8.0)
        self.assertGreater(ratio["errors"][0], 0.0)
        masked = postprocess.charge_ratio_target_combination(
            {"P": series([1.0])}, {"P": series([0.0])},
            {"P": series([0.0])}, {"P": 1.0},
        )
        self.assertEqual(masked["values"], [None])

        snapshot = reference.validate_vendored("COMPASS_2020_I1788430")
        predictions = {
            dataset_id: {
                "values": [point["value"] for point in dataset["points"]],
                "errors": [0.01 for _ in dataset["points"]],
            }
            for dataset_id, dataset in snapshot["datasets"].items()
        }
        goodness = postprocess.charge_ratio_goodness_of_fit(
            predictions, snapshot
        )
        self.assertEqual(
            set(goodness["independent_datasets"]), set(snapshot["datasets"])
        )
        self.assertTrue(all(
            result["chi2_correlated"] == 0.0
            for result in goodness["independent_datasets"].values()
        ))

    def test_hermes_projection_reference_lookup(self) -> None:
        measurement = self.registry["HERMES_2013_I1208547"]
        snapshot = reference.validate_vendored("HERMES_2013_I1208547")
        observable, projection = next(iter(snapshot["readable_projections"].items()))
        self.assertEqual(
            campaign._reference_path(measurement, observable, snapshot),
            projection["rivet_path"],
        )
        points = campaign._reference_points(measurement, snapshot, observable)
        self.assertEqual(points, projection["points"])
        source = (ROOT / "analyses/rivet/dis/HERMES_2013_I1208547.cc").read_text()
        self.assertIn("ProjectionHadronNumerator_", source)
        self.assertIn("ProjectionDISDenominator_", source)
        self.assertIn("count += counts", source)

    def test_modern_compass_slice_overlay_uses_bin_midpoint_without_means(self) -> None:
        measurement = self.registry["COMPASS_2025_I2840545"]
        snapshot = reference.validate_vendored("COMPASS_2025_I2840545")
        dataset = snapshot["datasets"]["hplus"]
        slice_spec = dataset["slices"][0]
        points = campaign._pp_reference_overlay_points(
            measurement, snapshot, Path(slice_spec["rivet_path"]).name
        )
        self.assertIsNotNone(points)
        first = dataset["points"][slice_spec["flat_bins"][0]]
        self.assertAlmostEqual(
            points[0]["plot_x"], .5*(first["z_low"]+first["z_high"])
        )


class TrancheControllerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.controller = load_controller()

    def test_frozen_totals_and_dry_run(self) -> None:
        config = self.controller.configuration()
        self.assertEqual(config["expected"]["pilot_events"], 64_900_000)
        self.assertEqual(config["expected"]["central_floor_events"], 618_200_000)
        self.assertEqual(config["expected"]["central_floor_shards"], 11_800)
        with contextlib.redirect_stdout(io.StringIO()) as output:
            self.controller.action_dry_run("paper", "all")
        payload = json.loads(output.getvalue())
        self.assertEqual(
            sum(item["logical_jobs"] for item in payload["campaigns"]), 3478
        )
        self.assertEqual(
            sum(item["total_events"] for item in payload["campaigns"]),
            64_994_600_000,
        )

    def test_event_plans_are_checksum_pinned_and_immutable(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "event-plan.json"
            payload = {"schema_version": 1, "entries": [{"events": 100}]}
            digest = self.controller.immutable_json(path, payload)
            self.assertEqual(digest, hashlib.sha256(
                self.controller.canonical_json(payload)
            ).hexdigest())
            self.assertEqual(self.controller.verify_immutable_json(path), payload)
            with self.assertRaises(self.controller.ControllerError):
                self.controller.immutable_json(
                    path, {"schema_version": 1, "entries": [{"events": 200}]}
                )

    def test_assessor_does_not_double_count_masked_bins(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            (directory / "postprocess").mkdir()
            (directory / "manifest.json").write_text("{}\n")
            (directory / "postprocess/summary.json").write_text(
                json.dumps({"masked_bins": {"primary": [0]}}) + "\n"
            )
            with (directory / "postprocess/central.csv").open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=(
                    "theory", "mc_stat", "data_stat", "data_systematic"
                ))
                writer.writeheader()
                writer.writerow({
                    "theory": "", "mc_stat": "", "data_stat": ".1",
                    "data_systematic": ".2",
                })
            original = self.controller.campaign_directory
            self.controller.campaign_directory = lambda _spec: directory
            try:
                result = self.controller.assess_one({"measurement": "synthetic"})
            finally:
                self.controller.campaign_directory = original
            self.assertEqual(result["masked_primary_bins"], 1)
            self.assertEqual(result["summary_masked_primary_bins"], 1)
            self.assertEqual(result["csv_masked_primary_bins"], 1)
            self.assertFalse(result["current_sample_passes"])

    def test_assessor_refuses_zero_information_extrapolation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            (directory / "postprocess").mkdir()
            (directory / "manifest.json").write_text("{}\n")
            (directory / "postprocess/summary.json").write_text(
                json.dumps({"masked_bins": {}}) + "\n"
            )
            with (directory / "postprocess/central.csv").open(
                "w", newline=""
            ) as stream:
                writer = csv.DictWriter(stream, fieldnames=(
                    "theory", "mc_stat", "data", "data_stat",
                    "data_systematic",
                ))
                writer.writeheader()
                writer.writerow({
                    "theory": "0", "mc_stat": "0", "data": ".2",
                    "data_stat": ".01", "data_systematic": ".02",
                })
            original = self.controller.campaign_directory
            self.controller.campaign_directory = lambda _spec: directory
            try:
                result = self.controller.assess_one({"measurement": "synthetic"})
            finally:
                self.controller.campaign_directory = original
            self.assertEqual(result["masked_primary_bins"], 1)
            self.assertEqual(result["zero_information_primary_bins"], 1)
            self.assertFalse(result["current_sample_passes"])


if __name__ == "__main__":
    unittest.main()
