#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import compass_sidis_postprocess as postprocess  # noqa: E402
import phenomenology_reference_data as reference  # noqa: E402
import run_experimental_campaign as experimental  # noqa: E402
import run_phenomenology_campaign as campaign  # noqa: E402


IDS = (
    "COMPASS_2009_I820721",
    "COMPASS_2017_I1444985",
    "COMPASS_2017_I1483098",
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


def series(
    values: list[float], variances: list[float] | None = None
) -> experimental.BinSeries:
    return experimental.BinSeries(
        [float(index) for index in range(len(values) + 1)],
        values,
        variances or [0.0 for _ in values],
    )


class ReferenceDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.snapshots = {
            identifier: reference.validate_vendored(identifier)
            for identifier in IDS
        }

    def test_complete_v1_inventories_and_checksums(self) -> None:
        expected = {
            "COMPASS_2009_I820721": ("10.17182/hepdata.55300.v1", 4),
            "COMPASS_2017_I1444985": ("10.17182/hepdata.76800.v1", 4),
            "COMPASS_2017_I1483098": ("10.17182/hepdata.77892.v1", 2),
        }
        for identifier, (doi, count) in expected.items():
            root = ROOT / "data" / "phenomenology" / identifier
            manifest = json.loads(
                (root / "source-manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(manifest["record_doi"], doi)
            self.assertEqual(manifest["record_version"], 1)
            self.assertEqual(len(manifest["tables"]), count)
            record = ROOT / manifest["record_path"]
            self.assertEqual(
                hashlib.sha256(record.read_bytes()).hexdigest(),
                manifest["record_sha256"],
            )
            for table in manifest["tables"]:
                source = ROOT / table["path"]
                self.assertEqual(
                    hashlib.sha256(source.read_bytes()).hexdigest(),
                    table["sha256"],
                )

    def test_2009_bins_covariance_systematics_and_frozen_audit(self) -> None:
        snapshot = self.snapshots["COMPASS_2009_I820721"]
        self.assertEqual(
            snapshot["binning"]["x"],
            [0.004, 0.006, 0.01, 0.02, 0.03, 0.04,
             0.06, 0.1, 0.15, 0.2, 0.3],
        )
        self.assertEqual(len(snapshot["point_order"]), 40)
        self.assertEqual(len(snapshot["statistical_correlation_blocks"]), 10)
        try:
            import numpy as np
        except ImportError as exc:  # pragma: no cover - required by validator
            self.fail(str(exc))
        covariance = np.asarray(snapshot["statistical_covariance"])
        self.assertTrue(np.allclose(covariance, covariance.T, atol=1.0e-12))
        self.assertGreaterEqual(float(np.linalg.eigvalsh(covariance)[0]), -1.0e-12)
        for dataset in snapshot["datasets"].values():
            self.assertTrue(str(dataset["rivet_path"]).startswith("/COMPASS_2009_"))
            for point in dataset["points"]:
                correlated = 0.08 * float(point["a1"])
                self.assertAlmostEqual(
                    point["systematic_correlated_8pct"], correlated
                )
                self.assertAlmostEqual(
                    point["systematic"] ** 2,
                    correlated**2 + point["systematic_uncorrelated"] ** 2,
                )
        audit = json.loads(
            (
                ROOT
                / "data/phenomenology/COMPASS_2009_I820721/"
                "paper-hepdata-discrepancy-audit.json"
            ).read_text(encoding="utf-8")
        )
        self.assertEqual(audit["numerical_authority"], "HEPData v1")
        self.assertEqual(audit["discrepancy_count"], 88)
        self.assertEqual(len(audit["discrepancies"]), 88)
        self.assertEqual(
            audit["paper_source"]["archive_sha256"],
            "ea66f1b938ed0976b14212867b77181c7bd5c38a1bb324e177a1b72aaa00a7b8",
        )

    def test_2017_sparse_cells_corrections_and_systematic_split(self) -> None:
        expected = {
            "COMPASS_2017_I1444985": (
                311, 4,
                "a6ea9d4f1151c727b02f5e2751ffdef81fd65f59407c3c78958a283c8c66deab",
            ),
            "COMPASS_2017_I1483098": (
                309, 2,
                "62eb2206c4e92bf7354002448b39a7d8164e2714497f2172775209423c2bbba4",
            ),
        }
        full_cartesian = 9 * 5 * 12
        expected_binning = {
            "x": [0.004, 0.01, 0.02, 0.03, 0.04, 0.06, 0.1, 0.14, 0.18, 0.4],
            "y": [0.1, 0.15, 0.2, 0.3, 0.5, 0.7],
            "z": [0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55,
                  0.6, 0.65, 0.7, 0.75, 0.85],
        }
        for identifier, (cells, species_count, cell_map_sha256) in expected.items():
            snapshot = self.snapshots[identifier]
            self.assertEqual(snapshot["binning"], expected_binning)
            self.assertEqual(snapshot["cell_count"], cells)
            self.assertLess(cells, full_cartesian)
            self.assertEqual(snapshot["slice_count"], 38)
            self.assertEqual(len(snapshot["valid_cell_map"]), cells)
            cell_map = json.dumps(
                snapshot["valid_cell_map"], separators=(",", ":")
            ).encode()
            self.assertEqual(
                hashlib.sha256(cell_map).hexdigest(), cell_map_sha256
            )
            self.assertEqual(len(snapshot["datasets"]), species_count)
            signatures = []
            for dataset in snapshot["datasets"].values():
                self.assertEqual(len(dataset["points"]), cells)
                self.assertEqual(len(dataset["slices"]), 38)
                self.assertTrue(
                    str(dataset["flat_rivet_path"]).endswith("_cells")
                )
                signatures.append(
                    [
                        (
                            point["x_low"], point["x_high"],
                            point["y_low"], point["y_high"],
                            point["z_low"], point["z_high"],
                        )
                        for point in dataset["points"]
                    ]
                )
                for point in dataset["points"]:
                    self.assertEqual(
                        set(point["corrections"]),
                        {"radiative_hadron", "radiative_dis", "dvm_hadron", "dvm_dis"},
                    )
                    self.assertAlmostEqual(
                        point["systematic_correlated_80pct"],
                        0.8 * point["systematic"],
                    )
                    self.assertAlmostEqual(
                        point["systematic_uncorrelated_60pct"],
                        0.6 * point["systematic"],
                    )
            self.assertTrue(all(item == signatures[0] for item in signatures[1:]))

    def test_reference_yoda_has_complete_flat_and_slice_inventory(self) -> None:
        expected_counts = {
            "COMPASS_2009_I820721": 4,
            "COMPASS_2017_I1444985": 4 * (1 + 38),
            "COMPASS_2017_I1483098": 2 * (1 + 38),
        }
        for identifier, expected_count in expected_counts.items():
            path = ROOT / "analyses" / "rivet" / "dis" / f"{identifier}.yoda.gz"
            with gzip.open(path, "rt", encoding="utf-8") as stream:
                paths = {
                    line.split(maxsplit=2)[2].strip()
                    for line in stream
                    if line.startswith("BEGIN YODA_ESTIMATE1D_V3 ")
                }
            snapshot = self.snapshots[identifier]
            expected_paths: set[str] = set()
            for dataset in snapshot["datasets"].values():
                primary = dataset.get("rivet_path") or dataset.get(
                    "flat_rivet_path"
                )
                expected_paths.add("/REF" + str(primary))
                expected_paths.update(
                    "/REF" + str(item["rivet_path"])
                    for item in dataset.get("slices", [])
                )
            self.assertEqual(len(paths), expected_count)
            self.assertEqual(paths, expected_paths)


class EstimatorTests(unittest.TestCase):
    def test_signed_nlo_addition_precedes_a1_and_retains_covariance(self) -> None:
        combined = postprocess.experimental.add_independent_series(
            [series([3.0], [0.4]), series([-1.0], [0.1])]
        )
        self.assertEqual(combined.values, [2.0])
        self.assertEqual(combined.variances, [0.5])

        ordinary_values = {
            "P:PP": 12.0, "P:PM": 8.0, "P:MP": 8.0, "P:MM": 12.0,
            "N:PP": 10.0, "N:PM": 6.0, "N:MP": 6.0, "N:MM": 10.0,
        }
        ordinary = {
            label: series([value], [1.0])
            for label, value in ordinary_values.items()
        }
        inverse_d = {
            label: series([2.0 * value], [4.0])
            for label, value in ordinary_values.items()
        }
        covariance = {
            label: series([0.0], [2.0]) for label in ordinary_values
        }
        result = postprocess.a1_deuteron(ordinary, inverse_d, covariance)
        self.assertAlmostEqual(result["values"][0], 3.7 / 9.0)
        self.assertGreater(result["errors"][0], 0.0)
        self.assertAlmostEqual(
            result["sigma_ll_over_d"].values[0], 3.7
        )
        self.assertAlmostEqual(result["sigma_uu"].values[0], 9.0)
        self.assertAlmostEqual(
            result["numerator_denominator_covariance"][0], 0.0
        )

    def test_a1_masks_nonpositive_denominator(self) -> None:
        labels = {
            f"{target}:{helicity}"
            for target in ("P", "N")
            for helicity in postprocess.HELICITIES
        }
        samples = {label: series([0.0]) for label in labels}
        result = postprocess.a1_deuteron(samples, samples, samples)
        self.assertEqual(result["values"], [None])
        self.assertEqual(result["errors"], [None])

    def test_isoscalar_sum_ratio_width_covariance_and_pn_order(self) -> None:
        numerators = {
            "P": series([4.0], [4.0]),
            "N": series([2.0], [1.0]),
        }
        denominators = {
            "P": series([10.0], [9.0]),
            "N": series([6.0], [4.0]),
        }
        covariance = {
            "P": series([0.0], [2.0]),
            "N": series([0.0], [1.0]),
        }
        result = postprocess.multiplicity_isoscalar(
            numerators, denominators, covariance, [0.05]
        )
        self.assertAlmostEqual(result["values"][0], 7.5)
        ratio_variance = 1.25 / 8.0**2 + 3.0**2 * 3.25 / 8.0**4 - 2 * 3.0 * 0.75 / 8.0**3
        self.assertAlmostEqual(
            result["errors"][0], math.sqrt(ratio_variance) / 0.05
        )
        swapped = postprocess.multiplicity_isoscalar(
            {"P": numerators["N"], "N": numerators["P"]},
            {"P": denominators["N"], "N": denominators["P"]},
            {"P": covariance["N"], "N": covariance["P"]},
            [0.05],
        )
        self.assertEqual(swapped["values"], result["values"])
        masked = postprocess.multiplicity_isoscalar(
            {"P": series([1.0]), "N": series([1.0])},
            {"P": series([-1.0]), "N": series([0.0])},
            {"P": series([0.0]), "N": series([0.0])},
            [0.1],
        )
        self.assertEqual(masked["values"], [None])

    def test_correlated_goodness_models(self) -> None:
        snapshot = reference.validate_vendored("COMPASS_2009_I820721")
        prediction = {
            species: {
                "values": [float(point["a1"]) for point in dataset["points"]],
                "errors": [0.0] * len(dataset["points"]),
            }
            for species, dataset in snapshot["datasets"].items()
        }
        result = postprocess.a1_goodness_of_fit(prediction, snapshot)
        self.assertEqual(result["points"], 40)
        self.assertAlmostEqual(result["chi2_correlated"], 0.0)

        small_snapshot = {
            "datasets": {
                "a": {"points": [{
                    "value": 1.0, "stat": 0.1,
                    "systematic_uncorrelated_60pct": 0.12,
                    "systematic_correlated_80pct": 0.16,
                }]},
                "b": {"points": [{
                    "value": 2.0, "stat": 0.2,
                    "systematic_uncorrelated_60pct": 0.18,
                    "systematic_correlated_80pct": 0.24,
                }]},
            }
        }
        equal = {
            "a": {"values": [1.0], "errors": [0.0]},
            "b": {"values": [2.0], "errors": [0.0]},
        }
        multiplicity = postprocess.multiplicity_goodness_of_fit(
            equal, small_snapshot
        )
        self.assertEqual(multiplicity["points"], 2)
        self.assertAlmostEqual(multiplicity["chi2_correlated"], 0.0)
        self.assertIn("record-wide 0.8", multiplicity["covariance"])

    def test_rank_one_covariance_matches_dense_inverse(self) -> None:
        try:
            import numpy as np
        except ImportError as exc:  # pragma: no cover - required by validator
            self.fail(str(exc))
        residual = [0.4, -0.2, 0.7]
        diagonal = [0.09, 0.16, 0.25]
        nuisance = [0.1, 0.2, -0.05]
        fast = postprocess._diagonal_plus_rank_one_result(
            residual, diagonal, nuisance, ["a", "b", "c"], "test"
        )
        covariance = np.diag(diagonal) + np.outer(nuisance, nuisance)
        expected = float(
            np.asarray(residual)
            @ np.linalg.inv(covariance)
            @ np.asarray(residual)
        )
        self.assertAlmostEqual(fast["chi2_correlated"], expected)
        self.assertAlmostEqual(
            sum(value**2 for value in fast["decorrelated_pulls"]),
            expected,
        )


class RegistryRunnerAndSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.registry = campaign.discover_pp_registry()

    def test_schema_v5_process_helicities_and_active_axes(self) -> None:
        polarized = self.registry["COMPASS_2009_I820721"]
        self.assertEqual(polarized["schema_version"], 5)
        self.assertEqual(polarized["process_kind"], "polarized_sidis")
        self.assertEqual(set(polarized["cards"]["helicities"]), set(postprocess.HELICITIES))
        for identifier in IDS[1:]:
            unpolarized = self.registry[identifier]
            self.assertEqual(unpolarized["schema_version"], 5)
            self.assertEqual(unpolarized["process_kind"], "unpolarized_sidis")
            self.assertEqual(set(unpolarized["cards"]["helicities"]), {"00"})
            self.assertFalse(unpolarized["pdf_ensembles"]["polarized"]["active"])
            self.assertTrue(unpolarized["pdf_ensembles"]["unpolarized"]["active"])

    def test_exact_central_jobs_and_default_statistics(self) -> None:
        expected = {
            "COMPASS_2009_I820721": 16,
            "COMPASS_2017_I1444985": 4,
            "COMPASS_2017_I1483098": 4,
        }
        for identifier, count in expected.items():
            measurement = self.registry[identifier]
            options = campaign._resolved_options(arguments(), measurement)
            jobs = campaign.build_job_matrix(measurement, options)
            self.assertEqual(len(jobs), count)
            self.assertEqual(options["events_by_contribution"], {
                "POSNLO": 300000, "NEGNLO": 30000,
            })
            self.assertEqual(measurement["campaign"]["smoke_events"], 100)
            self.assertEqual(options["shards"], 1)
            self.assertEqual(len({job["id"] for job in jobs}), count)
            self.assertEqual(len({job["seed"] for job in jobs}), count)

    def test_unpolarized_selector_rejection_and_paper_axis(self) -> None:
        measurement = self.registry["COMPASS_2017_I1444985"]
        with self.assertRaisesRegex(
            campaign.CampaignError, "central polarized-PDF"
        ):
            campaign._resolved_options(
                arguments(polarized_pdf_members="1"), measurement
            )
        options = campaign._resolved_options(arguments("paper"), measurement)
        points = [tuple(point) for point in options["variation_points"]]
        self.assertEqual(len(points), 103)
        self.assertTrue(all(point[0] == 0 for point in points))
        self.assertEqual(
            {point[1] for point in points if point[2] == 1.0}, set(range(101))
        )

    def test_explicit_00_cards_and_source_contracts(self) -> None:
        for identifier in IDS[1:]:
            measurement = self.registry[identifier]
            jobs = campaign.build_job_matrix(
                measurement, campaign._resolved_options(arguments(), measurement)
            )
            card = campaign._card_text(measurement, jobs[0])
            self.assertIn("FirstLongitudinalPolarization 0", card)
            self.assertIn("SecondLongitudinalPolarization 0", card)
            self.assertIn("PowhegMEDISNCPol", card)
            self.assertIn("COMPASSSIDISDiffPDF:Member 0", card)

        source_2009 = (
            ROOT / "analyses/rivet/dis/COMPASS_2009_I820721.cc"
        ).read_text(encoding="utf-8")
        self.assertIn("YieldOverD_", source_2009)
        self.assertIn("entry.second/std::sqrt(D)", source_2009)
        self.assertIn("dis.Q2 > 1.0", source_2009)
        self.assertIn("hadron.momentum > 10.0", source_2009)

        for identifier in IDS[1:]:
            source = (
                ROOT / f"analyses/rivet/dis/{identifier}.cc"
            ).read_text(encoding="utf-8")
            self.assertIn("dis.W2 > 25.0", source)
            self.assertIn("multiplicityDISCell", source)
            self.assertIn("findCell", source)
            self.assertIn("CovarianceProxy_", source)
        charged = (
            ROOT / "analyses/rivet/dis/COMPASS_2017_I1444985.cc"
        ).read_text(encoding="utf-8")
        self.assertIn("hadronKinematics(dis, particle, true)", charged)
        self.assertIn("PID::isHadron", charged)
        helper = (
            ROOT / "analyses/rivet/dis/COMPASSSIDIS.hh"
        ).read_text(encoding="utf-8")
        self.assertIn("bool pionMassAssumption", helper)

    def test_generated_sparse_header_and_reference_paths(self) -> None:
        header = (
            ROOT / "analyses/rivet/dis/COMPASSSIDISBinning.hh"
        ).read_text(encoding="utf-8")
        self.assertIn("pionHadronCells()", header)
        self.assertIn("kaonCells()", header)
        for identifier in IDS:
            measurement = self.registry[identifier]
            snapshot = reference.validate_vendored(identifier)
            for species, dataset in snapshot["datasets"].items():
                path = campaign._reference_path(
                    measurement, species, snapshot
                )
                self.assertEqual(
                    path,
                    dataset.get("rivet_path") or dataset.get("flat_rivet_path"),
                )


if __name__ == "__main__":
    unittest.main()
