#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import sys
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import compass_sidis_postprocess as postprocess  # noqa: E402
import phenomenology_reference_data as reference  # noqa: E402
import run_experimental_campaign as experimental  # noqa: E402
import run_phenomenology_campaign as campaign  # noqa: E402
import sidis_diagnostic_reference_data as diagnostic_reference  # noqa: E402


EXTERNAL_IDS = ("COMPASS_2013_I1236358", "COMPASS_2014_I1278730")
ALL_IDS = (*EXTERNAL_IDS, "HERMES_2013_I1111237")


def series(values: list[float], variances: list[float]) -> experimental.BinSeries:
    return experimental.BinSeries(
        [float(index) for index in range(len(values) + 1)], values, variances
    )


def arguments(**updates: object) -> argparse.Namespace:
    values = {
        "profile": "central", "tag": "unit", "jobs": 4, "shards": 1,
        "seed_base": None, "lo_events": None,
        "posnlo_events": None, "negnlo_events": None,
        "smoke": False, "polarized_pdf_members": None,
        "unpolarized_pdf_members": None, "scales": None,
        "families": None, "jet_kt_min_gev": None,
        "nominal_prediction": None, "plot_comparisons": False,
    }
    values.update(updates)
    return argparse.Namespace(**values)


class DiagnosticReferenceTests(unittest.TestCase):
    def test_paper_sources_reconstruct_complete_snapshots(self) -> None:
        snapshots = {
            identifier: reference.validate_vendored(identifier)
            for identifier in EXTERNAL_IDS
        }
        slopes = snapshots["COMPASS_2013_I1236358"]
        self.assertEqual(len(slopes["datasets"]), 2)
        self.assertEqual(
            sum(len(item["points"]) for item in slopes["datasets"].values()),
            368,
        )
        self.assertTrue(all(
            len(item["slices"]) == 23 for item in slopes["datasets"].values()
        ))
        moments = snapshots["COMPASS_2014_I1278730"]
        self.assertEqual(len(moments["datasets"]), 16)
        self.assertEqual(
            sum(len(item["points"]) for item in moments["datasets"].values()),
            480,
        )
        for dataset in moments["datasets"].values():
            for point in dataset["points"]:
                self.assertAlmostEqual(point["systematic"], 2.0 * point["stat"])

    def test_source_manifests_and_reference_yoda(self) -> None:
        try:
            import yoda
        except (ImportError, OSError):
            self.skipTest("YODA Python bindings are not active")
        expected_objects = {
            "COMPASS_2013_I1236358": 48,
            "COMPASS_2014_I1278730": 112,
        }
        for identifier, count in expected_objects.items():
            root = ROOT / "data/phenomenology" / identifier
            manifest = json.loads((root / "source-manifest.json").read_text())
            source = ROOT / manifest["source"]["path"]
            self.assertEqual(source.stat().st_size, manifest["source"]["bytes"])
            self.assertEqual(
                hashlib.sha256(source.read_bytes()).hexdigest(),
                manifest["source"]["sha256"],
            )
            with tempfile.TemporaryDirectory() as directory:
                output = str(Path(directory) / "reference.yoda.gz")
                with patch.dict(diagnostic_reference.REFERENCE_YODA_PATHS,
                                {identifier: output}):
                    reference.write_reference_yoda(
                        identifier, reference.validate_vendored(identifier))
                with gzip.open(output, "rt", encoding="utf-8") as stream:
                    objects = sum(
                        line.startswith("BEGIN YODA_ESTIMATE1D") for line in stream
                    )
            self.assertEqual(objects, count)

    def test_hermes_exact_bin_definition_has_no_pseudodata(self) -> None:
        registry = campaign.discover_all()
        measurement = registry["HERMES_2013_I1111237"]
        snapshot = campaign._measurement_snapshot(measurement)
        self.assertEqual(snapshot["binning"]["cells_per_target_species"], 900)
        self.assertEqual(len(snapshot["datasets"]), 24)
        self.assertEqual(len(snapshot["observables"]), 24)
        self.assertFalse(any(
            dataset["data_available"] for dataset in snapshot["datasets"].values()
        ))
        self.assertNotIn("points", next(iter(snapshot["datasets"].values())))


class DiagnosticEstimatorTests(unittest.TestCase):
    def test_signed_covariance_is_retained(self) -> None:
        numerator = {"P": series([4.0], [4.0]), "N": series([2.0], [1.0])}
        denominator = {"P": series([10.0], [10.0]), "N": series([6.0], [6.0])}
        positive = {"P": series([0.0], [5.0]), "N": series([0.0], [2.0])}
        negative = {"P": series([0.0], [1.0]), "N": series([0.0], [3.0])}
        result = postprocess.azimuthal_target_combination(
            numerator, denominator, positive, negative, {"P": .5, "N": .5}
        )
        self.assertAlmostEqual(result["values"][0], 3.0 / 8.0)
        self.assertAlmostEqual(
            result["numerator_denominator_covariance"][0], .75
        )
        self.assertGreater(result["errors"][0], 0.0)

    def test_exponential_fit_recovers_inverse_slope(self) -> None:
        edges = [.01 + .02 * index for index in range(36)] + [.7225]
        expected = .27
        values: list[float] = []
        variances: list[float] = []
        for low, high in zip(edges, edges[1:]):
            value = 1000.0 * expected * (math.exp(-low/expected)-math.exp(-high/expected))
            values.append(value)
            variances.append((.01 * value) ** 2)
        spectra = {
            "P": series(values, variances),
            "N": series([.8 * value for value in values],
                        [.8 ** 2 * variance for variance in variances]),
        }
        result = postprocess.pt2_slope_target_combination(
            spectra, {"P": .5, "N": .5}, edges, 1
        )
        self.assertAlmostEqual(result["values"][0], expected, places=9)
        self.assertGreater(result["errors"][0], 0.0)
        self.assertEqual(len(result["retained_pt2_bins"][0]), 36)


class DiagnosticRegistryTests(unittest.TestCase):
    def test_registry_and_central_plans(self) -> None:
        registry = campaign.discover_all()
        for identifier in ALL_IDS:
            self.assertIn(identifier, registry)
            plan = campaign._plan(registry[identifier], arguments())
            self.assertEqual(len(plan["jobs"]), 4)
            self.assertEqual(
                {job["target_component"] for job in plan["jobs"]}, {"P", "N"}
            )
            self.assertEqual(
                {job["contribution"] for job in plan["jobs"]},
                {"POSNLO", "NEGNLO"},
            )

    def test_rivet_sources_preserve_diagnostic_scope(self) -> None:
        compass = (
            ROOT / "analyses/rivet/dis/COMPASS_2014_I1278730.cc"
        ).read_text()
        hermes = (
            ROOT / "analyses/rivet/dis/HERMES_2013_I1111237.cc"
        ).read_text()
        self.assertIn("AzimuthalInputs_", compass)
        self.assertNotIn("SinPhi", compass)
        self.assertIn("900, 0.0, 900.0", hermes)
        self.assertIn("feynmanX", hermes)
        self.assertNotIn("2.0*SIDISAzimuthal::cosineHarmonic", hermes)


if __name__ == "__main__":
    unittest.main()
