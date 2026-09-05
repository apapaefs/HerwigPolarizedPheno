#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import sys
import tempfile
from unittest.mock import patch
import unittest
from pathlib import Path


DISPOL_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(DISPOL_ROOT / "scripts"))
import phenomenology_reference_data as reference  # noqa: E402
import sidis_tranche_reference_data as tranche_reference  # noqa: E402


SNAPSHOTS = {
    "HERMES_2007_I726689_LEGACY": (
        "data/experimental/HERMES_2007_I726689_LEGACY/reference.json",
        "e5e08390deecb21f3708147765eb965334ac4595112081e176692446e704f6d6",
    ),
    "COMPASS_2010_I843494": (
        "data/experimental/COMPASS_2010_I843494/reference.json",
        "7f2f4f7b62bea2fd2e37f43e5eb206464d36d66228b23ca10ebb1054e24de2df",
    ),
    "STAR_2019_I1708793": (
        "data/phenomenology/STAR_2019_I1708793/reference.json",
        "7f9f738a7bcb438016f159e2de96e53cdda6dee5ea895e8262fad6b99e7c8597",
    ),
    "PHENIX_2023_I2033856": (
        "data/phenomenology/PHENIX_2023_I2033856/reference.json",
        "8c6808346ad42d7a1e532874dad9cde05b3ccdb0811897130b46e9325c8306af",
    ),
}


def load(identifier: str) -> dict:
    return json.loads((DISPOL_ROOT / SNAPSHOTS[identifier][0]).read_text(encoding="utf-8"))


class PinnedReferenceTests(unittest.TestCase):
    def test_every_normalized_byte_and_raw_source_is_pinned(self) -> None:
        """The snapshot hashes lock every coordinate, edge, value and error."""
        for identifier, (relative, expected) in SNAPSHOTS.items():
            payload = (DISPOL_ROOT / relative).read_bytes()
            self.assertEqual(hashlib.sha256(payload).hexdigest(), expected)
            self.assertEqual(reference.validate_vendored(identifier), json.loads(payload))
            for source in reference.SOURCES[identifier]:
                if "raw" in source:
                    raw = (DISPOL_ROOT / source["raw"]).read_bytes()
                    self.assertEqual(hashlib.sha256(raw).hexdigest(), source["sha256"])

    def test_hermes_all_45_points_and_low_q2_mask(self) -> None:
        data = load("HERMES_2007_I726689_LEGACY")
        self.assertEqual(len(data["datasets"]), 19)
        self.assertEqual([len(item["points"]) for item in data["datasets"]],
                         [1, 1, 1, 1, 2, 2, 2, 2] + [3] * 11)
        points = [point for dataset in data["datasets"] for point in dataset["points"]]
        self.assertEqual(len(points), 45)
        self.assertEqual(
            [point["value"] for point in points],
            [0.0205,0.0138,0.0294,0.0347,0.0428,0.0931,0.0561,0.0724,
             0.0522,0.0592,0.0291,0.0613,0.0466,0.0549,0.0724,0.0382,
             0.0876,0.1521,0.0518,0.0631,0.1481,0.0414,0.0888,0.1717,
             0.0505,0.0999,0.1896,0.0433,0.0867,0.1875,0.0676,0.1047,
             0.2142,0.0646,0.1268,0.2581,0.0675,0.1719,0.287,0.0962,
             0.1126,0.3151,0.1307,0.3028,0.5274],
        )
        self.assertEqual(
            [point["stat"] for point in points],
            [0.0249,0.0167,0.0153,0.0173,0.0136,0.0294,0.0239,0.0199,
             0.0313,0.0156,0.036,0.0116,0.0179,0.0158,0.0263,0.0185,
             0.0147,0.0218,0.02,0.0143,0.0197,0.0226,0.0141,0.0189,
             0.027,0.0146,0.0189,0.0324,0.0154,0.0195,0.0354,0.0174,
             0.0222,0.0405,0.021,0.0278,0.0514,0.0278,0.0384,0.0739,
             0.0481,0.0643,0.1522,0.0939,0.1757],
        )
        self.assertEqual(
            [point["systematic_combined"] for point in points],
            [0.0026,0.0022,0.0027,0.0031,0.0032,0.0074,0.0031,0.0057,
             0.0031,0.0048,0.0029,0.0051,0.0028,0.0044,0.0069,0.0026,
             0.0053,0.0104,0.0033,0.0043,0.0098,0.0026,0.0055,0.0108,
             0.0029,0.0061,0.0116,0.0033,0.0056,0.0113,0.0039,0.0064,
             0.0127,0.0055,0.008,0.0151,0.006,0.0103,0.0165,0.0085,
             0.0077,0.0182,0.0108,0.0175,0.0298],
        )
        self.assertEqual(sum(point["nonperturbative_extrapolation"] for point in points), 8)
        self.assertTrue(all(point["nonperturbative_extrapolation"] ==
                            (point["q2_mean"] < 1.0) for point in points))
        self.assertEqual(data["scale_floor_gev"], 1.0)

    def test_compass_all_15_points(self) -> None:
        data = load("COMPASS_2010_I843494")
        points = data["points"]
        self.assertEqual(data["bin_edges"],
                         [0.004,0.005,0.006,0.008,0.01,0.02,0.03,0.04,
                          0.06,0.1,0.15,0.2,0.25,0.35,0.5,0.7])
        self.assertEqual([point["x_mean"] for point in points],
                         [0.0046,0.0055,0.007,0.009,0.0147,0.0247,0.0346,
                          0.0487,0.0765,0.122,0.172,0.222,0.29,0.405,0.568])
        self.assertEqual([point["q2_mean"] for point in points],
                         [1.1,1.2,1.37,1.59,2.14,3.24,4.36,6.05,9.42,
                          14.9,20.9,26.7,34.6,47.1,62.1])
        self.assertEqual([point["value"] for point in points],
                         [0.006,0.019,0.035,0.033,0.047,0.076,0.115,0.13,
                          0.172,0.218,0.286,0.446,0.453,0.594,0.855])
        self.assertEqual([point["stat"] for point in points],
                         [0.017,0.013,0.009,0.01,0.006,0.009,0.012,0.012,
                          0.013,0.017,0.024,0.032,0.033,0.049,0.094])
        self.assertEqual([point["systematic_combined"] for point in points],
                         [0.008,0.006,0.005,0.005,0.004,0.006,0.009,0.01,
                          0.012,0.015,0.02,0.03,0.032,0.043,0.068])

    def test_star_12_al_six_all_and_z_points(self) -> None:
        data = load("STAR_2019_I1708793")["channels"]
        expected = {
            "Wplus": (
                [-1.22,-0.71,-0.24,0.25,0.72,1.22],
                [-0.312,-0.251,-0.331,-0.412,-0.534,-0.482],
                [0.145,0.03,0.023,0.023,0.029,0.14],
                [0.012,0.011,0.007,0.005,0.013,0.014],
                [(0.25,0.016,0.042,0.011),(0.72,0.072,0.054,0.011),
                 (1.24,0.0,0.262,0.028)],
            ),
            "Wminus": (
                [-1.25,-0.75,-0.27,0.26,0.75,1.25],
                [0.241,0.26,0.281,0.239,0.385,0.205],
                [0.146,0.051,0.056,0.056,0.051,0.148],
                [0.002,0.004,0.001,0.002,0.004,0.002],
                [(0.27,-0.012,0.101,0.019),(0.74,-0.028,0.092,0.02),
                 (1.27,-0.147,0.26,0.038)],
            ),
        }
        for channel, (eta, values, stat, systematic, all_points) in expected.items():
            points = data[channel]["points"]
            self.assertEqual([point["eta"] for point in points], eta)
            self.assertEqual([point["value"] for point in points], values)
            self.assertEqual([point["stat"] for point in points], stat)
            self.assertEqual([point["systematic_combined"] for point in points], systematic)
            self.assertEqual(
                [(point["abs_eta"], point["value"], point["stat"],
                  point["systematic_combined"]) for point in data[channel]["all"]["points"]],
                all_points,
            )
        self.assertEqual(data["Zgamma"]["points"], [{
            "coordinate": 0.0, "point": 1, "stat": 0.07,
            "systematic_combined": 0.0,
            "systematics": {"published_negligible": 0.0}, "value": -0.04,
        }])

    def test_phenix_18_plus_18_plus_7_and_globals(self) -> None:
        data = load("PHENIX_2023_I2033856")
        channels = data["channels"]
        edges = [6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,12.0,
                 14.0,16.0,18.0,20.0,22.0,24.0,26.0,28.0,30.0]
        for name in ("inclusive_cross_section", "isolated_cross_section"):
            self.assertEqual(channels[name]["bin_edges"], edges)
            self.assertEqual(len(channels[name]["points"]), 18)
        self.assertEqual(
            [point["value"] for point in channels["inclusive_cross_section"]["points"]],
            [2322.0,1301.0,809.7,564.4,364.8,250.8,166.6,106.0,61.08,
             17.76,7.969,2.705,1.311,0.8033,0.4969,0.2487,0.1736,0.08643],
        )
        self.assertEqual(
            [point["value"] for point in channels["isolated_cross_section"]["points"]],
            [1088.0,640.3,404.7,279.5,185.5,133.1,97.66,64.92,35.89,
             13.3,5.856,2.837,1.374,0.791,0.4473,0.2639,0.1517,0.09666],
        )
        isolated_all = channels["isolated_all"]
        self.assertEqual(isolated_all["bin_edges"], [6.0,7.0,8.0,9.0,10.0,12.0,15.0,20.0])
        self.assertEqual([point["value"] for point in isolated_all["points"]],
                         [-0.002265,0.002301,0.004449,0.003846,0.01479,
                          0.01111,-0.02071])
        self.assertEqual(data["global_uncertainties"], {
            "all_polarization_relative": 0.066,
            "all_relative_luminosity_absolute": 0.00039,
            "cross_section_luminosity_relative": 0.1,
        })

    def test_generated_yoda_uses_reference_namespace_and_exact_edges(self) -> None:
        try:
            import yoda
        except (ImportError, OSError):
            self.skipTest("YODA Python bindings are not active")
        for identifier, relative in reference.REFERENCE_YODA_PATHS.items():
            if identifier == "COMPASS_2020_I1788430":
                # This generated product is intentionally not tracked. Exercise
                # its writer in temporary storage, even in a fresh checkout.
                snapshot = reference.validate_vendored(identifier)
                with tempfile.TemporaryDirectory() as directory:
                    output = str(Path(directory) / "reference.yoda.gz")
                    with patch.dict(tranche_reference.REFERENCE_YODA_PATHS,
                                    {identifier: output}):
                        reference.write_reference_yoda(identifier, snapshot)
                    objects = yoda.read(output)
            else:
                objects = yoda.read(str(DISPOL_ROOT / relative))
            self.assertTrue(objects)
            self.assertTrue(all(path.startswith("/REF/") for path in objects))
        phenix = yoda.read(str(DISPOL_ROOT / reference.REFERENCE_YODA_PATHS["PHENIX_2023_I2033856"]))
        self.assertEqual(
            list(phenix["/REF/PHENIX_2023_I2033856/d01-x01-y01"].xEdges()),
            [6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,12.0,14.0,
             16.0,18.0,20.0,22.0,24.0,26.0,28.0,30.0],
        )


if __name__ == "__main__":
    unittest.main()
