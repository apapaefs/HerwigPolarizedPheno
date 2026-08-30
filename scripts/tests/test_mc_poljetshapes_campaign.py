#!/usr/bin/env python3
"""Contracts for the combined particle-level shower-spin measurement."""

from __future__ import annotations

import argparse
import copy
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import run_experimental_campaign as experimental  # noqa: E402
import run_mc_poljetshapes_campaign as wrapper  # noqa: E402
import run_phenomenology_campaign as campaign  # noqa: E402


def arguments(**updates: object) -> argparse.Namespace:
    values: dict[str, object] = {
        "profile": "central", "jobs": 100, "shards": 1,
        "seed_base": 8207000, "lo_events": 500000,
        "posnlo_events": None, "negnlo_events": None, "smoke": False,
        "polarized_pdf_members": None, "unpolarized_pdf_members": None,
        "scales": None, "families": "nominal,shower_spin_off",
        "jet_kt_min_gev": None, "nominal_prediction": None,
        "plot_comparisons": True,
    }
    values.update(updates)
    return argparse.Namespace(**values)


def series(values: list[float], variances: list[float],
           edges: list[float] | None = None) -> experimental.BinSeries:
    return experimental.BinSeries(
        edges or [float(index) for index in range(len(values) + 1)],
        values, variances,
    )


def helicities(item: experimental.BinSeries) -> dict[str, experimental.BinSeries]:
    return {label: copy.deepcopy(item) for label in campaign.DENOMINATOR}


class MCPOLJETSHAPESCampaignTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.measurement = campaign.discover_pp_registry()["MC_POLJETSHAPES"]

    def test_particle_level_descriptor_and_snapshot(self) -> None:
        snapshot = campaign._measurement_snapshot(self.measurement)
        self.assertEqual(len(snapshot["observables"]), 75)
        self.assertEqual(snapshot["final_state"]["particles"],
                         "all stable visible particles")
        self.assertFalse(snapshot["final_state"]["truth_information"])
        self.assertEqual(snapshot["jet_definition"]["algorithm"], "anti-kT")
        self.assertEqual(snapshot["declustering"]["algorithm"],
                         "Cambridge/Aachen")
        self.assertEqual(snapshot["rates"]["jet_pt_thresholds_gev"],
                         [2, 3, 4, 5, 6, 8, 10])
        self.assertEqual(
            snapshot["declustering"]["working_points"]["hardshare_kt1"],
            {
                "z1_min": 0.25, "z1_max": 0.4, "z2_min": 0.35,
                "kt1_min_gev": 1.0, "kt2_min_gev": 1.0,
            },
        )
        self.assertFalse(
            snapshot["enrichment_proxies"]["truth_flavour_used"]
        )
        self.assertEqual(
            snapshot["resolved_radiation"]["fraction_edges"],
            [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0],
        )
        configured = self.measurement["channels"]["combined"]["raw_objects"]
        self.assertEqual(set(configured), set(snapshot["observables"]))
        channel = self.measurement["channels"]["combined"]
        self.assertEqual(
            set(channel["conditional_angular_observables"]),
            {"resolved_dphi31_vs_pt31", "resolved_dpsi34_vs_pt41"},
        )
        self.assertEqual(
            self.measurement["statistics_projection"][
                "candidate_events_per_helicity_family"
            ],
            [250000000, 500000000, 1000000000],
        )
        self.assertIn(
            "C2LL_resolved_dpsi34_vs_pt41",
            self.measurement["shard_block_covariance"][
                "prediction_observables"
            ],
        )
        self.assertTrue(all(
            definition["helicity_resolved"]
            for definition in snapshot["observables"].values()
        ))

    def test_rivet_numeric_option_canonicalization_is_exact(self) -> None:
        declared = (
            "MC_POLDIJETS:ETAMAX=1.5:LEVEL=HADRON:PTJ1MIN=5:"
            "PTJ2MIN=4:PTJ3MIN=2:PTJ4MIN=2:R=0.5"
        )
        canonical = (
            "MC_POLDIJETS:ETAMAX=1.5:LEVEL=HADRON:PTJ1MIN=5.0:"
            "PTJ2MIN=4.0:PTJ3MIN=2.0:PTJ4MIN=2.0:R=0.5"
        )
        self.assertTrue(
            experimental._rivet_instance_equivalent(declared, canonical)
        )
        self.assertFalse(
            experimental._rivet_instance_equivalent(
                declared,
                canonical.replace("PTJ2MIN=4.0", "PTJ2MIN=4.1"),
            )
        )

    def test_combined_plugin_instances_and_card_spin_delta(self) -> None:
        options = campaign._resolved_options(arguments(), self.measurement)
        jobs = campaign.build_job_matrix(self.measurement, options)
        self.assertEqual(len(jobs), 8)
        self.assertEqual(len({job["seed"] for job in jobs}), 8)
        for job in jobs:
            self.assertEqual(set(job["analysis_instances"]),
                             {"MC_POLDIJETS", "MC_POLJETSHAPES"})
            self.assertIn("LEVEL=HADRON",
                          job["analysis_instances"]["MC_POLDIJETS"])
            self.assertNotIn("LEVEL=", job["analysis_instances"]["MC_POLJETSHAPES"])
        nominal = next(job for job in jobs
                       if job["family"] == "nominal" and job["helicity"] == "PP")
        control = next(job for job in jobs
                       if job["family"] == "shower_spin_off" and job["helicity"] == "PP")
        nominal_card = campaign._card_text(self.measurement, nominal)
        control_card = campaign._card_text(self.measurement, control)
        common_card = (
            ROOT / self.measurement["cards"]["directory"]
            / self.measurement["cards"]["common"]
        ).read_text()
        for analysis in ("MC_POLDIJETS", "MC_POLJETSHAPES"):
            self.assertEqual(common_card.count(analysis + ":"), 1)
        self.assertIn("ShowerHandler:SpinCorrelations Yes", nominal_card)
        self.assertIn("ShowerHandler:SpinCorrelations No", control_card)
        normalized_nominal = nominal_card.replace(
            "set /Herwig/Shower/ShowerHandler:SpinCorrelations Yes\n", ""
        )
        normalized_control = control_card.replace(
            "set /Herwig/Shower/ShowerHandler:SpinCorrelations No\n", ""
        ).replace(control["stem"], nominal["stem"])
        self.assertEqual(normalized_nominal, normalized_control)

    def test_manifest_records_both_exact_analysis_instances(self) -> None:
        options = campaign._resolved_options(arguments(), self.measurement)
        plan = {"options": options}
        config = campaign._manifest_configuration(
            self.measurement,
            argparse.Namespace(tag="unit"),
            plan,
        )
        self.assertEqual(set(config["analysis_instances"]),
                         {"nominal", "shower_spin_off"})
        for instances in config["analysis_instances"].values():
            self.assertEqual(set(instances), {"MC_POLDIJETS", "MC_POLJETSHAPES"})
        changed = copy.deepcopy(self.measurement)
        changed["analysis"]["companions"][0]["analysis_options"]["R"] = 0.6
        self.assertNotEqual(campaign._signature(self.measurement),
                            campaign._signature(changed))

    def test_geometry_source_contract_is_particle_only_and_deterministic(self) -> None:
        source = (ROOT / self.measurement["analysis"]["source"]).read_text()
        for token in (
            "VisibleFinalState", "fastjet::cambridge_algorithm",
            "declusterHardBranch", "highestKtSplit", "signedPlaneAngle",
            "for (size_t i = 0; i + 2 < particles.size(); ++i)",
            "eventWeights[angleBin(angle)] +=",
            "LorentzTransform::mkFrameTransform(total)",
            "resolved_dphi31_vs_pt31", "resolved_dpsi34_vs_pt41",
            "proxies.constituents >= 8", "proxies.ptD < 0.45",
            "highestKtSplit(primary, 0.25",
        ):
            self.assertIn(token, source)
        for forbidden in ("hardPartonJets", "shower history", ".pid()"):
            self.assertNotIn(forbidden, source)

    def test_helicity_arithmetic_shapes_moments_and_rate_identities(self) -> None:
        angle_edges = [-math.pi, -math.pi/2, 0.0, math.pi/2, math.pi]
        angle = {
            "PP": series([3, 1, 1, 3], [0.04]*4, angle_edges),
            "PM": series([2, 2, 2, 2], [0.04]*4, angle_edges),
            "MP": series([2, 2, 2, 2], [0.04]*4, angle_edges),
            "MM": series([3, 1, 1, 3], [0.04]*4, angle_edges),
        }
        denominator = helicities(series([10, 10], [0.1, 0.1]))
        ge3 = helicities(series([6, 4], [0.06, 0.04]))
        ge4 = helicities(series([3, 1], [0.03, 0.01]))
        prediction = campaign._mc_poljetshapes_prediction(
            {
                "angle": angle,
                "dijet_threshold_denominator": denominator,
                "ge3_threshold": ge3,
                "ge4_threshold": ge4,
            },
            ["angle"],
        )
        for prefix in ("SigmaPP", "SigmaPM", "SigmaMP", "SigmaMM",
                       "SigmaUU", "DeltaSigmaLL", "ALL"):
            self.assertIn(f"{prefix}_angle", prediction)
        shape = prediction["ShapeUU_angle"]
        area = sum(
            value * (shape["edges"][index + 1] - shape["edges"][index])
            for index, value in enumerate(shape["values"])
        )
        self.assertAlmostEqual(area, 1.0)
        self.assertIn("A2UU_angle", prediction)
        self.assertIn("A2LL_angle", prediction)
        self.assertEqual(prediction["R32_UU"]["values"], [0.6, 0.4])
        self.assertEqual(prediction["R43_UU"]["values"], [0.5, 0.25])
        for r32, veto in zip(prediction["R32_UU"]["values"],
                             prediction["ThirdJetVeto_UU"]["values"]):
            self.assertAlmostEqual(r32 + veto, 1.0)

    def test_conditional_resolved_radiation_moments(self) -> None:
        angle_bins = 8
        contents = [
            2.0 if math.cos(2.0 * (-math.pi + (index + 0.5)
                                    * 2.0 * math.pi / angle_bins)) > 0.0
            else 1.0
            for index in range(angle_bins)
        ]
        flat = series(
            contents + [1.0] * angle_bins,
            [0.01] * (2 * angle_bins),
            [float(index) for index in range(2 * angle_bins + 1)],
        )
        prediction = campaign._mc_poljetshapes_prediction(
            {"conditional": helicities(flat)},
            [],
            {
                "resolved_dphi31_vs_pt31": {
                    "raw_observable": "conditional",
                    "fraction_edges": [0.0, 0.5, 1.0],
                    "angle_bins": angle_bins,
                }
            },
            ["conditional"],
        )
        self.assertNotIn("SigmaUU_conditional", prediction)
        cosine = [
            math.cos(2.0 * (-math.pi + (index + 0.5)
                           * 2.0 * math.pi / angle_bins))
            for index in range(angle_bins)
        ]
        expected = 2.0 * sum(
            weight * value for weight, value in zip(cosine, contents)
        ) / sum(contents)
        self.assertAlmostEqual(
            prediction["C2UU_resolved_dphi31_vs_pt31"]["values"][0],
            expected,
        )
        self.assertAlmostEqual(
            prediction["C2UU_resolved_dphi31_vs_pt31"]["values"][1],
            0.0,
            places=12,
        )
        self.assertAlmostEqual(
            prediction["C2LL_resolved_dphi31_vs_pt31"]["values"][0],
            0.0,
        )

    def test_delete_one_shard_covariance_matrix(self) -> None:
        matrix = campaign._jackknife_matrix(
            [[1.0, 2.0], [2.0, 4.0], [3.0, 6.0]]
        )
        self.assertAlmostEqual(matrix["covariance"][0][0], 4.0 / 3.0)
        self.assertAlmostEqual(matrix["covariance"][0][1], 8.0 / 3.0)
        self.assertAlmostEqual(matrix["covariance"][1][1], 16.0 / 3.0)
        self.assertAlmostEqual(matrix["correlation"][0][1], 1.0)

    def test_multi_object_loader_reads_each_yoda_once(self) -> None:
        jobs = [
            {
                "id": "s1", "output_yoda": "one.yoda", "events": 100,
                "analysis_instances": {"TEST": "TEST:OPT=1"},
            },
            {
                "id": "s2", "output_yoda": "two.yoda", "events": 300,
                "analysis_instances": {"TEST": "TEST:OPT=1"},
            },
        ]
        loaded = (
            {
                "first": series([1.0], [0.1]),
                "second": series([2.0], [0.2]),
            },
            {
                "first": series([3.0], [0.3]),
                "second": series([4.0], [0.4]),
            },
        )
        with mock.patch.object(
            experimental,
            "read_histogram_series_many",
            side_effect=loaded,
        ) as reader:
            combined = campaign._load_series_many(
                jobs,
                Path("/campaign"),
                {
                    "first": ("TEST", "A"),
                    "second": ("TEST", "B"),
                },
            )
        self.assertEqual(reader.call_count, 2)
        self.assertAlmostEqual(combined["first"].values[0], 2.5)
        self.assertAlmostEqual(combined["second"].values[0], 3.5)
        requested_paths = reader.call_args_list[0].args[1]
        self.assertEqual(
            requested_paths,
            {"first": "/TEST:OPT=1/A", "second": "/TEST:OPT=1/B"},
        )

    def test_shard_covariance_uses_one_multi_object_read_per_job(self) -> None:
        raw_names = (
            "jet3_pt", "dijet_threshold_denominator",
            "ge3_threshold", "ge4_threshold",
        )
        output_names = ("SigmaUU_jet3_pt", "R32_UU", "R43_UU")
        measurement = copy.deepcopy(self.measurement)
        measurement["channels"] = {
            "combined": {
                "raw_objects": {
                    name: {
                        "analysis": "MC_POLJETSHAPES", "object": name,
                    }
                    for name in raw_names
                },
                "conditional_angular_observables": {},
                "support_observables": [],
            }
        }
        measurement["shard_block_covariance"] = {
            "raw_observables": list(raw_names),
            "prediction_observables": list(output_names),
        }
        groups = {}
        raw_samples = {}
        predictions = {}
        for family_index, family in enumerate(
            ("nominal", "shower_spin_off")
        ):
            key = (family, "combined", 0, 0, 1.0, "off")
            groups[key] = {}
            raw_samples[key] = {name: {} for name in raw_names}
            for helicity_index, helicity in enumerate(campaign.DENOMINATOR):
                jobs = []
                shard_payloads = {name: [] for name in raw_names}
                for shard in range(1, 4):
                    shift = 0.1 * (
                        family_index + helicity_index + shard
                    )
                    payload = {
                        "jet3_pt": series(
                            [2.0 + shift, 1.0 + shift], [0.01, 0.01]
                        ),
                        "dijet_threshold_denominator": series(
                            [10.0 + shift, 9.0 + shift], [0.02, 0.02]
                        ),
                        "ge3_threshold": series(
                            [6.0 + shift, 4.0 + shift], [0.02, 0.02]
                        ),
                        "ge4_threshold": series(
                            [3.0 + shift, 1.0 + shift], [0.01, 0.01]
                        ),
                    }
                    jobs.append(
                        {
                            "shard": shard,
                            "events": 100,
                            "_payload": payload,
                        }
                    )
                    for name, item in payload.items():
                        shard_payloads[name].append(item)
                groups[key][helicity] = jobs
                for name in raw_names:
                    raw_samples[key][name][helicity] = (
                        experimental.combine_shard_series(
                            shard_payloads[name], [100, 100, 100]
                        )
                    )
            predictions[key] = campaign._mc_poljetshapes_prediction(
                raw_samples[key], []
            )

        def load_one(jobs, _campaign_dir, requested):
            self.assertEqual(len(jobs), 1)
            return {
                name: jobs[0]["_payload"][name]
                for name in requested
            }

        with mock.patch.object(
            campaign, "_load_series_many", side_effect=load_one
        ) as reader:
            covariance, rows = campaign._mc_poljetshapes_shard_covariance(
                measurement,
                groups,
                Path("/unused"),
                raw_samples,
                predictions,
            )
        self.assertEqual(reader.call_count, 24)
        self.assertEqual(covariance["families"]["nominal"]["status"],
                         "complete")
        self.assertEqual(
            set(covariance["spin_on_minus_off"]), set(output_names)
        )
        self.assertIsNotNone(
            covariance["spin_on_minus_off"]["R43_UU"]["chi2"]
        )
        self.assertTrue(rows)

    def test_bounded_statistics_projection_and_command(self) -> None:
        base = series([4.0, 4.0, 4.0, 4.0], [0.001]*4,
                      [-math.pi, -math.pi/2, 0, math.pi/2, math.pi])
        objects = {
            name: helicities(copy.deepcopy(base))
            for name in (
                "dpsi12_j1_hardshare_kt1",
                "dpsi12_j2_hardshare_kt1",
                "dpsi12_j1_hardshare_kt2",
                "dpsi12_j2_hardshare_kt2",
                "resolved_dpsi34",
            )
        }
        nominal = campaign._mc_poljetshapes_prediction(
            objects, list(objects)
        )
        control = copy.deepcopy(nominal)
        keys = {
            "nominal": ("nominal", "combined", 0, 0, 1.0, "off"),
            "shower_spin_off": (
                "shower_spin_off", "combined", 0, 0, 1.0, "off"
            ),
        }
        manifest = {
            "configuration": {
                "families": ["nominal", "shower_spin_off"],
                "lo_events": 500000,
            }
        }
        assessment, differences, ranking = campaign._mc_poljetshapes_assessment(
            self.measurement,
            {keys["nominal"]: nominal, keys["shower_spin_off"]: control},
            {keys["nominal"]: objects, keys["shower_spin_off"]: objects},
            manifest,
        )
        self.assertEqual(len(assessment["tiers"]), 3)
        self.assertLessEqual(
            assessment["bounded_recommendation_events_per_helicity_family"],
            1000000000,
        )
        self.assertIn("--seed-base 8407000", assessment["production_command"])
        self.assertIn(
            "mc_poljetshapes_conditional_spin_",
            assessment["production_command"],
        )
        self.assertIn("--jobs 100", assessment["production_command"])
        self.assertIn("A2UU_dpsi12_j1_hardshare_kt1", differences)
        self.assertTrue(ranking)

    def test_plot_ratio_propagates_both_independent_sample_errors(self) -> None:
        ratio = campaign._independent_ratio(
            {
                "edges": [0.0, 1.0, 2.0, 3.0],
                "values": [2.0, 1.0, 3.0],
                "errors": [0.2, 0.1, 0.3],
            },
            {
                "edges": [0.0, 1.0, 2.0, 3.0],
                "values": [1.0, 0.0, 2.0],
                "errors": [0.1, 0.2, 0.4],
            },
        )
        self.assertAlmostEqual(ratio["values"][0], 2.0)
        self.assertAlmostEqual(ratio["errors"][0], math.sqrt(0.08))
        self.assertIsNone(ratio["values"][1])
        self.assertIsNone(ratio["errors"][1])
        self.assertAlmostEqual(ratio["values"][2], 1.5)
        self.assertAlmostEqual(ratio["errors"][2], math.sqrt(0.1125))

    def test_plot_ratio_selection_excludes_signed_and_moment_objects(self) -> None:
        nominal_key = campaign._variation_id(
            ("nominal", "combined", 0, 0, 1.0, "off")
        )
        control_key = campaign._variation_id(
            ("shower_spin_off", "combined", 0, 0, 1.0, "off")
        )
        prediction = {
            "edges": [0.0, 1.0], "values": [2.0], "errors": [0.2]
        }
        control = {
            "edges": [0.0, 1.0], "values": [1.0], "errors": [0.1]
        }
        ratios = campaign._mc_poljetshapes_plot_ratios(
            self.measurement,
            {
                "variations": {
                    nominal_key: {
                        "SigmaUU_dpsi12_j1_loose": prediction,
                        "ShapeUU_dpsi12_j1_loose": prediction,
                        "DeltaSigmaLL_dpsi12_j1_loose": prediction,
                        "ALL_dpsi12_j1_loose": prediction,
                        "A2UU_dpsi12_j1_loose": prediction,
                    },
                    control_key: {
                        "SigmaUU_dpsi12_j1_loose": control,
                        "ShapeUU_dpsi12_j1_loose": control,
                        "DeltaSigmaLL_dpsi12_j1_loose": control,
                        "ALL_dpsi12_j1_loose": control,
                        "A2UU_dpsi12_j1_loose": control,
                    },
                }
            },
        )
        self.assertEqual(
            set(ratios),
            {
                "SigmaUU_dpsi12_j1_loose",
                "ShapeUU_dpsi12_j1_loose",
            },
        )
        self.assertEqual(ratios["SigmaUU_dpsi12_j1_loose"]["values"], [2.0])

    def test_focused_gallery_links_overlays_differences_and_ratios(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            campaign_dir = Path(temporary) / "campaign"
            output = campaign_dir / "plots" / "html"
            output.mkdir(parents=True)
            (output / "index.html").write_text(
                "<html><body><h1>Complete gallery</h1></body></html>",
                encoding="utf-8",
            )
            assets = (
                Path(
                    "MC_POLJETSHAPES/"
                    "C2UU_resolved_dphi31_vs_pt31.png"
                ),
                Path(
                    "MC_POLJETSHAPES/COMPARISON/"
                    "OnMinusOff_C2UU_resolved_dphi31_vs_pt31.png"
                ),
                Path(
                    "MC_POLJETSHAPES/"
                    "A2UU_dpsi12_j1_hardshare_kt1.png"
                ),
                Path(
                    "MC_POLJETSHAPES/COMPARISON/"
                    "OnMinusOff_A2UU_dpsi12_j1_hardshare_kt1.png"
                ),
                Path(
                    "MC_POLJETSHAPES/"
                    "ShapeUU_eeec_squeezed_j1_hardshare_kt1.png"
                ),
                Path(
                    "ratios/MC_POLJETSHAPES/"
                    "ShapeUU_eeec_squeezed_j1_hardshare_kt1.png"
                ),
                Path(
                    "MC_POLJETSHAPES/"
                    "A2UU_dpsi12_j1_hardshare_kt1_nconst_high.png"
                ),
                Path(
                    "MC_POLJETSHAPES/COMPARISON/"
                    "OnMinusOff_A2UU_dpsi12_j1_hardshare_kt1_nconst_high.png"
                ),
                Path("MC_POLJETSHAPES/A2UU_dpsi12_j1_loose.png"),
                Path(
                    "MC_POLJETSHAPES/COMPARISON/"
                    "OnMinusOff_A2UU_dpsi12_j1_loose.png"
                ),
                Path("MC_POLJETSHAPES/ShapeUU_dpsi12_j1_loose.png"),
                Path("ratios/MC_POLJETSHAPES/ShapeUU_dpsi12_j1_loose.png"),
                Path("MC_POLJETSHAPES/SigmaUU_jet3_pt.png"),
                Path("ratios/MC_POLJETSHAPES/SigmaUU_jet3_pt.png"),
                Path("MC_POLJETSHAPES/ALL_jet4_pt.png"),
                Path(
                    "MC_POLJETSHAPES/COMPARISON/"
                    "OnMinusOff_ALL_jet4_pt.png"
                ),
            )
            for relative in assets:
                path = output / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(b"nonempty")
            postprocess = campaign_dir / "postprocess"
            postprocess.mkdir()
            (postprocess / "shard-block-covariance.json").write_text(
                json.dumps(
                    {
                        "families": {
                            "nominal": {"status": "complete", "blocks": 50},
                            "shower_spin_off": {
                                "status": "complete", "blocks": 50,
                            },
                        },
                        "spin_on_minus_off": {
                            "R43_UU": {"chi2": 2.0, "rank": 4}
                        },
                    }
                ),
                encoding="utf-8",
            )
            focus = campaign._write_mc_poljetshapes_focus_index(
                output, self.measurement
            )
            rendered = focus.read_text(encoding="utf-8")
            root = (output / "index.html").read_text(encoding="utf-8")
            self.assertIn("Inclusive splitting-plane moments", rendered)
            self.assertIn("Inclusive splitting-plane shapes", rendered)
            self.assertIn("Resolved-radiation plane moments", rendered)
            self.assertIn("Hard energy-sharing declusterings", rendered)
            self.assertIn("Hard-sharing squeezed EEEC", rendered)
            self.assertIn("Particle-level gluon-enriched proxies", rendered)
            self.assertIn("Resolved third- and fourth-jet spectra", rendered)
            self.assertIn("Third-jet pT: UU cross section", rendered)
            self.assertIn("Fourth-jet pT: A_LL", rendered)
            self.assertIn("Spin on / spin off", rendered)
            self.assertIn("propagate both independent Monte Carlo errors", rendered)
            self.assertIn("Shard-block covariance audit", rendered)
            self.assertIn("R43_UU", rendered)
            self.assertIn(
                "../ratios/MC_POLJETSHAPES/ShapeUU_dpsi12_j1_loose.png",
                rendered,
            )
            self.assertIn(
                "../ratios/MC_POLJETSHAPES/SigmaUU_jet3_pt.png",
                rendered,
            )
            self.assertIn(
                "../MC_POLJETSHAPES/COMPARISON/OnMinusOff_ALL_jet4_pt.png",
                rendered,
            )
            self.assertIn("focus/index.html", root)

    def test_focused_gallery_promotes_resolved_jet_spectra(self) -> None:
        focused = {
            stem
            for section in campaign.MC_POLJETSHAPES_FOCUS_SECTIONS
            for stem, _, _ in section["plots"]
        }
        self.assertTrue(
            {
                "SigmaUU_jet3_pt",
                "DeltaSigmaLL_jet3_pt",
                "ALL_jet3_pt",
                "SigmaUU_jet4_pt",
                "DeltaSigmaLL_jet4_pt",
                "ALL_jet4_pt",
            }.issubset(focused)
        )
        self.assertTrue(
            {
                "C2UU_resolved_dphi31_vs_pt31",
                "C2LL_resolved_dpsi34_vs_pt41",
                "A2UU_dpsi12_j1_hardshare_kt1",
                "ShapeUU_eeec_squeezed_j1_hardshare_kt1",
                "A2UU_dpsi12_j1_hardshare_kt1_nconst_high",
                "A2UU_dpsi12_j1_hardshare_kt1_ptd_low",
            }.issubset(focused)
        )

    def test_ratio_plot_script_gets_linear_scale_and_unity_line(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            script = Path(temporary) / "ratio.py"
            script.write_text(
                "ax_yScale = 'log'\n"
                "\nlegend_handles = dict() # keep track of handles for the legend\n",
                encoding="utf-8",
            )
            campaign._configure_mc_poljetshapes_ratio_script(script)
            configured = script.read_text(encoding="utf-8")
            self.assertIn("ax_yScale = 'linear'", configured)
            self.assertIn("ax.axhline(1.0", configured)
            self.assertIn("spin-on/spin-off ratio presentation", configured)

    def test_dedicated_runner_pins_measurement(self) -> None:
        self.assertEqual(
            wrapper.pinned_arguments(["full", "--tag", "example"]),
            ["full", "--measurement", "MC_POLJETSHAPES", "--tag", "example"],
        )
        with self.assertRaises(campaign.CampaignError):
            wrapper.pinned_arguments(
                ["full", "--measurement", "MC_POLDIJETS", "--tag", "bad"]
            )

    def test_plot_contract_covers_overlays_moments_and_differences(self) -> None:
        plot = (ROOT / self.measurement["analysis"]["plot"]).read_text()
        self.assertIn("Shape(UU|PP|PM|MP|MM)", plot)
        self.assertIn("(A2UU|A2LL)", plot)
        self.assertIn("C2(UU|LL|PP|PM|MP|MM)", plot)
        self.assertIn("COMPARISON/OnMinusOff", plot)
        self.assertEqual(
            self.measurement["families"]["nominal"]["plot_options"]["LineColor"],
            "#CC3311",
        )
        self.assertEqual(
            self.measurement["families"]["shower_spin_off"]["plot_options"]["LineColor"],
            "#0077BB",
        )
        self.assertIn("ShapeUU_", campaign.MC_POLJETSHAPES_RATIO_PREFIXES)
        workflow = (ROOT / "docs" / "mc-poljetshapes-workflow.md").read_text()
        normalized_workflow = " ".join(workflow.split())
        self.assertIn("plots/html/focus/index.html", workflow)
        self.assertIn("third- and fourth-jet", normalized_workflow)
        self.assertIn("run_mc_poljetshapes_campaign.py plot", workflow)

    def test_conditional_moments_are_primary_plot_objects(self) -> None:
        snapshot = campaign._measurement_snapshot(self.measurement)
        self.assertEqual(
            campaign._reference_path(
                self.measurement,
                "C2UU_resolved_dphi31_vs_pt31",
                snapshot,
            ),
            "/MC_POLJETSHAPES/C2UU_resolved_dphi31_vs_pt31",
        )


if __name__ == "__main__":
    unittest.main()
