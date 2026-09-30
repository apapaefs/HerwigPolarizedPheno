#!/usr/bin/env python3
"""Independent target mixtures share generation inputs but retain output semantics."""
from __future__ import annotations

import copy
import contextlib
import io
import json
import math
import sys
import tarfile
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import run_experimental_campaign as campaign


def measurement() -> dict:
    return {
        "schema_version": 2, "id": "TEST", "analysis": {"name": "TEST"},
        "reference": {}, "campaign": {},
        "cards": {
            "components": {"P": {}, "N": {}},
            "helicities": {key: [] for key in ("PP", "PM", "MP", "MM")},
            "orders": {"POSNLO": "PositiveNLO", "NEGNLO": "NegativeNLO"},
        },
        "raw_observables": {
            key: {"ordinary": "ordinary", "weighted": "weighted", "covariance": "covariance"}
            for key in ("P", "D")
        },
        "outputs": {
            key: {
                **{name: f"/TEST/{name}_{key}" for name in ("a1", "apar", "parity_pp_mm", "parity_pm_mp")},
                "target": "proton" if key == "P" else "deuteron",
                "observable": f"A1{key.lower()}",
                "target_model": "free proton" if key == "P" else "free proton/neutron impulse approximation",
            }
            for key in ("P", "D")
        },
        "combination": {
            "unpolarized": {key: .25 for key in ("PP", "PM", "MP", "MM")},
            "longitudinal": {"PP": .25, "PM": -.25, "MP": -.25, "MM": .25},
            "target_components": {
                "unpolarized": {"P": 1., "N": 0.},
                "longitudinal": {"P": 1., "N": 0.},
            },
        },
        "physics": {"target_model": "free proton"}, "diagnostics": ["accepted"],
    }


def mixed_measurement() -> dict:
    result = measurement()
    result["outputs"]["D"].update({
        "target_components": {
            "unpolarized": {"P": .5, "N": .5},
            "longitudinal": {"P": .4625, "N": .4625},
        }, "d_state_factor": .925,
    })
    return result


FAMILY = {
    "components": ["P", "N"], "helicities": ["PP", "PM", "MP", "MM"],
    "orders": ["POSNLO", "NEGNLO"], "postprocess": "helicity_asymmetry",
    "label": "NLO+PS",
}
VALUES = {"P": [14., 8., 8., 10.], "N": [2., 4., 3., 1.]}
VARIANCES = {"P": [4., 1., 2., 3.], "N": [5., 2., 4., 6.]}


def load_series(groups, directory, analysis, name, family, component, helicities, orders, *variation):
    factor = {"ordinary": 1., "weighted": 2. if component == "P" else 3.,
              "covariance": math.sqrt(2. if component == "P" else 3.), "accepted": 1.}[name]
    if tuple(helicities) == ("00",):
        return {"00": campaign.BinSeries([0., 1.], [sum(VALUES[component]) / 4.], [1.])}
    return {
        helicity: campaign.BinSeries([0., 1.], [factor * value], [factor**2 * variance])
        for helicity, value, variance in zip(helicities, VALUES[component], VARIANCES[component])
    }


class FakeBin:
    def __init__(self):
        self.value = None
        self.errors = {}

    def setVal(self, value):
        self.value = value

    def setErr(self, low, high, label):
        self.errors[label] = (low, high)


class FakeEstimate:
    def __init__(self, edges, path):
        self.edges, self.path = edges, path
        self.annotations = {}
        self.bins = [FakeBin() for _ in edges[1:]]

    def setTitle(self, value):
        self.annotations["Title"] = value

    def setAnnotation(self, key, value):
        self.annotations[key] = value

    def bin(self, index):
        return self.bins[index - 1]


FAKE_YODA = types.SimpleNamespace(BinnedEstimate1D=FakeEstimate)


class OutputTargetTests(unittest.TestCase):
    def test_ratio_uses_cross_sections_and_component_covariance(self):
        descriptor = mixed_measurement()
        with mock.patch.object(campaign, "_load_component_series", side_effect=load_series):
            results, diagnostics = campaign._load_asymmetry_family_products(
                {}, Path("."), descriptor, "nominal", FAMILY
            )
        proton, deuteron = results["P"], results["D"]
        self.assertAlmostEqual(proton["sigma_uu"].values[0], 10.)
        self.assertAlmostEqual(proton["sigma_ll"].values[0], 2.)
        self.assertAlmostEqual(proton["a1_values"][0], .4)
        self.assertAlmostEqual(deuteron["sigma_uu"].values[0], 6.25)
        self.assertAlmostEqual(deuteron["sigma_ll"].values[0], .4625)
        self.assertAlmostEqual(deuteron["a1_values"][0], .4625 / 6.25)
        # P/N A1 values are .4 and -1.2: averaging their ratios gives -.37,
        # whereas the cross-section-weighted target estimator is +.074.
        self.assertNotAlmostEqual(deuteron["a1_values"][0], .925 * (.4 - 1.2) / 2.)
        numerator_variance = .4625**2 * (4. * 10. + 9. * 17.) / 16.
        denominator_variance = .5**2 * (10. + 17.) / 16.
        covariance = .4625 * .5 * (2. * (4.-1.-2.+3.) + 3. * (5.-2.-4.+6.)) / 16.
        expected_variance = (numerator_variance / 6.25**2
                             + .4625**2 * denominator_variance / 6.25**4
                             - 2. * .4625 * covariance / 6.25**3)
        self.assertAlmostEqual(deuteron["a1_errors"][0]**2, expected_variance)
        self.assertAlmostEqual(proton["parity_pp_mm_values"][0], 4. / 24.)
        self.assertAlmostEqual(deuteron["parity_pp_mm_values"][0], 5. / 27.)
        self.assertEqual(diagnostics["accepted"].values, proton["sigma_uu"].values)

    def test_direct00_uses_each_output_target_and_global_diagnostics(self):
        descriptor = mixed_measurement()
        family = {**FAMILY, "helicities": ["00"], "postprocess": "direct_unpolarized"}
        with mock.patch.object(campaign, "_load_component_series", side_effect=load_series):
            results, diagnostics = campaign._load_direct_unpolarized_products(
                {}, Path("."), descriptor, "unpolarized", family
            )
        self.assertEqual(results["P"].values, [10.])
        self.assertEqual(results["D"].values, [6.25])
        self.assertEqual(diagnostics["accepted"].values, [10.])
        objects, closure = campaign._direct_unpolarized_objects(
            FAKE_YODA, descriptor, results, diagnostics,
            {key: {"sigma_uu": value} for key, value in results.items()}, {},
        )
        self.assertEqual(closure["P"]["values"], [0.])
        self.assertEqual(closure["D"]["values"], [0.])
        self.assertEqual(objects[0].annotations["Target"], "proton")
        self.assertNotIn("DeuteronDStateFactor", objects[0].annotations)
        self.assertEqual(objects[2].annotations["DeuteronDStateFactor"], "0.925")

    def test_asymmetry_annotations_are_per_output(self):
        descriptor = mixed_measurement()
        with mock.patch.object(campaign, "_load_component_series", side_effect=load_series):
            results, _ = campaign._load_asymmetry_family_products({}, Path("."), descriptor, "nominal", FAMILY)
        objects = campaign._asymmetry_family_objects(FAKE_YODA, descriptor, {}, results, {},
                                                    campaign._family_annotations(descriptor, "nominal", FAMILY))
        by_path = {item.path: item for item in objects}
        proton, deuteron = by_path["/TEST/a1_P"].annotations, by_path["/TEST/a1_D"].annotations
        self.assertEqual(proton["Observable"], "A1p")
        self.assertEqual(proton["TargetComponents"], "P")
        self.assertEqual(proton["TargetModel"], "free proton")
        self.assertNotIn("DeuteronDStateFactor", proton)
        self.assertEqual(deuteron["Observable"], "A1d")
        self.assertEqual(deuteron["TargetComponents"], "P,N")
        self.assertEqual(deuteron["DeuteronDStateFactor"], "0.925")

    def test_diagnostic_annotations_use_the_global_default_target(self):
        descriptor = mixed_measurement()
        annotations = campaign._family_annotations(descriptor, "nominal", FAMILY)
        self.assertEqual(annotations["TargetComponents"], "P")
        direct_family = {**FAMILY, "postprocess": "direct_unpolarized", "helicities": ["00"]}
        direct_annotations = campaign._family_annotations(descriptor, "unpolarized", direct_family)
        self.assertEqual(direct_annotations["TargetComponents"], "P")
        objects = campaign._asymmetry_family_objects(
            FAKE_YODA, descriptor, {}, {},
            {"accepted": campaign.BinSeries([0., 1.], [1.], [1.])}, annotations,
        )
        self.assertEqual(objects[0].annotations["TargetComponents"], "P")
        self.assertNotIn("DeuteronDStateFactor", objects[0].annotations)

    def test_missing_wrong_and_nonfinite_coefficients_are_rejected(self):
        for channel, component, value in (("unpolarized", "N", float("nan")),
                                         ("longitudinal", "P", float("inf")),
                                         ("longitudinal", "N", "invalid")):
            for location in ("global", "output"):
                with self.subTest(channel=channel, value=value, location=location):
                    descriptor = mixed_measurement()
                    target = (descriptor["combination"] if location == "global" else descriptor["outputs"]["D"])
                    target["target_components"][channel][component] = value
                    with self.assertRaises(campaign.CampaignError):
                        campaign._validate_measurement(descriptor, Path("TEST.json"))
        for mutation in (lambda target: target["unpolarized"].pop("N"),
                         lambda target: target["longitudinal"].update({"X": 0.}),
                         lambda target: target.update({"extra": {}})):
            descriptor = mixed_measurement()
            mutation(descriptor["outputs"]["D"]["target_components"])
            with self.assertRaises(campaign.CampaignError):
                campaign._validate_measurement(descriptor, Path("TEST.json"))

    def test_single_target_default_remains_compatible(self):
        descriptor = measurement()
        descriptor["cards"].pop("components")
        descriptor["combination"].pop("target_components")
        family = {**FAMILY, "components": [""]}
        campaign._validate_measurement(descriptor, Path("TEST.json"))
        self.assertEqual(campaign._target_component_coefficients(descriptor, family, "unpolarized", "P"), {"": 1.})

    def test_pdf_variation_postprocess_preserves_output_target_metadata(self):
        descriptor = mixed_measurement()
        descriptor["reference"]["snapshot"] = "reference.json"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            (root / "reference.json").write_text("{}")
            campaign.atomic_write_json(root / campaign.MANIFEST_NAME, {
                "configuration": {"variation_points": [[0, 0, 1.], [1, 0, 1.]]},
            })
            with mock.patch.object(campaign, "DISPOL_ROOT", root), \
                    mock.patch.object(campaign, "campaign_directory", return_value=root), \
                    mock.patch.object(campaign, "_assert_manifest_signatures_current"), \
                    mock.patch.object(campaign, "_require_complete_matrix", return_value={}), \
                    mock.patch.object(campaign, "campaign_family_specs", return_value={"nominal": FAMILY}), \
                    mock.patch.object(campaign, "_load_component_series", side_effect=load_series), \
                    mock.patch.object(campaign, "_import_yoda", return_value=FAKE_YODA), \
                    mock.patch.object(campaign, "_asymmetry_family_objects", wraps=campaign._asymmetry_family_objects) as objects, \
                    contextlib.redirect_stdout(io.StringIO()):
                campaign.postprocess_campaign(types.SimpleNamespace(tag="test", dry_run=True), descriptor)
            variation = objects.call_args_list[1].args[1]
            output = variation["outputs"]["D"]
            self.assertEqual(output["target_components"], descriptor["outputs"]["D"]["target_components"])
            self.assertEqual(output["target"], "deuteron")
            self.assertEqual(output["observable"], "A1d")
            self.assertEqual(output["d_state_factor"], .925)
            self.assertEqual(output["a1"], "/TEST/VARIATIONS/p001-u000-mu1/a1_D")
            self.assertEqual(descriptor["outputs"]["D"]["a1"], "/TEST/a1_D")


def tar_fixture():
    payload = io.BytesIO()
    sources, datasets = [], []
    with tarfile.open(fileobj=payload, mode="w") as archive:
        for key, value in (("P", .1), ("D", -.2)):
            member = f"{key}.dat"
            raw = f"0.2 3.0 {value} 0.02 0.03\n".encode("ascii")
            info = tarfile.TarInfo(member)
            info.size = len(raw)
            archive.addfile(info, io.BytesIO(raw))
            sources.append({"id": key, "archive_member": member,
                            "member_sha256": campaign.sha256_bytes(raw), "expected_rows": 1})
            datasets.append({
                "id": key, "observable": f"A1{key.lower()}", "selection": key,
                "rivet_path": f"/REF/TEST/d14-x01-y0{1 if key == 'P' else 2}",
                "plot_axis": "x_mean", "bin_edges": [.1, .3],
                "provenance": {"hepdata_table": "Table 14"},
                "points": [{"bin": 1, "x_mean": .2, "q2_mean": 3., "value": value,
                            "stat": .02, "systematic_combined": .03, "systematics": {"sys": .03}}],
            })
    return payload.getvalue(), sources, {"schema_version": 2, "measurement": "TEST",
                                         "provenance": {}, "datasets": datasets}


class MultiReferenceTests(unittest.TestCase):
    def test_each_member_and_snapshot_value_is_verified(self):
        payload, sources, snapshot = tar_fixture()
        descriptor = {"id": "TEST", "reference": {"datasets": sources}}
        rows = campaign._multi_tar_reference_rows(descriptor, snapshot, payload)
        self.assertEqual(rows["P"][0]["value"], .1)
        self.assertEqual(rows["D"][0]["value"], -.2)
        for mutate in (lambda item: item["reference"]["datasets"][1].update({"member_sha256": "0" * 64}),
                       lambda item: item["reference"]["datasets"][1].update({"id": "P"}),
                       lambda item: item["reference"]["datasets"][1].update({"id": "unknown"})):
            bad = copy.deepcopy(descriptor)
            mutate(bad)
            with self.assertRaises(campaign.CampaignError):
                campaign._multi_tar_reference_rows(bad, snapshot, payload)
        snapshot["datasets"][1]["points"][0]["value"] += .01
        with self.assertRaises(campaign.CampaignError):
            campaign._multi_tar_reference_rows(descriptor, snapshot, payload)

    def test_fetch_validates_common_archive_and_writes_both_datasets(self):
        payload, sources, snapshot = tar_fixture()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            reference_path = root / "reference.json"
            reference_path.write_text(json.dumps(snapshot))
            destination = root / "reference.yoda"
            descriptor = {
                "id": "TEST", "analysis": {"reference_yoda": str(destination)},
                "reference": {"format": "tar-five-column-multidataset", "datasets": sources,
                              "snapshot": str(reference_path), "source_url": "https://example.test/data",
                              "source_sha256": campaign.sha256_bytes(payload)},
            }
            with mock.patch.object(campaign, "CAMPAIGN_ROOT", root), \
                    mock.patch.object(campaign, "DISPOL_ROOT", root), \
                    mock.patch.object(campaign.urllib.request, "urlopen", side_effect=lambda *args, **kwargs: io.BytesIO(payload)), \
                    mock.patch.object(campaign, "write_reference_yoda") as writer:
                cached = campaign.fetch_reference_data(descriptor)
                self.assertEqual(cached.read_bytes(), payload)
                self.assertEqual(writer.call_args.args, (snapshot, destination))
                validated = campaign.load_json(cached.parent / "validation.json")
                self.assertEqual(set(validated["rows"]), {"P", "D"})
                descriptor["reference"]["source_sha256"] = "0" * 64
                with self.assertRaisesRegex(campaign.CampaignError, "Source checksum mismatch"):
                    campaign.fetch_reference_data(descriptor)

    def test_x_reference_yoda_writer_is_generic_and_preserves_uncertainties(self):
        _, _, snapshot = tar_fixture()
        with mock.patch.object(campaign, "_import_yoda", return_value=FAKE_YODA), \
                mock.patch.object(campaign, "_write_yoda_objects") as writer:
            campaign.write_reference_yoda(snapshot, Path("reference.yoda"))
        objects = writer.call_args.args[1]
        self.assertEqual(len(objects), 2)
        self.assertEqual(objects[1].path, "/REF/TEST/d14-x01-y02")
        self.assertEqual(objects[1].edges, [.1, .3])
        self.assertEqual(objects[1].annotations["PlotAxis"], "x_mean")
        self.assertEqual(objects[1].bins[0].value, -.2)
        self.assertEqual(objects[1].bins[0].errors, {"stat": (-.02, .02), "sys": (-.03, .03)})

    def test_q2_reference_writer_keeps_legacy_delegation(self):
        legacy = {"measurement": "LEGACY", "datasets": [{"q2_edges": [1., 2.]}]}
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary) / "legacy.yoda"
            delegate = mock.Mock(return_value=destination)
            module = types.SimpleNamespace(write_reference_yoda=delegate)
            with mock.patch.dict(sys.modules, {"phenomenology_reference_data": module}):
                campaign.write_reference_yoda(legacy, destination)
            delegate.assert_called_once_with("LEGACY", legacy)

    def test_overlay_uses_declared_x_means_and_retains_legacy_q2_default(self):
        _, _, snapshot = tar_fixture()
        points = campaign._fixed_reference_overlay_points({}, snapshot, "d14-x01-y02")
        self.assertEqual(points[0]["plot_x"], .2)
        del snapshot["datasets"][1]["plot_axis"]
        points = campaign._fixed_reference_overlay_points({}, snapshot, "d14-x01-y02")
        self.assertEqual(points[0]["plot_x"], 3.)


if __name__ == "__main__":
    unittest.main()
