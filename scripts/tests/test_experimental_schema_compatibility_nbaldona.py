"""Exact E142 schema migration and fixed-target CLI compatibility."""
from __future__ import annotations

import contextlib
import copy
import io
import json
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"scripts"))
import run_experimental_campaign as campaign


class SchemaCompatibilityTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)
        self.current = campaign.get_measurement("E142_1996_I424108")
        self.historical = copy.deepcopy(self.current)
        self.historical.update(json.loads((Path(__file__).parent/"fixtures/e142_legacy_output_schema.json").read_text()))
        self.source = self.directory/"resolved-measurement.json"
        campaign.atomic_write_json(self.source, self.historical)
        self.recorded = campaign.measurement_signature(self.historical)
        self.manifest = {"configuration": {"measurement_signature": self.recorded}}

    def check(self, current=None):
        return campaign._assert_manifest_signatures_current(
            self.manifest, current or self.current, campaign_dir=self.directory,
        )

    def test_exact_migration_reuses_samples_and_preserves_generation_record(self):
        original_manifest = copy.deepcopy(self.manifest)
        original_source = self.source.read_bytes()
        report = self.check()
        self.assertEqual(report["mode"], "e142_output_schema_migration_v1")
        self.assertEqual(report["generation_measurement_signature"], self.recorded)
        self.assertEqual(report["current_measurement_signature"], campaign.measurement_signature(self.current))
        self.assertEqual(self.manifest, original_manifest)
        self.assertEqual(self.source.read_bytes(), original_source)

    def test_changed_coefficients_paths_and_other_settings_are_rejected(self):
        mutations = (
            lambda m: m["outputs"]["Neutron_W050"]["target_components"]["longitudinal"].update(N=.9),
            lambda m: m["outputs"]["He3_W050"].update(a1="/OTHER/path"),
            lambda m: m["physics"].update(target_model="different physics"),
            lambda m: m["campaign"].update(default_shards=999),
            lambda m: m["outputs"]["Neutron_W050"].update(acceptance="different cuts"),
        )
        for mutate in mutations:
            with self.subTest(mutation=mutate):
                current = copy.deepcopy(self.current)
                mutate(current)
                with self.assertRaisesRegex(campaign.CampaignError, "more than the exact"):
                    self.check(current)

    def test_tampered_history_cannot_be_used_to_relabel_generation(self):
        historical = copy.deepcopy(self.historical)
        historical["title"] = "tampered title"
        campaign.atomic_write_json(self.source, historical)
        current = copy.deepcopy(self.current)
        current["title"] = historical["title"]
        with self.assertRaisesRegex(campaign.CampaignError, "historical signature cannot be reproduced"):
            self.check(current)

    def test_missing_invalid_and_nonhistorical_descriptors_are_rejected(self):
        for payload in ("not JSON", "[]", json.dumps(self.current)):
            self.source.write_text(payload)
            with self.subTest(payload=payload[:20]), self.assertRaises(campaign.CampaignError):
                self.check()
        self.source.unlink()
        with self.assertRaises(campaign.CampaignError):
            self.check()

    def test_changed_source_fingerprint_is_not_accepted(self):
        current_signature = campaign.measurement_signature(self.current)
        with mock.patch.object(campaign, "measurement_signature", side_effect=[current_signature, "changed-source"]):
            with self.assertRaisesRegex(campaign.CampaignError, "historical signature cannot be reproduced"):
                self.check()

    def test_stale_comparison_signature_is_still_rejected(self):
        self.manifest["configuration"].update(comparisons=True, comparison_signature="stale")
        with self.assertRaisesRegex(campaign.CampaignError, "comparison-family"):
            self.check()

    def test_plot_requires_reprocessed_outputs_after_migration(self):
        campaign.atomic_write_json(self.directory/campaign.MANIFEST_NAME, self.manifest)
        with mock.patch.object(campaign, "campaign_directory", return_value=self.directory):
            with self.assertRaisesRegex(campaign.CampaignError, "run postprocess before plotting"):
                campaign.plot_campaign(types.SimpleNamespace(tag="test", dry_run=True), self.current)

    def test_postprocess_records_new_signature_without_changing_generation_identity(self):
        manifest = {**self.manifest, "history": [], "jobs": []}
        campaign.atomic_write_json(self.directory/campaign.MANIFEST_NAME, manifest)
        original_source = self.source.read_bytes()
        def write_yoda(yoda, objects, path):
            path.write_text("test output")
        with mock.patch.object(campaign, "campaign_directory", return_value=self.directory), \
                mock.patch.object(campaign, "_require_complete_matrix", return_value={}), \
                mock.patch.object(campaign, "_load_asymmetry_family_products", return_value=({}, {})), \
                mock.patch.object(campaign, "_import_yoda", return_value=None), \
                mock.patch.object(campaign, "_asymmetry_family_objects", return_value=[]), \
                mock.patch.object(campaign, "_write_yoda_objects", side_effect=write_yoda), \
                mock.patch.object(campaign, "build_campaign_monitor_payload", return_value={}), \
                mock.patch.object(campaign, "write_campaign_monitor_files"), \
                contextlib.redirect_stdout(io.StringIO()):
            campaign.postprocess_campaign(types.SimpleNamespace(tag="test", dry_run=False), self.current)
        updated = campaign.load_json(self.directory/campaign.MANIFEST_NAME)
        summary = campaign.load_json(self.directory/"postprocess/summary.json")
        self.assertEqual(updated["configuration"], manifest["configuration"])
        self.assertEqual(updated["postprocess"]["measurement_signature"], campaign.measurement_signature(self.current))
        for record in (updated["postprocess"], summary):
            self.assertEqual(record["input_compatibility"]["generation_measurement_signature"], self.recorded)
            self.assertEqual(record["input_compatibility"]["mode"], "e142_output_schema_migration_v1")
        self.assertEqual(self.source.read_bytes(), original_source)

    def test_sign_report_distinguishes_compatible_from_identical(self):
        campaign.atomic_write_json(self.directory/campaign.MANIFEST_NAME, {**self.manifest, "jobs": []})
        with mock.patch.object(campaign, "campaign_directory", return_value=self.directory), \
                mock.patch.object(campaign, "_audit_campaign_nlo_signs", return_value={}), \
                contextlib.redirect_stdout(io.StringIO()):
            report = campaign.check_campaign_nlo_signs(types.SimpleNamespace(tag="test"), self.current)
        self.assertFalse(report["generation_signature_current"])
        self.assertTrue(report["postprocessing_compatible"])
        self.assertEqual(report["input_compatibility"]["generation_measurement_signature"], self.recorded)

    def test_primary_cli_forwards_plot_metadata_permission_for_postprocess(self):
        import run_phenomenology_campaign as primary
        for command in ("postprocess", "check-nlo-signs"):
            args = primary.make_parser().parse_args([
                command, "--measurement", "COMPASS_2010_I843494", "--tag", "test", "--allow-plot-metadata-refresh",
            ])
            self.assertTrue(campaign.make_parser().parse_args(primary._legacy_arguments(args)).allow_plot_metadata_refresh)




if __name__ == "__main__":
    unittest.main()
