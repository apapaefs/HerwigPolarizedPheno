from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import build_results_browser as browser

ID = "COMPASS_2014_I1278730"
NO_PROCESSES = {"available": True, "controllers": [], "worker_directories": []}


class ResultsBrowserTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.registry = browser.registry(ROOT)
        cls.kind, cls.descriptor = cls.registry[ID]
        cls.signature = browser.measurement_signature(ROOT, cls.descriptor, cls.kind)

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="latest-results-test-")
        self.addCleanup(self.temporary.cleanup)
        self.base = Path(self.temporary.name)
        self.root = self.base / "campaigns"
        self.output = self.base / "browser"

    def campaign(self, tag, created="2026-09-01T10:00:00+00:00", *, success=True,
                 plotted=True, smoke=False, signature=None, identifier=ID, index="plots/html/index.html"):
        path = self.root / identifier / tag
        path.mkdir(parents=True)
        descriptor = self.registry[identifier][1]
        kind = self.registry[identifier][0]
        d = {"measurement": identifier, "tag": tag, "created_at": created,
             "updated_at": created, "status": "complete" if success else "running",
             "configuration": {"smoke": smoke, "measurement_signature": signature or browser.measurement_signature(ROOT, descriptor, kind)},
             "jobs": [{"status": "success" if success else "queued", "events": 1000}],
             "runtime": {"provenance": {"source_control": {"commit": "recorded-source"}}}}
        if plotted:
            d["plots"] = {"index": index}
            d["postprocess"] = {"prediction": "postprocess/prediction.yoda", "summary": "postprocess/summary.json",
                                "predictions": [{"path": "postprocess/prediction.yoda"}]}
            for name, value in [(index, "<html>Existing plots</html>"),
                                (str(Path(index).parent / identifier / 'cos1_#? ü.png'), "PNG fixture"),
                                (str(Path(index).parent / identifier / 'cos1_#? ü.pdf'), "PDF fixture"),
                                ("postprocess/prediction.yoda", "Prediction fixture"),
                                ("postprocess/summary.json", json.dumps({"masked_bins": {"first": [3, 4], "second": []}, "goodness_of_fit": {"chi2": 7.5, "points": 8}})),
                                ("postprocess/central.csv", "observable,value\nfixture,1\n")]:
                target = path / name; target.parent.mkdir(parents=True, exist_ok=True); target.write_text(value)
        (path / "manifest.json").write_text(json.dumps(d))
        return path, d

    def inspect(self, **kwargs):
        return browser.inspect_analysis(ROOT, [self.root], ID, self.kind, self.descriptor,
                                        NO_PROCESSES, "both", kwargs.get("include_smoke", False), kwargs.get("tag"))[0]

    def build(self, **kwargs):
        return browser.build(ROOT, [self.root], self.output, identifiers=[ID], observation=NO_PROCESSES, **kwargs)

    def test_all_registry_metadata_and_signature_parity(self):
        import run_phenomenology_campaign as pheno
        import run_experimental_campaign as experimental
        self.assertEqual(set(self.registry), set(pheno.discover_all()))
        native = pheno.discover_pp_registry()
        fixed = experimental.discover_registry()
        for identifier, (kind, descriptor) in self.registry.items():
            with self.subTest(identifier=identifier):
                expected = (experimental.measurement_signature(fixed[identifier]) if kind == "experimental"
                            else browser.open_charm.signature(ROOT, descriptor) if kind == "reference"
                            else pheno._signature(native[identifier]))
                self.assertEqual(browser.measurement_signature(ROOT, descriptor, kind), expected)
                item = browser.metadata(ROOT, identifier, kind, descriptor)
                self.assertEqual(item["id"], identifier)
                if kind == "experimental":
                    self.assertEqual(browser.comparison_signature(ROOT, descriptor), experimental.comparison_signature(fixed[identifier]))
        self.assertEqual(browser.metadata(ROOT, "MC_POLDIJETS", *self.registry["MC_POLDIJETS"])["reference_entries"], 0)

    def test_additional_reference_and_raw_source_corruption_rejected(self):
        kind, descriptor = self.registry["HERMES_2007_I726689"]
        changed = copy.deepcopy(descriptor)
        changed["reference"]["additional_snapshots"][0]["sha256"] = "0" * 64
        with self.assertRaises(browser.common.ExportError):
            browser.measurement_signature(ROOT, changed, kind)

        snapshot = descriptor["reference"]["additional_snapshots"][0]["path"]
        source = json.loads((ROOT / snapshot).read_text(encoding="utf-8"))["source_files"][0]["path"]
        damaged = self.base / "damaged-reference"
        damaged.write_bytes((ROOT / source).read_bytes() + b"\n")
        original = browser.source_file
        with patch.object(browser, "source_file", side_effect=lambda repo, relative:
                          damaged if relative == source else original(repo, relative)):
            with self.assertRaises(browser.common.ExportError):
                browser.measurement_signature(ROOT, descriptor, kind)

    def test_creation_time_not_old_replot_or_lexical_tag_controls_latest(self):
        old, manifest = self.campaign("zzz-older")
        manifest["updated_at"] = "2026-12-01T00:00:00Z"
        (old / "manifest.json").write_text(json.dumps(manifest))
        self.campaign("aaa-newer", "2026-09-02T10:00:00Z")
        self.assertEqual(self.inspect()["selected"]["tag"], "aaa-newer")

    def test_channel_reference_counts_and_authorities(self):
        for identifier, entries, groups in [("PHENIX_2023_I2033856", 43, 3), ("STAR_2019_I1708793", 19, 5)]:
            item = browser.metadata(ROOT, identifier, *self.registry[identifier])
            self.assertEqual(item["reference_entries"], entries)
            self.assertEqual(item["dataset_count"], groups)
            self.assertIn("HEPData", item["authority"])
            self.assertTrue(any("hepdata" in url for _, url in item["links"]))

    def test_smoke_excluded_unless_explicitly_enabled(self):
        self.campaign("production")
        self.campaign("smoke", "2026-09-03T00:00:00Z", smoke=True)
        item = self.inspect()
        self.assertEqual(item["selected"]["tag"], "production")
        self.assertEqual(item["rejected"][0]["reason"], "Smoke test excluded")
        self.assertEqual(self.inspect(include_smoke=True)["selected"]["tag"], "smoke")

    def test_stopped_new_run_does_not_hide_older_ready_result(self):
        self.campaign("finished")
        self.campaign("new-incomplete", "2026-09-03T00:00:00Z", success=False, plotted=False)
        item = self.inspect()
        self.assertEqual(item["selected"]["tag"], "finished")
        self.assertEqual(item["latest"]["tag"], "new-incomplete")
        self.assertIn("Stopped", item["latest"]["activity"])
        self.assertEqual(item["availability"], "partial")
        self.assertIn("newer campaign", item["issues"][0])

    def test_incompatible_history_never_exports_old_physics_plots(self):
        self.campaign("old-physics", signature="different-definition")
        self.campaign("new-incomplete", "2026-09-03T00:00:00Z", success=False, plotted=False)
        item = self.inspect()
        self.assertIsNone(item["selected"])
        self.assertEqual(item["plots"], [])
        self.assertFalse(item["history"][1]["compatible"])
        self.build()
        self.assertFalse(list(self.output.rglob("*.png")))

    def test_explicit_pin_and_invalid_pin(self):
        self.campaign("old")
        self.campaign("new", "2026-09-04T00:00:00Z")
        item = self.inspect(tag="old")
        self.assertEqual(item["selected"]["tag"], "old")
        self.assertIn("Explicitly selected", item["issues"][0])
        with self.assertRaises(browser.common.ExportError):
            self.inspect(tag="absent")
        with self.assertRaises(browser.common.ExportError):
            self.build(selections={ID: "../../escape"})

    def test_partial_plots_are_labelled_and_masks_are_read(self):
        self.campaign("partial", success=False)
        item = self.inspect()
        self.assertEqual(item["availability"], "partial")
        self.assertEqual(len(item["plots"]), 1)
        self.assertEqual(item["diagnostics"]["total_masked_bins"], 2)

    def test_legacy_experimental_products_follow_recorded_paths(self):
        identifier = "HERMES_2007_I726689"
        path, manifest = self.campaign("legacy", identifier=identifier, index="plots/index.html")
        (path / "postprocess/prediction.yoda").rename(path / "postprocess/HERMES-analyzed.yoda")
        (path / "postprocess/central.csv").rename(path / "postprocess/summary.csv")
        manifest["postprocess"] = {"yoda": "postprocess/HERMES-analyzed.yoda", "summary_csv": "postprocess/summary.csv",
                                   "summary_json": "postprocess/summary.json",
                                   "predictions": [{"family": "nominal", "yoda": "postprocess/HERMES-analyzed.yoda"}]}
        (path / "manifest.json").write_text(json.dumps(manifest))
        browser.build(ROOT, [self.root], self.output, identifiers=[identifier], observation=NO_PROCESSES)
        metadata = json.loads((self.output / "analyses" / identifier / "metadata.json").read_text())
        self.assertEqual(metadata["availability"], "ready")
        self.assertTrue((self.output / metadata["downloads"]["Numerical results (CSV)"]).is_file())

    def test_missing_recorded_prediction_family_is_partial(self):
        path, manifest = self.campaign("missing-family")
        manifest["postprocess"] = {"predictions": [{"yoda": "postprocess/absent-comparison.yoda"}]}
        (path / "manifest.json").write_text(json.dumps(manifest))
        item = self.inspect()
        self.assertEqual(item["availability"], "partial")
        self.assertTrue(any("absent-comparison" in issue for issue in item["issues"]))

    def test_missing_creation_time_uses_explicitly_labelled_fallback(self):
        self.campaign("no-date", created=None)
        self.assertIn("mtime", self.inspect()["selected"]["ordering_basis"])

    def test_malformed_manifest_and_source_symlink_are_reported(self):
        good, _ = self.campaign("good")
        bad = self.root / ID / "broken"; bad.mkdir(); (bad / "manifest.json").write_text('{')
        linked = self.root / ID / "linked"; linked.symlink_to(good, target_is_directory=True)
        item = self.inspect()
        self.assertEqual(item["selected"]["tag"], "good")
        self.assertEqual(len(item["rejected"]), 2)

    def test_escaped_plot_index_is_never_copied(self):
        path, manifest = self.campaign("bad-index")
        manifest["plots"]["index"] = "../../secret/index.html"
        (path / "manifest.json").write_text(json.dumps(manifest))
        item = self.inspect()
        self.assertIsNone(item["selected"])
        self.assertIn("escapes", item["rejected"][0]["reason"])

    def test_presentation_only_refresh_reconstructs_original_signature(self):
        plot = self.descriptor["analysis"]["plot"]
        historical = (ROOT / plot).read_bytes() + b'\n# old presentation\n'
        old_signature = browser.measurement_signature(ROOT, self.descriptor, self.kind, historical)
        refresh = {"mode": "presentation_only_plot_metadata_refresh", "current_measurement_signature": self.signature,
                   "generation_measurement_signature": old_signature, "plot_path": plot,
                   "current_plot_sha256": browser.common.digest(ROOT / plot), "generation_source_commit": 'a' * 40,
                   "generation_plot_sha256": browser.hashlib.sha256(historical).hexdigest()}
        manifest = {"configuration": {"measurement_signature": old_signature}, "plots": {"presentation_only_refresh": refresh}}
        with patch.object(browser.subprocess, 'check_output', return_value=historical):
            self.assertTrue(browser.compatible(ROOT, self.descriptor, self.kind, manifest, self.signature)[0])
        with patch.object(browser.subprocess, 'check_output', return_value=b'wrong historical bytes'):
            self.assertFalse(browser.compatible(ROOT, self.descriptor, self.kind, manifest, self.signature)[0])

    def test_actual_controller_observation_and_portable_unknown_state(self):
        self.campaign("unfinished", success=False, plotted=False)
        candidate = self.inspect()["latest"]
        observed = {**NO_PROCESSES, "controllers": [[ID, "unfinished"]]}
        self.assertEqual(browser.activity(candidate, observed), "Controller active")
        self.assertIn("not checked", browser.activity(candidate, {**NO_PROCESSES, "available": False}))

    def test_archive_contains_only_viewing_payload_and_all_plot_paths_resolve(self):
        path, _ = self.campaign("complete")
        (path / 'logs').mkdir(); (path / 'logs/secret.log').write_text('not for export')
        before = browser.common.digest(path / 'manifest.json')
        report = self.build(archive=True)
        self.assertEqual(report['plots'], 1)
        self.assertEqual(before, browser.common.digest(path / 'manifest.json'))
        page = (self.output / 'index.html').read_text()
        payload = page.split('<script id="results-data" type="application/json">')[1].split('</script>')[0]
        data = json.loads(payload)
        for item in data['items']:
            for plot in item['plots']:
                for relative in plot['files'].values(): self.assertTrue((self.output / relative).is_file())
            for relative in item['downloads'].values(): self.assertTrue((self.output / relative).is_file())
        with tarfile.open(report['archive']) as archive:
            self.assertTrue(all(m.isfile() for m in archive.getmembers()))
            self.assertFalse(any(Path(m.name).suffix in {'.log', '.yoda', '.run', '.so', '.py'} for m in archive.getmembers()))
        inventory = json.loads((self.output / 'bundle.json').read_text())
        for record in inventory['files']:
            self.assertEqual(record['sha256'], browser.common.digest(self.output / record['path']))

    def test_no_output_overwrites_and_dry_run_is_read_only(self):
        self.campaign('complete')
        self.build(dry_run=True)
        self.assertFalse(self.output.exists())
        self.build()
        with self.assertRaises(browser.common.ExportError): self.build()
        with self.assertRaises(browser.common.ExportError):
            browser.build(ROOT, [self.root], self.root / 'output', identifiers=[ID])

    def test_json_embedding_is_safe_and_nonfinite_diagnostics_are_explicit(self):
        value = browser.compact({'chi2': float('nan')})
        self.assertIn('non-finite', value['chi2'])
        page = browser.index_html({'title': '</script><script>bad()</script>'})
        self.assertNotIn('<script>bad()', page)
        self.assertIn('\\u003c/script\\u003e', page)


if __name__ == '__main__':
    unittest.main()
