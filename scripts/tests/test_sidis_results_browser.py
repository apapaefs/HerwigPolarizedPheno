from __future__ import annotations

from html.parser import HTMLParser
import json
from pathlib import Path
import sys
import tarfile
import tempfile
import unittest
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import build_sidis_results_browser as browser

IDENTIFIER = "COMPASS_2014_I1278730"


class LinkParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []

    def handle_starttag(self, tag, attrs):
        self.links.extend(value for key, value in attrs if key in {"href", "src"})


def assert_links(testcase, output):
    for path in output.rglob("*.html"):
        parser = LinkParser()
        parser.feed(path.read_text())
        for address in parser.links:
            parsed = urlsplit(address)
            if parsed.scheme in {"https", "http"}:
                continue
            testcase.assertFalse(parsed.scheme, address)
            target = (path.parent / unquote(parsed.path)).resolve()
            testcase.assertTrue(browser.contained(target, output), address)
            testcase.assertTrue(target.is_file(), f"Broken link in {path}: {address}")


class BrowserTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.descriptor = browser.read_json(ROOT / "config/phenomenology" / f"{IDENTIFIER}.json")
        cls.signature = browser.measurement_signature(ROOT, cls.descriptor)

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="sidis-browser-test-")
        self.addCleanup(self.temporary.cleanup)
        self.base = Path(self.temporary.name)
        self.campaign_root = self.base / "campaigns"
        self.campaign = self.campaign_root / IDENTIFIER / "fixture"
        self.output = self.base / "results"
        self.selection = {IDENTIFIER: "fixture"}
        self.manifest = {
            "measurement": IDENTIFIER, "tag": "fixture", "status": "postprocessed",
            "configuration": {"measurement_signature": self.signature, "profile": "central"},
            "jobs": [{"status": "success", "events": 100}],
            "plots": {"index": "plots/html/index.html", "created_at": "2026-09-04"},
            "runtime": {"provenance": {"source_control": {"commit": "fixture-commit"}}},
        }
        self.put("manifest.json", json.dumps(self.manifest))
        # These are packaging fixtures, never physics results or preview images.
        self.put("plots/html/index.html", "<html>Original gallery</html>")
        self.put("plots/html/COMPASS/azimuth #? ü.png", "fixture PNG")
        self.put("plots/html/COMPASS/azimuth #? ü.pdf", "fixture PDF")
        self.put("postprocess/prediction.yoda", "fixture analyzed YODA")
        self.put("postprocess/summary.json", '{"dense_covariance": "not exported"}')
        self.put("postprocess/central.csv", "observable,value\nfixture,1\n")

    def put(self, name, text):
        path = self.campaign / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return path

    def build(self, **kwargs):
        return browser.build(ROOT, self.campaign_root, self.output, self.selection, **kwargs)

    def test_compact_archive_roundtrip_and_every_internal_link(self):
        for name in ("logs/huge.log", "work/events.hepmc", "yoda/shard.yoda",
                     "runs/sample.run", "build/plugin.so", "plots/html/plot.py",
                     "plots/html/plot__data.py", "plots/html/plot.dat", "plots/html/large.log"):
            self.put(name, "do not copy" * 100)
        before = browser.digest(self.campaign / "manifest.json")
        report = self.build(archive=True)
        self.assertEqual(report["plots"], 1)
        self.assertEqual(report["issues"], [])
        self.assertEqual(before, browser.digest(self.campaign / "manifest.json"))
        assert_links(self, self.output)
        inventory = browser.read_json(self.output / "bundle.json")
        for record in inventory["files"]:
            path = self.output / record["path"]
            self.assertEqual(record["sha256"], browser.digest(path))
            self.assertEqual(record["bytes"], path.stat().st_size)
        with tarfile.open(report["archive"]) as archive:
            for member in archive.getmembers():
                self.assertTrue(member.isfile())
                self.assertNotIn("..", Path(member.name).parts)
                self.assertEqual(Path(member.name).parts[0], "results")
                self.assertNotIn(Path(member.name).suffix, {".log", ".py", ".dat", ".run", ".so", ".yoda"})
                self.assertNotIn("summary.json", member.name)
                path = self.base / "unpacked" / member.name
                path.parent.mkdir(parents=True, exist_ok=True)
                with archive.extractfile(member) as payload:
                    path.write_bytes(payload.read())
        assert_links(self, self.base / "unpacked/results")
        metadata = browser.read_json(self.output / "analyses" / IDENTIFIER / "metadata.json")
        self.assertEqual(metadata["campaign"]["source_control"]["commit"], "fixture-commit")
        self.assertEqual(metadata["campaign"]["requested_events_in_successful_shards"], 100)

    def test_format_choices_have_no_broken_downloads(self):
        for fmt in ("png", "pdf"):
            self.output = self.base / fmt
            self.build(formats=fmt)
            assert_links(self, self.output)
            other = "pdf" if fmt == "png" else "png"
            self.assertFalse(list(self.output.rglob("*." + other)))

    def test_dry_run_does_not_create_output(self):
        report = self.build(dry_run=True, archive=True)
        self.assertEqual(report["plots"], 1)
        self.assertFalse(self.output.exists())
        self.assertFalse(self.output.with_suffix(".tar.gz").exists())

    def test_no_overwrite_of_bundle_or_archive(self):
        self.output.mkdir()
        with self.assertRaisesRegex(browser.ExportError, "overwrite"):
            self.build()
        self.output = self.base / "other"
        self.output.with_suffix(".tar.gz").write_text("existing archive")
        with self.assertRaisesRegex(browser.ExportError, "overwrite"):
            self.build(archive=True)

    def test_no_output_inside_input_or_enclosing_input(self):
        for output in (self.campaign / "export", self.base):
            with self.assertRaises(browser.ExportError):
                browser.build(ROOT, self.campaign_root, output, self.selection)

    def test_missing_campaign_explicitly_labelled(self):
        self.selection = {IDENTIFIER: "absent"}
        with self.assertRaisesRegex(browser.ExportError, "manifest is missing"):
            self.build()
        self.assertFalse(self.output.exists())
        self.build(allow_incomplete=True)
        self.assertIn("Incomplete", (self.output / "index.html").read_text())
        self.assertNotIn("fixture PNG", (self.output / "index.html").read_text())
        assert_links(self, self.output)

    def test_failed_shard_not_mistaken_for_completion(self):
        self.manifest["jobs"][0]["status"] = "failed"
        self.put("manifest.json", json.dumps(self.manifest))
        with self.assertRaisesRegex(browser.ExportError, "0/1 shards"):
            self.build()
        report = self.build(allow_incomplete=True)
        self.assertTrue(report["issues"])

    def test_identity_and_signature_mismatch_always_refused(self):
        self.manifest["tag"] = "other"
        self.put("manifest.json", json.dumps(self.manifest))
        with self.assertRaisesRegex(browser.ExportError, "identity mismatch"):
            self.build(allow_incomplete=True)
        self.manifest["tag"] = "fixture"
        self.manifest["configuration"]["measurement_signature"] = "wrong"
        self.put("manifest.json", json.dumps(self.manifest))
        with self.assertRaisesRegex(browser.ExportError, "signature mismatch"):
            self.build(allow_incomplete=True)

    def test_partial_plot_pair_is_reported(self):
        self.put("plots/html/COMPASS/unfinished.png", "partial image")
        with self.assertRaisesRegex(browser.ExportError, "Incomplete PNG/PDF"):
            self.build()

    def test_symlink_cannot_smuggle_logs_or_external_files(self):
        log = self.put("logs/secret.log", "not a plot")
        (self.campaign / "plots/html/leak.png").symlink_to(log)
        with self.assertRaisesRegex(browser.ExportError, "escapes source tree|symlink"):
            self.build(allow_incomplete=True)

    def test_path_traversal_in_tag_refused(self):
        self.selection = {IDENTIFIER: "../../other"}
        with self.assertRaisesRegex(browser.ExportError, "Invalid analysis/tag"):
            self.build()

    def test_source_classification_and_diagnostic_caveat(self):
        self.assertEqual(len(browser.CATALOG), 13)
        _, hermes = browser.metadata(ROOT, "HERMES_2013_I1111237")
        self.assertFalse(hermes["data_available"])
        self.assertEqual(hermes["reference_entries"], 0)
        self.assertEqual(hermes["reference_metadata"]["binning"]["cells_per_target_species"], 900)
        _, compass = browser.metadata(ROOT, IDENTIFIER)
        self.assertEqual(compass["authority"], "Paper tables")
        self.assertEqual(compass["reference_entries"], 480)
        self.assertEqual(compass["reference_metadata"]["selection"]["pt_gev"], [0.1, 1.0])
        _, multiplicity = browser.metadata(ROOT, "HERMES_2013_I1208547")
        self.assertEqual(multiplicity["authority"], "Official archive + HEPData cross-checks")

    def test_signature_parity_with_actual_runner_for_all_sidis(self):
        import run_phenomenology_campaign as runner
        registry = runner.discover_all()
        self.assertEqual(set(browser.CATALOG), {k for k, v in registry.items() if "sidis" in v["process_kind"]})
        for identifier in browser.CATALOG:
            with self.subTest(identifier=identifier):
                descriptor = browser.read_json(ROOT / "config/phenomenology" / f"{identifier}.json")
                self.assertEqual(browser.measurement_signature(ROOT, descriptor), runner._signature(registry[identifier]))


if __name__ == "__main__":
    unittest.main()
