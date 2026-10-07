"""The compact aperture gallery is portable and never reuses damaged output."""
from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import hermes_acceptance_control as gallery


class AcceptanceGalleryTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="hermes-acceptance-gallery-")
        self.addCleanup(temporary.cleanup)
        self.campaign = Path(temporary.name)
        self.source = self.campaign / "postprocess/summary.json"
        self.source.parent.mkdir()
        rows = []
        for target, selection in (("P", "Q2GT1"), ("D", "D_Q2GT1")):
            for selected in (f"Born{target}_X05", selection):
                primary = {"selection": selected, "bin": 1, "bin_low": .0212,
                           "bin_high": .0296, "a_parallel": .1,
                           "a_parallel_stat": .02, "a1": .2, "a1_stat": .03}
                ring = dict(primary, selection="RingControl_" + selected,
                            a_parallel=.12, a1=.21,
                            ring_minus_primary_a_parallel=.02,
                            ring_minus_primary_a_parallel_stat=.01)
                rows.extend((primary, ring))
        self.summary = {"measurement": gallery.MEASUREMENT,
                        "acceptance": "rectangle_intersect_polar_ring", "bins": rows}
        self.source.write_text(json.dumps(self.summary), encoding="utf-8")

    def test_complete_cache_reuses_and_corruption_creates_preserved_sibling(self):
        first = gallery.ensure_gallery(self.campaign, self.campaign / "plots")
        directory = self.campaign / first["directory"]
        original = {path.name: path.read_bytes() for path in directory.iterdir()}
        self.assertTrue(first["created"])
        self.assertTrue(gallery._valid_cache(directory, first["cache_key"]))
        self.assertEqual(set(original), gallery.OUTPUTS | {gallery.MANIFEST})
        directory.relative_to(self.campaign / "plots" / gallery.MEASUREMENT)
        with mock.patch("matplotlib.pyplot.subplots", side_effect=AssertionError("Cache rerendered")):
            repeated = gallery.ensure_gallery(self.campaign, self.campaign / "plots")
        self.assertEqual(repeated["directory"], first["directory"])
        self.assertFalse(repeated["created"])
        damaged = directory / "AcceptanceControl_P.png"
        damaged.write_bytes(b"")
        second = gallery.ensure_gallery(self.campaign, self.campaign / "plots")
        recovered = self.campaign / second["directory"]
        self.assertNotEqual(recovered, directory)
        self.assertEqual(recovered.name, directory.name + "-001")
        self.assertTrue(gallery._valid_cache(recovered, second["cache_key"]))
        self.assertEqual(damaged.read_bytes(), b"")
        for name, expected in original.items():
            if name != damaged.name:
                self.assertEqual((directory / name).read_bytes(), expected)

    def test_wrong_measurement_and_pre_v5_acceptance_fail_before_writing(self):
        for key, value, message in (("measurement", "COMPASS_2017_I1501480", "HERMES_2007"),
                                    ("acceptance", "polar_ring_control", "v5 rectangle/ring")):
            summary = copy.deepcopy(self.summary)
            summary[key] = value
            self.source.write_text(json.dumps(summary), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, message):
                gallery.ensure_gallery(self.campaign, self.campaign / "plots")
            self.assertFalse((self.campaign / "plots").exists())


if __name__ == "__main__":
    unittest.main()
