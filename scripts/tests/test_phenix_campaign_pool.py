"""Normalization and input-integrity checks for independent PHENIX pooling."""
import copy
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import combine_PHENIX_2023_I2033856 as pool


def group(events, value, variance):
    return {"events": events, "objects": {
        name: pool.exp.BinSeries([6, 7], [value], [variance]) for name in pool.OBJECTS}}


def manifest():
    jobs = []
    for i, (channel, helicity) in enumerate(pool.itertools.product(
            ("photon", "fragmentation"), pool.HELICITIES)):
        jobs.append({"id": str(i), "output_yoda": f"yoda/{i}.yoda", "seed": 100+i,
                     "channel": channel, "helicity": helicity, "hard_pt_range_gev": [4, 8],
                     "family": "nominal", "order": "LO", "mpi": "on", "scale": 1.0,
                     "polarized_pdf_member": 0, "unpolarized_pdf_member": 0,
                     "analysis_instance": pool.MEASUREMENT, "status": "success", "events": 100,
                     "shard": 1})
    return {"status": "complete", "jobs": jobs, "configuration": {
        "shards": 1, "hard_pt_ranges_gev": [[4, 8]], "lo_events": 100}}


class PoolingTests(unittest.TestCase):
    def test_plot_only_reuses_saved_yoda_without_loading_campaigns(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            with mock.patch.object(pool, "plot_results") as plot, \
                    mock.patch.object(pool, "collect") as collect:
                self.assertEqual(pool.main(["--plot-only", "--output-dir", temp]), 0)
                plot.assert_called_once_with(directory.resolve(), None)
                collect.assert_not_called()

    def test_native_yoda_roundtrip_and_standard_paths(self):
        try:
            yoda = pool.exp._import_yoda()
        except pool.exp.CampaignError:
            self.skipTest("Native YODA is not available in this Python environment")
        paths = ("d01-x01-y01", "d01-x01-y02", "d02-x01-y01")
        reference = {"channels": {name: {"rivet_path": f"/{pool.MEASUREMENT}/{stem}"}
                                 for name, stem in zip(pool.OBJECTS, paths)}}
        prediction = {name: {"edges": [6, 7, 8], "values": [0.234567891, -0.01],
                             "errors": [0.0034567891, 0.004]} for name in pool.OBJECTS}
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "prediction.yoda"
            pool.write_yoda(path, prediction, reference)
            objects = yoda.read(str(path))
            self.assertEqual(set(objects), {f"/{pool.MEASUREMENT}/{p}" for p in paths})
            for obj in objects.values():
                self.assertEqual(list(obj.xEdges()), [6, 7, 8])
                for b, value, error in zip(obj.bins(), [0.234567891, -0.01], [0.0034567891, 0.004]):
                    self.assertTrue(math.isclose(b.val(), value, rel_tol=5.1e-7))
                    self.assertTrue(math.isclose(b.errAvg("stat"), error, rel_tol=5.1e-7))

    def test_unequal_counts_use_event_fractions_and_squared_variances(self):
        key = ("photon", "PP", 4, 8)
        combined = pool.pool_groups([{key: group(100, 2, 4)}, {key: group(300, 6, 16)}])
        self.assertEqual(combined[key]["events"], 400)
        series = combined[key]["objects"]["isolated_all"]
        self.assertEqual(series.values, [5])
        self.assertEqual(series.variances, [9.25])

    def test_three_runs_and_associativity(self):
        key = ("photon", "PP", 4, 8)
        parts = [{key: group(n, v, var)} for n, v, var in [(10, 2, 1), (20, 3, 2), (30, 4, 3)]]
        direct = pool.pool_groups(parts)
        staged = pool.pool_groups([pool.pool_groups(parts[:2]), parts[2]])
        self.assertAlmostEqual(direct[key]["objects"]["isolated_all"].values[0],
                               staged[key]["objects"]["isolated_all"].values[0])
        self.assertAlmostEqual(direct[key]["objects"]["isolated_all"].variances[0],
                               staged[key]["objects"]["isolated_all"].variances[0])

    def test_loop_is_added_before_asymmetry_not_event_averaged_with_tree(self):
        base = {("photon", h, 4, 8): group(10000, 3 if h in ("PP", "MM") else 1, 1)
                for h in pool.HELICITIES}
        loop = {("photon", h, 4, 8): group(10, 1 if h in ("PP", "MM") else 3, 4)
                for h in pool.HELICITIES}
        combined = pool.add_loop_groups(base, [loop])
        prediction = pool.predict(combined)
        self.assertEqual(len(combined), 8)
        self.assertEqual(prediction["inclusive_cross_section"]["values"], [4])
        self.assertAlmostEqual(prediction["inclusive_cross_section"]["errors"][0], math.sqrt(5)/2)
        self.assertEqual(prediction["isolated_all"]["values"], [0])
        self.assertAlmostEqual(prediction["isolated_all"]["errors"][0], math.sqrt(5)/8)
        with self.assertRaisesRegex(ValueError, "already present"):
            pool.add_loop_groups(combined, [loop])

    def test_loop_repetitions_are_pooled_only_with_other_loop_repetitions(self):
        key = ("photon", "PP", 4, 8)
        result = pool.add_loop_groups({key: group(1000, 20, 1)},
                                     [{key: group(10, 2, 4)}, {key: group(30, 6, 16)}])
        self.assertEqual(result[("gammajet_loop", "PP", 4, 8)]["objects"]["isolated_all"].values, [5])
        self.assertEqual(result[key]["objects"]["isolated_all"].values, [20])
        with self.assertRaisesRegex(ValueError, "different channel"):
            pool.pool_groups([{key: group(10, 1, 1)}, {}])

    def test_loop_only_job_matrix_still_checks_complete_and_independent_shards(self):
        m = manifest()
        m["jobs"] = [j for j in m["jobs"] if j["channel"] == "photon"]
        pool.validate_jobs(m, set(), channels=("photon",))
        with self.assertRaisesRegex(ValueError, "Incomplete channel"):
            pool.validate_jobs(m, set())
        m["jobs"].pop()
        with self.assertRaisesRegex(ValueError, "Incomplete channel"):
            pool.validate_jobs(m, set(), channels=("photon",))

    def test_runtime_audit_is_bound_to_exact_artifacts(self):
        wrap = lambda h: {"runtime": {"provenance": {"artifacts": {"HwShower": {"sha256": h}}}}}
        audit = {"schema_version": 1, "baseline_artifacts": {"HwShower": "old"},
                 "loop_artifacts": {"HwShower": "new"}, "changed_artifacts": ["HwShower"],
                 "evidence": {"reason": "Recorded partner fallback replay",
                              "partner_replay": {"status": "PASS", "sha256": "record",
                                                 "all_compared_outputs_match": True}}}
        self.assertEqual(pool.validate_loop_runtime(wrap("old"), wrap("new"), audit), ["HwShower"])
        self.assertEqual(pool.validate_loop_runtime(wrap("old"), wrap("old"), None), [])
        for new, record in (("new", None), ("different", audit)):
            with self.assertRaises(ValueError):
                pool.validate_loop_runtime(wrap("old"), wrap(new), record)

    def test_generated_card_comparison_preserves_every_other_setting(self):
        obj = pool.LOOP_OBJECT
        tree = "read Common.in\ninsert ME 0 /Herwig/MatrixElements/MEGammaJet\nsaverun Tree EventGenerator\n"
        loop = (f"create Herwig::MEPP2GammaJetLoop {obj} HwMEHadron.so\n"
                f"set {obj}:AmplitudeMode Massive\nset {obj}:ScalePreFactor 1.0\n"
                "read Common.in\n" + f"insert ME 0 {obj}\nsaverun Loop EventGenerator\n")
        self.assertEqual(pool.card_commands(tree), pool.card_commands(loop, loop=True))
        self.assertNotEqual(pool.card_commands(tree), pool.card_commands(loop + "set Cut 5\n", loop=True))
        with self.assertRaisesRegex(ValueError, "Massive"):
            pool.card_commands(loop.replace("Massive", "Massless"), loop=True)

    def test_add_disjoint_slices_then_form_asymmetry(self):
        groups = {(c, h, lo, hi): group(100, 3 if h in ("PP", "MM") else 1, 1)
                  for c, h, (lo, hi) in pool.itertools.product(
                      ("photon", "fragmentation"), pool.HELICITIES, [(4, 8), (8, 12)])}
        prediction = pool.predict(groups)
        self.assertEqual(prediction["inclusive_cross_section"]["values"], [8])
        self.assertEqual(prediction["inclusive_cross_section"]["errors"], [1])
        self.assertEqual(prediction["isolated_all"]["values"], [0.5])
        self.assertAlmostEqual(prediction["isolated_all"]["errors"][0], math.sqrt(5)/16)

    def test_completed_and_independent_jobs_required(self):
        first = manifest()
        seen = set()
        pool.validate_jobs(first, seen)
        with self.assertRaisesRegex(ValueError, "Repeated final seed"):
            pool.validate_jobs(first, seen)
        for mutate in (lambda m: m.update(status="failed"),
                       lambda m: m["jobs"].pop(),
                       lambda m: m["jobs"][0].update(status="failed"),
                       lambda m: m["jobs"][0].update(events=99),
                       lambda m: m["jobs"][0].update(scale=2)):
            m = copy.deepcopy(first)
            mutate(m)
            with self.assertRaises(ValueError):
                pool.validate_jobs(m, set())

    def test_partial_yoda_and_nonstatistical_errors_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            job = manifest()["jobs"][0]
            job["output_yoda"] = "test.yoda"
            parts = ["BEGIN YODA_COUNTER_V3 /RAW/_EVTCOUNT\n---\n100 100 100\nEND YODA_COUNTER_V3\n"]
            for obj in pool.OBJECTS.values():
                parts.append(f"BEGIN YODA_ESTIMATE1D_V3 /{pool.MEASUREMENT}/{obj}\n"
                             '---\nEdges(A1): [6,7]\nErrorLabels: ["stats"]\n'
                             'nan --- ---\n2 -3 3\nnan --- ---\nEND YODA_ESTIMATE1D_V3\n')
            path = directory / "test.yoda"
            content = "\n".join(parts)
            path.write_text(content)
            _, count, objects = pool.read_shard((directory, job))
            self.assertEqual(count, 100)
            self.assertEqual(objects["isolated_all"].variances, [9])
            empty = content.replace('ErrorLabels: ["stats"]\n', '').replace('nan --- ---', 'nan').replace('2 -3 3', '0')
            path.write_text(empty)
            _, _, objects = pool.read_shard((directory, job))
            self.assertEqual(objects["isolated_all"].variances, [0])
            path.write_text(empty.replace('\n0\n', '\n1\n'))
            with self.assertRaisesRegex(ValueError, "without errors"):
                pool.read_shard((directory, job))
            path.write_text(content.replace("100 100 100", "99 99 99"))
            with self.assertRaisesRegex(ValueError, "Incomplete event count"):
                pool.read_shard((directory, job))
            path.write_text(content.replace('["stats"]', '["systematic"]'))
            with self.assertRaisesRegex(ValueError, "only stats"):
                pool.read_shard((directory, job))


if __name__ == "__main__":
    unittest.main()
