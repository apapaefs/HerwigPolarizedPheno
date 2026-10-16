"""Saved-campaign checks must be strict, independent and non-destructive."""
import ast
import contextlib
import copy
import io
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import check_he3_saved_campaign as check
import run_experimental_campaign as campaign
from rivet_mkhtml_safe import RESIDUAL_LABEL


class SavedCampaignCheckTests(unittest.TestCase):
    def make_source(self, root):
        measurement = campaign.get_measurement('E154_1997_I443170')
        source = root / 'source'
        source.mkdir()
        jobs = campaign.build_job_matrix(measurement, 1, 1, 1, 123)
        for job in jobs:
            path = source / job['output_yoda']
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text('a deliberately tiny saved input')
            job['status'] = 'success'
        manifest = {'measurement': measurement['id'], 'tag': 'test', 'jobs': jobs,
                    'configuration': {'shards': 1, 'measurement_signature': campaign.measurement_signature(measurement)},
                    'runtime': {'tools': {'rivet-mkhtml': '/unused/original'}}, 'history': [],
                    'postprocess': {'summary_json': 'postprocess/summary.json', 'yoda': 'postprocess/prediction.yoda'}}
        campaign.atomic_write_json(source / 'manifest.json', manifest)
        campaign.atomic_write_json(source / 'resolved-measurement.json', measurement)
        campaign.atomic_write_json(source / 'postprocess/summary.json', {'order_sign_audit': {'order_coefficients': {'POSNLO': 1, 'NEGNLO': -1}}})
        (source / 'postprocess/prediction.yoda').write_text('baseline')
        return source

    def test_plan_is_read_only_and_accounts_for_all_shards(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = self.make_source(root)
            before = {p: p.read_bytes() for p in source.rglob('*') if p.is_file()}
            with contextlib.redirect_stdout(io.StringIO()) as stream, mock.patch.object(campaign, '_import_yoda') as native:
                self.assertEqual(check.main(['--campaign', str(source), '--output-root', str(root/'check'), '--dry-run']), 0)
            plan = json.loads(stream.getvalue())[0]
            self.assertEqual(plan['jobs'], 8)
            self.assertEqual(plan['input_files'], 12)
            native.assert_not_called()
            self.assertFalse((root/'check').exists())
            self.assertEqual(before, {p: p.read_bytes() for p in source.rglob('*') if p.is_file()})

    def test_path_escape_existing_output_and_nested_output_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = self.make_source(root)
            (source/'escape').symlink_to(root/'outside')
            (root/'outside').write_text('outside')
            for path in ('../outside', str(root/'outside'), 'escape'):
                with self.assertRaises(campaign.CampaignError):
                    check.contained_file(source, path)
            with self.assertRaises(campaign.CampaignError):
                check.inspect(source, source/'nested-check')
            with contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(check.main(['--campaign',str(source),'--output-root',str(root)]), 2)

    def test_failed_reprocessing_preserves_original_and_records_failure(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = self.make_source(root)
            plan = check.inspect(source, root/'validation')
            before = {p: p.read_bytes() for p in source.rglob('*') if p.is_file()}
            original_campaign_root = campaign.CAMPAIGN_ROOT
            tool = root/'rivet-mkhtml'
            tool.write_text('# test tool')
            with mock.patch.object(campaign, 'postprocess_campaign', side_effect=campaign.CampaignError('intentional failure')), \
                    contextlib.redirect_stdout(io.StringIO()), self.assertRaisesRegex(campaign.CampaignError, 'intentional failure'):
                check.run_check(plan, tool)
            self.assertEqual(campaign.CAMPAIGN_ROOT, original_campaign_root)
            self.assertEqual(before, {p: p.read_bytes() for p in source.rglob('*') if p.is_file()})
            report = campaign.load_json(plan['destination']/'validation.json')
            self.assertEqual(report['status'], 'failed')
            self.assertTrue(report['original_inputs_unchanged'])
            copied = plan['destination']/plan['manifest']['jobs'][0]['output_yoda']
            self.assertFalse(copied.is_symlink())
            self.assertNotEqual(copied.stat().st_ino, (source/plan['manifest']['jobs'][0]['output_yoda']).stat().st_ino)
            shim = plan['destination']/'validation-controls/rivet-mkhtml-validation.py'
            compile(shim.read_text(), str(shim), 'exec')
            # Rivet filters both MC and reference YODA with this same expression.
            options = ast.literal_eval(ast.parse(shim.read_text()).body[1].value)
            pattern = options[options.index('--match') + 1]
            for output in plan['measurement']['outputs'].values():
                self.assertRegex(output['a1'], pattern)
                self.assertRegex('/REF' + output['a1'], pattern)
                self.assertNotRegex('/RAW' + output['a1'], pattern)
                self.assertNotRegex(output['a1'] + '_unrelated', pattern)


    def test_scientific_comparison_rejects_changed_errors_values_masks_and_shape(self):
        check.compare({'value': [1., None, float('nan')]}, {'value': [1.+1e-12, None, float('nan')]})
        for expected, actual in (({'stat': .1}, {'stat': .2}), ([None], [0.]), ([1.], [1.,2.]), ({'value':1}, {'other':1})):
            with self.assertRaises(campaign.CampaignError):
                check.compare(expected, actual)

    def summary_fixture(self):
        row = {'selection': 'He3_W050', 'bin': 1, 'x_low': .03, 'x_high': .04,
               'a1': -.02, 'a1_stat': .0006, 'a_parallel': -.013,
               'a_parallel_stat': .0004, 'sigma_uu_pb': 72000., 'sigma_ll_pb': -970.}
        old = {'bins': [row], 'uncertainties': {'mc': [.0006]}, 'masked_bins': [],
               'diagnostics': {'accepted': [10.]}, 'variation_points': [[0, 0, 1.]]}
        new = copy.deepcopy(old)
        new['bins'][0].update(axis='x', bin_low=.03, bin_high=.04, supported=True,
                              acceptance='published-mean-Q2 window')
        return old, new

    def test_summary_schema_upgrade_preserves_physics_and_inputs(self):
        old, new = self.summary_fixture()
        originals = copy.deepcopy((old, new))
        check.compare(check.summary_science(old), check.summary_science(new))
        self.assertEqual((old, new), originals)
        new['bins'][0]['acceptance'] = 'expanded explanatory wording'
        check.compare(check.summary_science(old), check.summary_science(new))

    def test_summary_schema_upgrade_rejects_changed_or_unknown_scientific_fields(self):
        old, new = self.summary_fixture()
        for changes in ({'a1': -.03}, {'a1_stat': .001}, {'a1': None},
                        {'x_low': .02}, {'bin_low': .02}, {'axis': 'q2'},
                        {'supported': False}, {'supported': 1}, {'extra_physics': 1.}):
            with self.subTest(changes=changes), self.assertRaises(campaign.CampaignError):
                changed = copy.deepcopy(new)
                changed['bins'][0].update(changes)
                check.compare(check.summary_science(old), check.summary_science(changed))
        changed = copy.deepcopy(new)
        del changed['bins'][0]['a1_stat']
        with self.assertRaises(campaign.CampaignError):
            check.compare(check.summary_science(old), check.summary_science(changed))
        changed = copy.deepcopy(new)
        changed['masked_bins'] = [1]
        with self.assertRaises(campaign.CampaignError):
            check.compare(check.summary_science(old), check.summary_science(changed))

    def test_numeric_reader_does_not_execute_generated_code(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary); path=root/'data.py'
            path.write_text("x = {'a':[1., -2., nan, inf]}\nraise RuntimeError('must not execute')\n")
            result=check.numeric_assignments(path)['x']['a']
            self.assertEqual(result[:2], [1.,-2.])
            self.assertTrue(math.isnan(result[2])); self.assertTrue(math.isinf(result[3]))

    def test_plot_checks_detect_wrong_curve_errors_residual_or_label(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary); script=root/'plot.py'
            good_script = (f"ax.set_ylabel({RESIDUAL_LABEL!r})\n"
                           "_theory_mc = [0.1]\n_theory_pdf = None\n_theory_scale_dn = None\n_theory_scale_up = None\n")
            good_data = "yvals={'Data':[.2], 'MC':[.3]}\nyerrs={'Data':[[.05],[.05]]}\nratio0_yvals={'Data':[0.], 'MC':[2.]}\n"
            for suffix in ('.png','.pdf'): script.with_suffix(suffix).write_bytes(b'fixture')
            data=script.with_name('plot__data.py')
            expected={'edges':[0.,1.], 'bins':[{'value':.3,'errors':{'stat':[-.1,.1]}}]}
            reference=[{'value':.2}]; bands={'monte_carlo':[.1]}
            script.write_text(good_script); data.write_text(good_data)
            check.check_plot(script, expected, reference, bands)
            variants=[(good_script,good_data.replace("'MC':[2.]", "'MC':[1.5]")),
                      (good_script.replace('[0.1]', '[0.2]'),good_data),
                      (good_script.replace(repr(RESIDUAL_LABEL),repr('MC/data')),good_data),
                      (good_script,good_data.replace("'MC':[.3]", "'MC':[.4]"))]
            for source,numbers in variants:
                script.write_text(source);data.write_text(numbers)
                with self.assertRaises(campaign.CampaignError):
                    check.check_plot(script,expected,reference,bands)


if __name__ == '__main__':
    unittest.main()
