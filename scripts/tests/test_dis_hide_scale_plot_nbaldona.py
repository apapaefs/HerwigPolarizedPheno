"""Presentation-only scale filtering, including preservation of full results."""
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

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import run_experimental_campaign as campaign


class Bin:
    def __init__(self, errors):
        self.errors = dict(errors)
    def sources(self):
        return list(self.errors)
    def rmSource(self, name):
        del self.errors[name]


class Estimate:
    def __init__(self, bins):
        self.contents = bins
    def clone(self):
        return copy.deepcopy(self)
    def bins(self, overflows, masked):
        assert overflows and masked
        return self.contents


class HideScalePlotTests(unittest.TestCase):
    def test_filter_keeps_stat_exactly_and_does_not_change_source(self):
        original = Estimate([Bin({}), Bin({'stat': (-.03, .04),
                              'hard_scale': (-.3, .4), 'pdf68': (-.1, .1)}),
                             Bin({'stat': (-.01, .01)})])
        yoda = types.SimpleNamespace(read=lambda _: {'/test': original})
        with mock.patch.object(campaign, '_import_yoda', return_value=yoda), \
             mock.patch.object(campaign, '_write_yoda_objects') as write:
            campaign._write_no_scale_plot_copy(Path('source.yoda'), Path('copy.yoda'))
        result = write.call_args.args[1][0]
        self.assertEqual(result.contents[1].errors, {'stat': (-.03, .04), 'pdf68': (-.1, .1)})
        self.assertEqual(result.contents[2].errors, {'stat': (-.01, .01)})
        self.assertEqual(len(original.contents[1].errors), 3)

    def test_no_scale_input_remains_unchanged(self):
        original = Estimate([Bin({'stat': (-.03, .04), 'pdf68': (-.1, .1)})])
        yoda = types.SimpleNamespace(read=lambda _: {'/test': original})
        with mock.patch.object(campaign, '_import_yoda', return_value=yoda), \
             mock.patch.object(campaign, '_write_yoda_objects') as write:
            campaign._write_no_scale_plot_copy(Path('source.yoda'), Path('copy.yoda'))
        self.assertEqual(write.call_args.args[1][0].contents[0].errors,
                         original.contents[0].errors)

    def fixture(self, root):
        prediction = root/'postprocess'/'nominal.yoda'
        prediction.parent.mkdir()
        prediction.write_text('original full uncertainty input')
        (prediction.parent/'summary.json').write_text(json.dumps({'uncertainties': {
            'nominal': {'a1': {'monte_carlo': [.03], 'pdf_68': [.1],
                              'scale_down': [-.2], 'scale_up': [.2]}}}}))
        (root/'reference.json').write_text('{}')
        (root/'metadata.plot').write_text('metadata')
        (root/'logs').mkdir()
        (root/'plots').mkdir()
        (root/'plots'/'old.png').write_bytes(b'original full plot')
        manifest = {'runtime': {'tools': {'rivet-mkhtml': '/test/rivet-mkhtml'}},
                    'history': [], 'plots': {'index': 'plots/index.html'},
                    'postprocess': {'predictions': [{'family': 'nominal', 'label': 'Theory',
                                                    'yoda': 'postprocess/nominal.yoda'}]}}
        campaign.atomic_write_json(root/campaign.MANIFEST_NAME, manifest)
        measurement = {'id': 'TEST', 'reference': {'snapshot': str(root/'reference.json')},
                       'analysis': {'plot': str(root/'metadata.plot')},
                       'outputs': {'nominal': {'a1': '/TEST/d01'}}}
        return measurement, prediction

    def test_plot_only_pipeline_and_default_keep_separate_outputs(self):
        for stat_only in (False, True):
            with self.subTest(stat_only=stat_only), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                measurement, prediction = self.fixture(root)
                output = root/('plots-central' if stat_only else 'plots')
                args = campaign.make_parser().parse_args(['plot', '--measurement', 'TEST',
                    '--tag', 'test'] + (['--hide-scale-uncertainty'] if stat_only else []))
                def generate(command, *unused):
                    self.assertIn(str(output), command)
                    expected = root/'work'/'plot-central'/'nominal.yoda' if stat_only else prediction
                    self.assertIn(str(expected)+':Title=Theory', command)
                    output.mkdir(exist_ok=True)
                    (output/'d01.py').write_text("# fake plot script\nplt.savefig('test.png')\n")
                    (output/'d01.png').write_bytes(b'new plot')
                def index(*unused, **kwargs):
                    result = output/'index.html'
                    result.write_text('index')
                    return result
                with contextlib.ExitStack() as stack:
                    stack.enter_context(mock.patch.object(campaign, "DISPOL_ROOT", root))
                    for name, replacement in [
                        ('campaign_directory', mock.Mock(return_value=root)),
                        ('_assert_manifest_signatures_current', mock.Mock(return_value=None)),
                        ('ensure_reference_yoda', mock.Mock()),
                        ('analysis_environment', mock.Mock(return_value={})),
                        ('write_campaign_monitor_files', mock.Mock()),
                        ('build_campaign_monitor_payload', mock.Mock()),
                        ('existing_herwig_monitor', mock.Mock()),
                        ('_run_logged', mock.Mock(side_effect=generate)),
                        ('add_experimental_error_overlay', mock.Mock()),
                        ('plot_script_has_finite_y', mock.Mock(return_value=True)),
                        ('write_plot_indexes', mock.Mock(side_effect=index)),
                    ]:
                        stack.enter_context(mock.patch.object(campaign, name, replacement))
                    filtered = stack.enter_context(mock.patch.object(campaign, '_write_no_scale_plot_copy'))
                    overlay = stack.enter_context(mock.patch.object(campaign, 'add_theory_uncertainty_overlay'))
                    stack.enter_context(mock.patch.object(campaign.subprocess, 'run', return_value=types.SimpleNamespace(returncode=0)))
                    self.assertEqual(campaign.plot_campaign(args, measurement), output)
                bands = overlay.call_args.args[1]
                self.assertEqual(bands['monte_carlo'], [.03])
                self.assertEqual('scale_up' in bands, not stat_only)
                self.assertEqual(bands['pdf_68'], [.1])
                self.assertEqual(filtered.call_count, int(stat_only))
                self.assertEqual(prediction.read_text(), 'original full uncertainty input')
                updated = json.loads((root/campaign.MANIFEST_NAME).read_text())
                if stat_only:
                    self.assertEqual(updated['plots'], {'index': 'plots/index.html'})
                    self.assertEqual(updated['plots_central']['index'], 'plots-central/index.html')
                    self.assertEqual((root/'plots'/'old.png').read_bytes(), b'original full plot')

    def test_dry_run_no_writes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            measurement, _ = self.fixture(root)
            before = {p: p.read_bytes() for p in root.rglob('*') if p.is_file()}
            args = campaign.make_parser().parse_args(['plot', '--measurement', 'TEST',
                '--tag', 'test', '--hide-scale-uncertainty', '--dry-run'])
            with mock.patch.object(campaign, 'campaign_directory', return_value=root), \
                 mock.patch.object(campaign, '_assert_manifest_signatures_current', return_value=None), \
                 mock.patch.object(campaign, '_import_yoda') as yoda, \
                 contextlib.redirect_stdout(io.StringIO()) as captured:
                self.assertEqual(campaign.plot_campaign(args, measurement), root/'plots-central')
                self.assertIn('Scale-free plotting copy:', captured.getvalue())
                yoda.assert_not_called()
            self.assertEqual(before, {p: p.read_bytes() for p in root.rglob('*') if p.is_file()})

    def test_native_yoda_roundtrip(self):
        try:
            yoda = campaign._import_yoda()
        except campaign.CampaignError:
            self.skipTest('Native YODA unavailable locally')
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source, destination = root/'source.yoda', root/'copy.yoda'
            obj = campaign._estimate_from_values(yoda, [0., 1., 2.], '/TEST/d01',
                                                  [.1, None], [.03, None], {})
            obj.bin(1).setErr(-.2, .4, 'hard_scale')
            obj.bin(1).setErr(-.1, .1, 'pdf68')
            campaign._write_yoda_objects(yoda, [obj], source)
            before = source.read_bytes()
            campaign._write_no_scale_plot_copy(source, destination)
            result = yoda.read(str(destination))['/TEST/d01']
            self.assertEqual(set(result.bin(1).sources()), {'stat', 'pdf68'})
            self.assertAlmostEqual(result.bin(1).val(), .1)
            self.assertAlmostEqual(result.bin(1).errUp('stat'), .03)
            self.assertAlmostEqual(result.bin(1).errUp('pdf68'), .1)
            self.assertTrue(result.isMasked(2))
            self.assertEqual(before, source.read_bytes())


if __name__ == '__main__':
    unittest.main()
