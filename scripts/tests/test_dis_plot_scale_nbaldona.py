"""Select saved scale predictions without substituting nominal values/errors."""
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
import test_dis_hide_scale_plot_nbaldona as legacy


class Bin:
    def __init__(self, value, stat=.12):
        self.value = value
        self.errors = {"stat": (-stat, stat), "pdf68": (-.2, .2), "hard_scale": (-.3, .3)}
    def sources(self): return list(self.errors)
    def rmSource(self, name): del self.errors[name]
    def val(self): return self.value
    def errAvg(self, name): return sum(abs(x) for x in self.errors[name]) / 2


class Estimate:
    def __init__(self, path, values=(-.27, float('nan'), -.08)):
        self.path = path
        self.annotations = {}
        self.contents = [Bin(x) for x in values]
        self.masked = {2}
    def clone(self): return copy.deepcopy(self)
    def setPath(self, path): self.path = path
    def setAnnotation(self, key, value): self.annotations[key] = value
    def bins(self, overflows=False, masked=False):
        assert masked  # Must preserve masked bins when constructing the error array.
        return self.contents
    def isMasked(self, index): return index in self.masked


class SelectedScaleTests(unittest.TestCase):
    def test_selects_correct_scale_paths_statistics_and_masks_for_both_targets(self):
        measurement = {"analysis": {"name": "E142"}, "outputs": {
            "He3_W050": {"a1": "/E142/d01", "apar": "/E142/aparHe3"},
            "Neutron_W050": {"a1": "/E142/d03", "apar": "/E142/aparN"},
        }}
        for scale in (.5, 1., 2.):
            objects = {}
            for selection, outputs in measurement["outputs"].items():
                for observable, target in outputs.items():
                    path = (target if scale == 1 else
                            f"/E142/VARIATIONS/p000-u000-mu{campaign._scale_token(scale)}/{observable}_{selection}")
                    objects[path] = Estimate(path)
            objects['/E142/VARIATIONS/p001-u000-mu2/a1_He3_W050'] = Estimate('/unwanted-replica')
            before = copy.deepcopy(objects)
            with mock.patch.object(campaign, '_import_yoda', return_value=types.SimpleNamespace(read=lambda _:objects)), \
                    mock.patch.object(campaign, '_write_yoda_objects') as write:
                bands = campaign._write_selected_scale_plot_copy(measurement, Path('input'), Path('copy'), scale)
            copies = write.call_args.args[1]
            self.assertEqual({obj.path for obj in copies}, {'/E142/d01','/E142/d03','/E142/aparHe3','/E142/aparN'})
            for obj in copies:
                self.assertEqual(obj.masked, {2})
                self.assertEqual(obj.contents[0].val(), -.27)
                self.assertEqual(obj.contents[0].sources(), ['stat'])
                self.assertEqual(obj.annotations['HardScaleFactor'], str(scale))
                self.assertEqual(bands[Path(obj.path).name], {'monte_carlo': [.12, None, .12]})
            for key, original in objects.items():
                self.assertEqual(original.path, before[key].path)
                self.assertEqual(original.contents[0].errors, before[key].contents[0].errors)
                self.assertEqual(original.annotations, {})

    def test_missing_scale_or_statistical_error_never_falls_back_to_central(self):
        measurement = {'analysis': {'name':'TEST'}, 'outputs': {'S': {'a1':'/TEST/d01'}}}
        path='/TEST/VARIATIONS/p000-u000-mu0p5/a1_S'
        for objects in ({'/TEST/d01': Estimate('/TEST/d01')}, {path: Estimate(path)}):
            if path in objects: del objects[path].contents[0].errors['stat']
            with mock.patch.object(campaign, '_import_yoda', return_value=types.SimpleNamespace(read=lambda _:objects)), \
                    mock.patch.object(campaign, '_write_yoda_objects') as write, \
                    self.assertRaises(campaign.CampaignError):
                campaign._write_selected_scale_plot_copy(measurement, Path('input'), Path('copy'), .5)
            write.assert_not_called()

    def fixture(self, root):
        measurement, prediction = legacy.HideScalePlotTests().fixture(root)
        measurement['analysis']['name'] = 'TEST'
        path=root/campaign.MANIFEST_NAME
        manifest=json.loads(path.read_text())
        manifest['configuration']={'variation_points':[[0,0,1.],[0,0,.5],[0,0,2.]]}
        manifest['postprocess']['variations_yoda']='postprocess/variations.yoda'
        (root/'postprocess/variations.yoda').write_text('saved variation values and errors')
        campaign.atomic_write_json(path,manifest)
        return measurement,prediction

    def test_selected_plot_pipeline_uses_own_errors_and_preserves_full_results(self):
        for scale in (.5,1.,2.):
            with self.subTest(scale=scale), tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp); measurement,prediction=self.fixture(root)
                token=campaign._scale_token(scale); output=root/f'plots-mu{token}'
                args=campaign.make_parser().parse_args(['plot','--measurement','TEST','--tag','test','--plot-scale',str(scale)])
                before={p:p.read_bytes() for p in (root/'postprocess').iterdir()}
                def generate(command,*unused):
                    self.assertIn(str(output),command)
                    self.assertTrue(any(str(root/'work'/f'plot-mu{token}'/'nominal.yoda') in v and f'={scale:g}$' in v for v in command))
                    output.mkdir(exist_ok=True);(output/'d01.py').write_text("# fake script\nplt.savefig('test.png')\n");(output/'d01.png').write_bytes(b'plot')
                def index(*unused, **kwargs):
                    result=output/'index.html';result.write_text('index');return result
                with contextlib.ExitStack() as stack:
                    stack.enter_context(mock.patch.object(campaign,'DISPOL_ROOT',root))
                    for name,replacement in [
                        ('campaign_directory',mock.Mock(return_value=root)),
                        ('_assert_manifest_signatures_current',mock.Mock(return_value=None)),
                        ('ensure_reference_yoda', mock.Mock()),
                        ('analysis_environment',mock.Mock(return_value={})),
                        ('write_campaign_monitor_files',mock.Mock()),
                        ('build_campaign_monitor_payload',mock.Mock()),
                        ('existing_herwig_monitor',mock.Mock()),
                        ('_run_logged',mock.Mock(side_effect=generate)),
                        ('add_experimental_error_overlay',mock.Mock()),
                        ('plot_script_has_finite_y',mock.Mock(return_value=True)),
                        ('write_plot_indexes',mock.Mock(side_effect=index)),
                    ]:stack.enter_context(mock.patch.object(campaign,name,replacement))
                    selected=stack.enter_context(mock.patch.object(campaign,'_write_selected_scale_plot_copy',return_value={'d01':{'monte_carlo':[.12]}}))
                    overlay=stack.enter_context(mock.patch.object(campaign,'add_theory_uncertainty_overlay'))
                    stack.enter_context(mock.patch.object(campaign.subprocess,'run',return_value=types.SimpleNamespace(returncode=0)))
                    self.assertEqual(campaign.plot_campaign(args,measurement),output)
                source=prediction if scale==1 else root/'postprocess/variations.yoda'
                self.assertEqual(selected.call_args.args[1],source)
                self.assertEqual(overlay.call_args.args[1],{'monte_carlo':[.12]})  # NOT nominal .03!
                self.assertEqual(overlay.call_args.args[2],'nominal.yoda')
                self.assertEqual(before,{p:p.read_bytes() for p in (root/'postprocess').iterdir()})
                self.assertEqual((root/'plots/old.png').read_bytes(),b'original full plot')
                updated=json.loads((root/campaign.MANIFEST_NAME).read_text())
                self.assertEqual(updated['plots'],{'index':'plots/index.html'})
                self.assertEqual(updated[f'plots_mu{token}']['hard_scale_factor'],scale)
                self.assertEqual(updated[f'plots_mu{token}']['uncertainty_mode'],'mc-stat-only')

    def test_dry_run_no_yoda_no_writes(self):
        for scale in (.5,1.,2.):
            with tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);measurement,_=self.fixture(root)
                before={p:p.read_bytes() for p in root.rglob('*') if p.is_file()}
                args=campaign.make_parser().parse_args(['plot','--measurement','TEST','--tag','test','--plot-scale',str(scale),'--dry-run'])
                with mock.patch.object(campaign,'campaign_directory',return_value=root), \
                        mock.patch.object(campaign,'_assert_manifest_signatures_current',return_value=None), \
                        mock.patch.object(campaign,'_import_yoda') as yoda,contextlib.redirect_stdout(io.StringIO()) as output:
                    self.assertEqual(campaign.plot_campaign(args,measurement),root/f'plots-mu{campaign._scale_token(scale)}')
                yoda.assert_not_called();self.assertIn('MC statistical only:',output.getvalue())
                self.assertEqual(before,{p:p.read_bytes() for p in root.rglob('*') if p.is_file()})

    def test_invalid_combinations_and_absent_scale_are_rejected(self):
        with contextlib.redirect_stderr(io.StringIO()),self.assertRaises(SystemExit):
            campaign.make_parser().parse_args(['plot','--measurement','TEST','--tag','test','--plot-scale','2','--hide-scale-uncertainty'])
        for case in ('comparison','absent-scale','missing-file'):
            with tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);measurement,_=self.fixture(root)
                args=campaign.make_parser().parse_args(['plot','--measurement','TEST','--tag','test','--plot-scale','2','--dry-run']+(['--plot-comparisons'] if case=='comparison' else []))
                if case=='absent-scale':
                    p=root/campaign.MANIFEST_NAME;m=json.loads(p.read_text());m['configuration']['variation_points']=[[0,0,1.]];p.write_text(json.dumps(m))
                if case=='missing-file':(root/'postprocess/variations.yoda').unlink()
                with mock.patch.object(campaign,'campaign_directory',return_value=root), \
                        mock.patch.object(campaign,'_assert_manifest_signatures_current',return_value=None),self.assertRaises(campaign.CampaignError):
                    campaign.plot_campaign(args,measurement)

    def test_native_yoda_selected_scale_roundtrip(self):
        try:yoda=campaign._import_yoda()
        except campaign.CampaignError:self.skipTest('Native YODA unavailable locally')
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'variations.yoda';destination=root/'selected.yoda'
            path='/TEST/VARIATIONS/p000-u000-mu0p5/a1_S'
            obj=campaign._estimate_from_values(yoda,[0.,1.,2.],path,[-.27,None],[.12,None],{})
            obj.bin(1).setErr(-.3,.4,'hard_scale')
            campaign._write_yoda_objects(yoda,[obj],source);before=source.read_bytes()
            bands=campaign._write_selected_scale_plot_copy({'analysis':{'name':'TEST'},'outputs':{'S':{'a1':'/TEST/d01'}}},source,destination,.5)
            result=yoda.read(str(destination))['/TEST/d01']
            self.assertAlmostEqual(result.bin(1).val(),-.27)
            self.assertAlmostEqual(result.bin(1).errAvg('stat'),.12)
            self.assertEqual(set(result.bin(1).sources()),{'stat'})
            self.assertTrue(result.isMasked(2));self.assertEqual(source.read_bytes(),before)
            self.assertEqual(bands,{'d01':{'monte_carlo':[.12,None]}})


if __name__=='__main__':unittest.main()
