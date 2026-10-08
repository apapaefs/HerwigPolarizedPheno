"""Additive ratio publication preserves prior campaigns and verifies complete caches."""
from pathlib import Path
import copy
import hashlib
import json
import re
import sys
import tempfile
import unittest
from unittest import mock

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import hermes_data_ratio_gallery as gallery
import hermes_born_ratio_plots as born
import hermes_projection_ratio_plots as projections
import run_experimental_campaign as campaign
from scripts.tests.test_hermes_born_cell_plots import fixtures
from scripts.tests.test_hermes_projection_ratio_plots import inputs


class RatioGalleryTests(unittest.TestCase):
    def setUp(self):
        tmp=tempfile.TemporaryDirectory();self.addCleanup(tmp.cleanup)
        self.root=Path(tmp.name);self.campaign=self.root/'campaign';self.output=self.campaign/'plots'
        self.analysis=self.output/gallery.MEASUREMENT
        self.analysis.mkdir(parents=True);(self.campaign/'postprocess').mkdir()
        ref,summary=fixtures();projected,aref,integrated=inputs()
        summary.update(tag='fixture',acceptance='rectangle_intersect_polar_ring')
        summary['bins']+=projected['bins']
        self.summary=self.campaign/'postprocess/summary.json';self.summary.write_text(json.dumps(summary))
        self.bref=self.root/'born.json';self.bref.write_text(json.dumps(ref))
        self.aref=self.root/'a1.json';self.aref.write_text(json.dumps(aref))
        integrated['input_files']={'summary':{'sha256':gallery.digest(self.summary)},'reference':{'sha256':gallery.digest(self.bref)}}
        self.integrated=self.campaign/'integrated.json';self.integrated.write_text(json.dumps(integrated))
        self.addCleanup(mock.patch.stopall)
        mock.patch.object(gallery,'REFERENCE',self.aref).start()
        mock.patch.object(gallery,'BORN_REFERENCE',self.bref).start()
        mock.patch.object(born,'render',side_effect=lambda snapshot,output:self.fake(output,born.REQUIRED_OUTPUTS)).start()
        mock.patch.object(projections,'render',side_effect=lambda snapshot,output:self.fake(output,projections.REQUIRED_OUTPUTS)).start()
        (self.analysis/'old.png').write_bytes(b'preserved absolute figure')
        (self.analysis/'index.html').write_text('<html><body><img src="old.png">original analysis</body></html>')
        (self.output/'index.html').write_text('<html><body><a href="HERMES_2007_I726689/index.html">original root</a></body></html>')

    @staticmethod
    def fake(output,names):
        result=[]
        for name in sorted(names):
            path=output/name;path.write_bytes(b'valid fixture plot');result.append(path)
        return result

    def ensure(self):return gallery.ensure_campaign_plots(self.campaign,self.output,self.integrated)

    def test_complete_reuse_preserves_original_figures_and_numerical_inputs(self):
        protected={p:p.read_bytes() for p in (self.summary,self.integrated,self.analysis/'old.png')}
        first=self.ensure();second=self.ensure()
        self.assertTrue(first['created']);self.assertFalse(second['created'])
        self.assertEqual(first['directory'],second['directory'])
        self.assertEqual(len(gallery.expected_outputs()),127)
        for path,payload in protected.items():self.assertEqual(path.read_bytes(),payload)
        direct=json.loads((self.campaign/first['directory']/'ratios.json').read_text(encoding="utf-8"))
        self.assertEqual(len(direct['born_cells']['targets']['P']['cells']),45)

    def test_damaged_output_is_retained_and_rebuilt_in_a_new_revision(self):
        first=self.ensure();old=self.campaign/first['directory']
        damaged=old/'A1_P_Q2GT1_ratio.pdf';damaged.write_bytes(b'corrupt')
        second=self.ensure()
        self.assertNotEqual(first['directory'],second['directory'])
        self.assertEqual(damaged.read_bytes(),b'corrupt')
        self.assertTrue(gallery._valid_cache(self.campaign/second['directory'],second['cache_key']))

    def test_manifest_cannot_silently_drop_a_required_output(self):
        first=self.ensure();directory=self.campaign/first['directory']
        path=directory/gallery.MANIFEST;record=json.loads(path.read_text(encoding="utf-8"))
        name='A1_P_Q2GT1_ratio.pdf';record['outputs'].pop(name);(directory/name).unlink()
        path.write_text(json.dumps(record))
        self.assertFalse(gallery._valid_cache(directory,first['cache_key']))

    def test_stale_or_wrong_campaign_reconstruction_is_rejected(self):
        original=json.loads(self.integrated.read_text(encoding="utf-8"))
        for kind in ('summary','reference','campaign'):
            item=copy.deepcopy(original)
            if kind=='campaign':item['campaign_tag']='wrong'
            else:item['input_files'][kind]['sha256']='wrong'
            self.integrated.write_text(json.dumps(item))
            with self.subTest(kind=kind),self.assertRaises(ValueError):self.ensure()

    def test_additive_indexes_are_idempotent_and_portable(self):
        info=self.ensure();gallery.append_to_existing_indexes(self.campaign,self.output,info)
        original={p:p.read_bytes() for p in (self.output/'index.html',self.analysis/'index.html')}
        gallery.append_to_existing_indexes(self.campaign,self.output,info)
        for path,payload in original.items():self.assertEqual(path.read_bytes(),payload)
        self.assertEqual((self.analysis/'old.png').read_bytes(),b'preserved absolute figure')
        for page in self.output.rglob('*.html'):
            for link in re.findall(r'(?:href|src)="([^"]+)"',page.read_text(encoding="utf-8")):
                resolved=(page.parent/link).resolve()
                self.assertTrue(resolved.is_file(),(page,link))
                boundary=self.analysis if page.is_relative_to(self.analysis) else self.output
                self.assertTrue(resolved.is_relative_to(boundary.resolve()),(page,link))
        self.assertIn('original analysis',(self.analysis/'index.html').read_text(encoding="utf-8"))
        self.assertIn('original root',(self.output/'index.html').read_text(encoding="utf-8"))

    def test_both_a1_figures_are_embedded_with_original_plot_identifiers(self):
        for rivet_id in ('d14-x01-y01','d14-x01-y02'):
            (self.analysis/f'{rivet_id}.pdf').write_bytes(b'preserved original PDF')
        info=self.ensure();gallery.append_to_existing_indexes(self.campaign,self.output,info)
        for page in (self.output/'index.html',self.analysis/'index.html'):
            text=page.read_text(encoding='utf-8')
            for target,rivet_id in (('P','d14-x01-y01'),('D','d14-x01-y02')):
                self.assertEqual(text.count(f'id="{rivet_id}-ratio"'),1)
                self.assertRegex(text,rf'<img src="[^"]*/A1_{target}_Q2GT1_ratio.png"')
                self.assertIn(f'A1_{target}_Q2GT1_ratio_fullrange.pdf',text)
                self.assertIn(f'{rivet_id}.pdf',text)
            self.assertEqual(text.count('Original absolute plot'),2)
        for rivet_id in ('d14-x01-y01','d14-x01-y02'):
            self.assertEqual((self.analysis/f'{rivet_id}.pdf').read_bytes(),b'preserved original PDF')

    def test_auto_plot_hook_is_HERMES_only_and_uses_current_reconstruction(self):
        with mock.patch.object(gallery,'ensure_campaign_plots',return_value={'ok':True}) as ensure:
            self.assertIsNone(campaign._hermes_data_ratio_plot_gallery(self.campaign,self.output,{'id':'OTHER'},{}))
            self.assertIsNone(campaign._hermes_data_ratio_plot_gallery(self.campaign,self.output,{'id':gallery.MEASUREMENT},None))
            info=campaign._hermes_data_ratio_plot_gallery(self.campaign,self.output,{'id':gallery.MEASUREMENT},{'directory':'derived/current'})
            self.assertEqual(info,{'ok':True})
            ensure.assert_called_once_with(self.campaign,self.output,self.campaign/'derived/current/integrated.json')

    def test_normal_gallery_writer_embeds_ratios_without_touching_existing_plot(self):
        info=self.ensure();script=self.analysis/'old.py';script.write_text('# original')
        campaign.write_plot_indexes(self.output,{'id':gallery.MEASUREMENT},[script],data_ratios=info)
        self.assertIn('Additional MC/data ratio figures',(self.analysis/'index.html').read_text(encoding="utf-8"))
        self.assertEqual((self.analysis/'old.png').read_bytes(),b'preserved absolute figure')


if __name__=='__main__':unittest.main()
