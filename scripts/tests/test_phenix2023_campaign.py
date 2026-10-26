"""Independent event classes must be added before any spin ratio is formed."""
import hashlib
import argparse
from collections import defaultdict
import contextlib
import copy
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import run_phenomenology_campaign as campaign
import run_experimental_campaign as experimental
import run_PHENIX_2023_I2033856 as wrapper
import run_PHENIX_2023_I2033856_MEGammaJet_Only as baseline_wrapper

ID='PHENIX_2023_I2033856'

def objects(vals, variance=1., helicities=('PP','PM','MP','MM')):
    return {name:{h:experimental.BinSeries([6.,7.],[v],[variance])
                 for h,v in zip(helicities,vals)}
            for name in ('inclusive_cross_section','isolated_cross_section','isolated_all')}

class PHENIXCampaignTests(unittest.TestCase):
    def setUp(self):
        self.m=campaign.discover_pp_registry()[ID]

    def test_sum_cross_sections_before_ratio_with_independent_errors(self):
        direct=('nominal','photon',0,0,1.,'on')
        frag=('nominal','fragmentation',0,0,1.,'on')
        raw={direct:objects([12.,8.,8.,12.],2.),frag:objects([1.,9.,9.,1.],3.)}
        summed,pred=campaign._combine_phenix_channels(self.m,raw)
        # Direct A=.2; fragmentation A=-.8; the combined A is -8/60,
        # not the unweighted average -.3. Variances add before propagation.
        self.assertAlmostEqual(pred[direct]['isolated_all']['values'][0],-8./60.)
        self.assertAlmostEqual(pred[direct]['inclusive_cross_section']['values'][0],15.)
        self.assertAlmostEqual(pred[direct]['inclusive_cross_section']['errors'][0]**2,1.25)
        self.assertEqual(summed[direct]['isolated_all']['PP'].variances,[5.])
        gradient=[(s*60.+8.)/60.**2 for s in (1,-1,-1,1)]
        expected=5.*sum(g*g for g in gradient)
        self.assertAlmostEqual(pred[direct]['isolated_all']['errors'][0]**2,expected)
        self.assertAlmostEqual(pred[direct]['Component_photon_isolated_all']['values'][0],.2)
        self.assertAlmostEqual(pred[direct]['Component_fragmentation_isolated_all']['values'][0],-.8)
        self.assertEqual(raw[direct]['isolated_all']['PP'].values,[12.])

    def test_zero_total_yield_is_masked(self):
        a=('nominal','photon',0,0,1.,'on'); b=('nominal','fragmentation',0,0,1.,'on')
        _,pred=campaign._combine_phenix_channels(self.m,{a:objects([0,0,0,0],0.),b:objects([0,0,0,0],0.)})
        self.assertEqual(pred[a]['isolated_all']['values'],[None])
        self.assertEqual(pred[a]['isolated_all']['errors'],[None])

    def test_missing_component_and_mismatched_bins_are_rejected(self):
        a=('nominal','photon',0,0,1.,'on'); b=('nominal','fragmentation',0,0,1.,'on')
        with self.assertRaisesRegex(campaign.CampaignError,'Missing additive'):
            campaign._combine_phenix_channels(self.m,{a:objects([1,1,1,1])})
        raw={a:objects([1,1,1,1]),b:objects([1,1,1,1])}
        raw[b]['isolated_all']['PP'].edges=[6.,8.]
        with self.assertRaisesRegex(experimental.CampaignError,'different bin edges'):
            campaign._combine_phenix_channels(self.m,raw)

    def test_pdf_scale_family_points_and_closure_remain_separate(self):
        raw={}
        for family,pol,scale,mpi in [('nominal',0,1.,'on'),('nominal',7,1.,'on'),
                                     ('nominal',0,.5,'on'),('mpi_off',0,1.,'off')]:
            for channel in self.m['channels']:
                raw[(family,channel,pol,0,scale,mpi)]=objects([3,1,1,3])
        for channel in self.m['channels']:
            raw[('unpolarized_closure',channel,0,0,1.,'on')]=objects([2.],helicities=['00'])
        combined,pred=campaign._combine_phenix_channels(self.m,raw)
        self.assertEqual(len(pred),4);self.assertEqual(len(combined),5)
        self.assertEqual(combined[('unpolarized_closure','photon',0,0,1.,'on')]['isolated_all']['00'].values,[4.])
        self.assertEqual(campaign._prediction_channels(self.m),['photon'])

    def test_channel_map_cannot_double_count(self):
        bad=copy.deepcopy(self.m);bad['additive_channels']['photon'].append('photon')
        with self.assertRaisesRegex(campaign.CampaignError,'partition'):
            campaign._validate_descriptor(bad,ROOT/f'config/phenomenology/{ID}.json')

    def test_job_totals_unique_seeds_and_component_cards(self):
        args=campaign.make_parser().parse_args(['full','--measurement',ID,'--tag','test',
            '--shards','100','--lo-events','3000000','--jobs','192','--dry-run'])
        plan=campaign._plan(self.m,args)
        self.assertEqual(plan['logical_jobs'],8);self.assertEqual(plan['shard_jobs'],800)
        self.assertEqual(sum(j['events'] for j in plan['jobs']),24000000)
        self.assertEqual(len({j['seed'] for j in plan['jobs']}),800)
        for channel,me in [('photon','MEGammaJet'),('fragmentation','MEQCD2to2')]:
            job=next(j for j in plan['jobs'] if j['channel']==channel)
            card=campaign._card_text(self.m,job)
            self.assertIn('MatrixElements[0] /Herwig/MatrixElements/'+me,card)
            if channel=='fragmentation': self.assertIn('JetKtCut:MinKT 4.0*GeV',card)

    def test_wrapper_external_root_and_dry_run_write_nothing(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(campaign,'CAMPAIGN_ROOT',campaign.CAMPAIGN_ROOT):
            output=io.StringIO()
            with contextlib.redirect_stdout(output):
                code=wrapper.main(['--output-root',directory,'full','--tag','dry','--dry-run'])
            plan=json.loads(output.getvalue())
            self.assertEqual(code,0)
            self.assertEqual(plan['campaign_directory'],str(Path(directory)/ID/'dry'))
            self.assertEqual(list(Path(directory).iterdir()),[])
        with contextlib.redirect_stderr(io.StringIO()),self.assertRaises(SystemExit):
            wrapper.main(['full','--measurement','OTHER','--tag','bad','--dry-run'])

    def test_baseline_reproduces_original_source_and_cards(self):
        # Pinned against the pre-expansion source and original five input cards.
        expected={'analyses/rivet/pp/MEGammaJet-Only/PHENIX_2023_I2033856.cc': 'd11053b0565c99d736a348ed52a74ae568f339610782ed863539d07ee2ebd3a5', 'cards/phenomenology/PHENIX_2023_I2033856_MEGammaJet-Only/PHENIX_2023_I2033856-Common.in': '519c1a8606cbbb9c842d2dbdda4282c008aee4c8a54e6fde92ce102f6c58b561', 'cards/phenomenology/PHENIX_2023_I2033856_MEGammaJet-Only/PHENIX_2023_I2033856_photon_MM-LO.in': 'c1455f2992f9246a4cdde3d1c3182e6a19cb7a5056d1df26abee385065cb203f', 'cards/phenomenology/PHENIX_2023_I2033856_MEGammaJet-Only/PHENIX_2023_I2033856_photon_MP-LO.in': '01a64c140ab7fd03fbc31259dcab2a3e2f6b236ec1887fa77686e7e390306642', 'cards/phenomenology/PHENIX_2023_I2033856_MEGammaJet-Only/PHENIX_2023_I2033856_photon_PM-LO.in': '1b4c65c215158ecd76e6df8ec1b4f51bd9d7af68174d20d5cc69278610018127', 'cards/phenomenology/PHENIX_2023_I2033856_MEGammaJet-Only/PHENIX_2023_I2033856_photon_PP-LO.in': '5ec239c4e8aa1a78e47888511ad794e15074fbe16fb4dc7552ec507ec0c0dc8a'}
        for relative,digest in expected.items():
            self.assertEqual(hashlib.sha256((ROOT/relative).read_bytes()).hexdigest(),digest,relative)

    def test_baseline_has_four_jobs_and_shared_validated_data(self):
        baseline=campaign.discover_pp_registry()[baseline_wrapper.MEASUREMENT]
        self.assertEqual(set(baseline['channels']),{'photon'})
        self.assertNotIn('additive_channels',baseline)
        self.assertEqual(campaign._measurement_snapshot(baseline),campaign._measurement_snapshot(self.m))
        self.assertNotEqual(campaign._signature(baseline),campaign._signature(self.m))
        args=campaign.make_parser().parse_args(['full','--measurement',baseline['id'],'--tag','baseline',
            '--shards','3','--lo-events','10','--dry-run'])
        plan=campaign._plan(baseline,args)
        self.assertEqual(plan['logical_jobs'],4);self.assertEqual(plan['shard_jobs'],12)
        self.assertEqual(sum(j['events'] for j in plan['jobs']),40)
        self.assertEqual([j['events'] for j in plan['jobs'][:3]],[4,3,3])
        for job in plan['jobs']:
            card=campaign._card_text(baseline,job)
            self.assertIn('MEGammaJet',card);self.assertNotIn('MEQCD2to2',card)
        with tempfile.TemporaryDirectory() as directory, patch.object(campaign,'CAMPAIGN_ROOT',campaign.CAMPAIGN_ROOT):
            output=io.StringIO()
            with contextlib.redirect_stdout(output):
                baseline_wrapper.main(['--output-root',directory,'full','--tag','dry','--dry-run'])
            result=json.loads(output.getvalue())
            self.assertEqual(result['logical_jobs'],4)
            self.assertEqual(Path(result['campaign_directory']),Path(directory)/baseline['id']/'dry')
            self.assertEqual(list(Path(directory).iterdir()),[])

    def test_unequal_shards_average_but_distinct_processes_add(self):
        # Already normalized estimates of a cross section: shard averaging
        # uses event fractions, and variances use their squares.
        B=experimental.BinSeries
        shards=[B([6,7],[4.],[16.]),B([6,7],[8.],[64.])]
        combined=experimental.combine_shard_series(shards,[1,3])
        self.assertEqual(combined.values,[7.])
        self.assertEqual(combined.variances,[37.])
        direct=('nominal','photon',0,0,1.,'on')
        frag=('nominal','fragmentation',0,0,1.,'on')
        _,pred=campaign._combine_phenix_channels(self.m,{
            direct:objects([7.]*4,37.),frag:objects([2.]*4,1.)})
        self.assertEqual(pred[direct]['inclusive_cross_section']['values'],[9.])
        self.assertAlmostEqual(pred[direct]['inclusive_cross_section']['errors'][0]**2,38./4)
        self.assertEqual(pred[direct]['isolated_all']['values'],[0.])

    def test_baseline_does_not_require_fragmentation(self):
        baseline=campaign.discover_pp_registry()[baseline_wrapper.MEASUREMENT]
        direct=('nominal','photon',0,0,1.,'on')
        _,pred=campaign._combine_phenix_channels(baseline,{direct:objects([12.,8.,8.,12.])})
        self.assertEqual(pred[direct]['isolated_all']['values'],[.2])
        self.assertFalse(any(name.startswith('Component_') for name in pred[direct]))

    def sliced_plan(self, edges='4,8,12,16,20', measurement=None):
        measurement = measurement or self.m
        args = campaign.make_parser().parse_args([
            'full', '--measurement', measurement['id'], '--tag', 'sliced',
            '--lo-events', '10', '--shards', '3', '--seed-base', '80000000',
            f'--hard-pt-edges={edges}', '--dry-run'])
        return campaign._plan(measurement, args)

    def test_slices_have_separate_cards_seeds_and_complete_helicities(self):
        for measurement, components in [(self.m, 2), (
                campaign.discover_pp_registry()[baseline_wrapper.MEASUREMENT], 1)]:
            plan = self.sliced_plan(measurement=measurement)
            count = 4 * 4 * components * 3
            self.assertEqual(plan['logical_jobs'], 4 * 4 * components)
            self.assertEqual(plan['shard_jobs'], count)
            self.assertEqual(sum(j['events'] for j in plan['jobs']), 4 * 4 * components * 10)
            for field in ('id', 'seed', 'output_yoda'):
                self.assertEqual(len({j[field] for j in plan['jobs']}), count)
            self.assertEqual(len({j['stem'] for j in plan['jobs']}), 4 * 4 * components)
            self.assertEqual([j['seed'] for j in plan['jobs']],
                             list(range(80000000, 80000000 + count)))
            for job in plan['jobs']:
                low, high = job['hard_pt_range_gev']
                cut = 'PhotonKtCut' if job['channel'] == 'photon' else 'JetKtCut'
                card = campaign._card_text(measurement, job)
                self.assertIn(f'set /Herwig/Cuts/{cut}:MinKT {low:g}*GeV', card)
                self.assertIn(f'set /Herwig/Cuts/{cut}:MaxKT {high:g}*GeV', card)
                self.assertLess(card.index(f'{cut}:MaxKT'), card.index('saverun '))
            # Immutable configuration survives serialization, including the tail sentinel.
            self.assertEqual(json.loads(json.dumps(plan)), plan)

    def test_open_tail_and_invalid_edges(self):
        plan = self.sliced_plan('4,8,12,16,20,inf')
        self.assertEqual(plan['options']['hard_pt_ranges_gev'][-1], [20., None])
        for job in plan['jobs']:
            if job['hard_pt_range_gev'][1] is None:
                card = campaign._card_text(self.m, job)
                self.assertIn(':MinKT 20*GeV', card)
                self.assertNotIn(':MaxKT', card)
        for edges in ('4', '', '4,4,7', '4,3', '0,4', '-1,4', '4,nan',
                      'nan,7', '4,inf,16', '4,-inf', 'a,7', '4,'):
            with self.subTest(edges=edges), self.assertRaises(campaign.CampaignError):
                self.sliced_plan(edges)
        star = campaign.discover_pp_registry()['STAR_2019_I1708793']
        with self.assertRaisesRegex(campaign.CampaignError, 'only valid for PHENIX'):
            self.sliced_plan(measurement=star)

    def test_slices_add_after_shard_averaging_before_asymmetry(self):
        jobs = []
        fixture = {}
        # First slice: unequal shards -> same-helicity sigma=7, opposite=3.
        # Second slice: same-helicity sigma=1, opposite=9.
        # The slice asymmetries are +0.4 and -0.8; the sum is -0.2.
        for h in ('PP', 'PM', 'MP', 'MM'):
            same = h in ('PP', 'MM')
            for n, bounds, events, value, variance in [
                (0, [4.,8.], 1, 4. if same else 3., 16.),
                (1, [4.,8.], 3, 8. if same else 3., 64.),
                (2, [8.,12.], 20, 1. if same else 9., 1.),
            ]:
                path = f'{h}-{n}.yoda'
                jobs.append({'helicity': h, 'events': events, 'output_yoda': path,
                             'hard_pt_range_gev': bounds})
                fixture[path] = experimental.BinSeries([6.,7.], [value], [variance])
        def read(path, requested):
            return {name: fixture[path.name] for name in requested}
        requested = {name: ('PHENIX_2023_I2033856', name) for name in (
            'inclusive_cross_section', 'isolated_cross_section', 'isolated_all')}
        with patch.object(experimental, 'read_histogram_series_many', side_effect=read) as reader:
            loaded = {h: campaign._load_series_many(
                [j for j in jobs if j['helicity'] == h], Path('/unused'), requested)
                for h in ('PP', 'PM', 'MP', 'MM')}
        self.assertEqual(reader.call_count, 12)  # one read per file, not per observable
        raw = {name: {h: loaded[h][name] for h in loaded} for name in requested}
        self.assertEqual(raw['isolated_all']['PP'].values, [8.])
        self.assertEqual(raw['isolated_all']['PM'].values, [12.])
        self.assertEqual(raw['isolated_all']['PP'].variances, [38.])
        prediction = campaign._phenix_prediction(raw)
        self.assertAlmostEqual(prediction['isolated_all']['values'][0], -.2)
        self.assertAlmostEqual(prediction['inclusive_cross_section']['values'][0], 10.)
        gradient = [(s * 40. + 8.) / 40.**2 for s in (1.,-1.,-1.,1.)]
        self.assertAlmostEqual(prediction['isolated_all']['errors'][0]**2,
                               38. * sum(g*g for g in gradient))

    def test_overlapping_or_mixed_slices_are_rejected_before_loading(self):
        for bounds in [([4.,8.],[7.,10.]), ([4.,None],[7.,10.]), ([4.,7.],None)]:
            jobs = [{'hard_pt_range_gev': b} if b else {} for b in bounds]
            with self.subTest(bounds=bounds), self.assertRaises(campaign.CampaignError):
                campaign._load_series_many(jobs, Path('/unused'), {'yield': ('A','Y')})

    def test_missing_slice_is_rejected(self):
        plan = self.sliced_plan('4,8,12')
        manifest = {'configuration': plan['options'], 'jobs': plan['jobs']}
        for job in manifest['jobs']:
            job['status'] = 'success'
        manifest['jobs'] = [j for j in manifest['jobs'] if not (
            j['channel'] == 'photon' and j['helicity'] == 'PP'
            and j['hard_pt_range_gev'] == [8.,12.])]
        with patch.object(experimental, '_reconcile_manifest_outputs'), \
                patch.object(experimental, '_nonempty', return_value=True), \
                self.assertRaisesRegex(campaign.CampaignError, 'Incomplete hard-pT'):
            campaign._logical_groups(manifest, Path('/unused'))

    def test_wrapper_sliced_dry_run_writes_nothing(self):
        with tempfile.TemporaryDirectory() as directory, \
                patch.object(campaign, 'CAMPAIGN_ROOT', campaign.CAMPAIGN_ROOT):
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                code = wrapper.main(['--output-root', directory, 'full', '--tag', 'slice-dry',
                    '--hard-pt-edges', '4,8,12,16,20', '--lo-events', '3000000',
                    '--shards', '100', '--jobs', '192', '--dry-run'])
            plan = json.loads(output.getvalue())
            self.assertEqual(code, 0)
            self.assertEqual(plan['shard_jobs'], 3200)
            self.assertEqual(sum(j['events'] for j in plan['jobs']), 96000000)
            self.assertEqual(list(Path(directory).iterdir()), [])

    def budget_plan(self, *extra, measurement=None):
        measurement = measurement or self.m
        args = campaign.make_parser().parse_args([
            'full', '--measurement', measurement['id'], '--tag', 'budgets',
            '--lo-events', '10', '--shards', '3', '--seed-base', '81000000',
            '--hard-pt-edges', '4,8,inf', '--dry-run', *extra])
        return campaign._plan(measurement, args), args

    def test_component_slice_budgets_and_remainders(self):
        plan, args = self.budget_plan('--hard-pt-events', 'fragmentation=31,11')
        grouped = defaultdict(list)
        for job in plan['jobs']:
            grouped[(job['channel'], job['helicity'], job['hard_pt_range_gev'][0])].append(job)
        for (channel, helicity, low), jobs in grouped.items():
            expected = ([11,10,10] if low == 4. else [4,4,3]) if channel == 'fragmentation' else [4,3,3]
            self.assertEqual([job['events'] for job in jobs], expected)
        self.assertEqual(sum(j['events'] for j in plan['jobs']), 248)
        self.assertEqual(plan['shard_jobs'], 48)
        self.assertEqual([j['seed'] for j in plan['jobs']], list(range(81000000,81000048)))
        for field in ('id', 'output_yoda'):
            self.assertEqual(len({j[field] for j in plan['jobs']}), 48)
        self.assertEqual(plan['options']['hard_pt_events_by_channel'], {'fragmentation':[31,11]})
        self.assertEqual(json.loads(json.dumps(plan)), plan)
        config = campaign._manifest_configuration(self.m, args, plan)
        repeated, repeated_args = self.budget_plan('--hard-pt-events', 'fragmentation=31,11')
        self.assertEqual(config, campaign._manifest_configuration(self.m, repeated_args, repeated))
        changed, changed_args = self.budget_plan('--hard-pt-events', 'fragmentation=32,11')
        self.assertNotEqual(config, campaign._manifest_configuration(self.m, changed_args, changed))
        self.assertEqual(config['lo_events'], 10)  # Default for unspecified components.

    def test_both_components_and_original_variant_accept_overrides(self):
        plan, _ = self.budget_plan('--hard-pt-events', 'fragmentation=31,11',
                                   '--hard-pt-events', 'photon=13,14')
        self.assertEqual(sum(j['events'] for j in plan['jobs']), 4*(31+11+13+14))
        reverse, _ = self.budget_plan('--hard-pt-events', 'photon=13,14',
                                      '--hard-pt-events', 'fragmentation=31,11')
        self.assertEqual(plan, reverse)
        baseline = campaign.discover_pp_registry()[baseline_wrapper.MEASUREMENT]
        original, _ = self.budget_plan('--hard-pt-events', 'photon=13,14', measurement=baseline)
        self.assertEqual(sum(j['events'] for j in original['jobs']), 4*(13+14))
        self.assertEqual(original['shard_jobs'], 24)
        with self.assertRaisesRegex(campaign.CampaignError, 'available channels: photon'):
            self.budget_plan('--hard-pt-events', 'fragmentation=13,14', measurement=baseline)

    def test_invalid_budgets_are_rejected_before_runtime(self):
        for spec in ('31,11', 'unknown=31,11', 'fragmentation=31',
                     'fragmentation=31,11,7', 'fragmentation=0,11',
                     'fragmentation=-1,11', 'fragmentation=1.5,11',
                     'fragmentation=inf,11', 'fragmentation=31,', 'fragmentation='):
            with self.subTest(spec=spec), self.assertRaises(campaign.CampaignError):
                self.budget_plan('--hard-pt-events', spec)
        with self.assertRaisesRegex(campaign.CampaignError, 'Repeated'):
            self.budget_plan('--hard-pt-events', 'photon=31,11', '--hard-pt-events', 'photon=31,11')
        with self.assertRaisesRegex(campaign.CampaignError, 'smoke'):
            self.budget_plan('--hard-pt-events', 'photon=31,11', '--smoke')
        args = campaign.make_parser().parse_args([
            'full','--measurement',ID,'--tag','bad','--hard-pt-events','photon=31,11'])
        with self.assertRaisesRegex(campaign.CampaignError, 'requires --hard-pt-edges'):
            campaign._plan(self.m, args)
        for identifier in ('STAR_2019_I1708793','COMPASS_2010_I843494'):
            with contextlib.redirect_stderr(io.StringIO()) as error, \
                    patch.object(experimental, 'main') as legacy:
                code = campaign.main(['full','--measurement',identifier,'--tag','bad',
                                      '--hard-pt-events','photon=31,11','--dry-run'])
            self.assertEqual(code, 2)
            self.assertIn('only valid for PHENIX', error.getvalue())
            legacy.assert_not_called()

    def test_production_variable_budget_dry_run_writes_nothing(self):
        with tempfile.TemporaryDirectory() as directory, \
                patch.object(campaign, 'CAMPAIGN_ROOT', campaign.CAMPAIGN_ROOT), \
                patch.object(campaign, '_runtime', side_effect=AssertionError('No runtime in dry-run')):
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                code = wrapper.main(['--output-root',directory,'full','--tag','lowpt',
                    '--lo-events','600000','--shards','100','--jobs','192',
                    '--hard-pt-edges','4,8,12,16,20,inf',
                    '--hard-pt-events','fragmentation=60000000,6000000,600000,600000,600000',
                    '--dry-run'])
            plan = json.loads(output.getvalue())
            self.assertEqual(code, 0)
            self.assertEqual(plan['shard_jobs'], 4000)
            self.assertEqual(sum(j['events'] for j in plan['jobs']), 283200000)
            for channel, total in [('photon',12000000),('fragmentation',271200000)]:
                self.assertEqual(sum(j['events'] for j in plan['jobs'] if j['channel']==channel), total)
            self.assertEqual(list(Path(directory).iterdir()), [])

    def test_default_budget_is_preserved_without_overrides(self):
        plan, args = self.budget_plan()
        self.assertNotIn('hard_pt_events_by_channel', plan['options'])
        self.assertNotIn('lo_events_by_channel', plan['options'])
        self.assertEqual(sum(j['events'] for j in plan['jobs']), 160)

    def test_unequal_component_slice_budgets_preserve_cross_sections(self):
        plan, _ = self.budget_plan('--hard-pt-events', 'fragmentation=31,11')
        fixture = {}
        for job in plan['jobs']:
            same = job['helicity'] in ('PP','MM')
            low = job['hard_pt_range_gev'][0]
            # Already normalized shard estimates of fixed physical cross sections.
            rates = {'photon': {4.: (12.,8.), 8.: (3.,1.)},
                     'fragmentation': {4.: (1.,9.), 8.: (2.,4.)}}
            value = rates[job['channel']][low][0 if same else 1]
            fixture[job['output_yoda']] = experimental.BinSeries(
                [6.,7.], [value], [value**2/job['events']])
        def read(path, requested):
            return {name: fixture[str(path)] for name in requested}
        requested = {name: (ID, name) for name in (
            'inclusive_cross_section','isolated_cross_section','isolated_all')}
        samples = {}
        with patch.object(experimental, 'read_histogram_series_many', side_effect=read):
            for channel in ('photon','fragmentation'):
                loaded = {h: campaign._load_series_many(
                    [j for j in plan['jobs'] if j['channel']==channel and j['helicity']==h],
                    Path('.'), requested) for h in ('PP','PM','MP','MM')}
                samples[('nominal',channel,0,0,1.,'on')] = {
                    name: {h: loaded[h][name] for h in loaded} for name in requested}
        _, predictions = campaign._combine_phenix_channels(self.m, samples)
        prediction = predictions[('nominal','photon',0,0,1.,'on')]
        self.assertAlmostEqual(prediction['inclusive_cross_section']['values'][0], 20.)
        self.assertAlmostEqual(prediction['isolated_all']['values'][0], -.1)
        same_var = 12.**2/10 + 3.**2/10 + 1.**2/31 + 2.**2/11
        opposite_var = 8.**2/10 + 1.**2/10 + 9.**2/31 + 4.**2/11
        expected_variance = 2*same_var*(1.1/80)**2 + 2*opposite_var*(-.9/80)**2
        self.assertAlmostEqual(prediction['isolated_all']['errors'][0]**2, expected_variance)

if __name__=='__main__':unittest.main()
