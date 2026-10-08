"""Additive, provenance-checked HERMES ratio galleries; existing plots are retained."""
from __future__ import annotations
import argparse
import csv
import hashlib
import html
import json
import os
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
MEASUREMENT = 'HERMES_2007_I726689'
REFERENCE = ROOT/'data/experimental'/MEASUREMENT/'reference.json'
BORN_REFERENCE = REFERENCE.with_name('born-apar-reference.json')
MANIFEST = 'cache-manifest.json'
SETTINGS = {'version': 1, 'ratio': 'MC/data', 'focus_limits': [-1.,3.],
            'mc_uncertainty': 'MC statistical / absolute data central value',
            'data_uncertainty': 'Separate relative statistical and outer errors about unity',
            'zero_compatible': 'Keep nonzero signed ratios; use hollow markers',
            'zero_data': 'Undefined; no epsilon substitution',
            'overflow': 'Explicit triangles; central values labelled; full-range companions'}
SOURCES = ['hermes_data_ratio_gallery.py','hermes_ratio_common.py','hermes_born_ratio_plots.py',
           'hermes_projection_ratio_plots.py','hermes_born_cell_plots.py']
BEGIN = '<!-- BEGIN HERMES DATA RATIOS -->'
END = '<!-- END HERMES DATA RATIOS -->'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def input_hashes(summary, integrated):
    paths={'summary.json':summary,'integrated.json':integrated,
           'reference.json':REFERENCE,'born-apar-reference.json':BORN_REFERENCE}
    paths.update({f'scripts/{name}':ROOT/'scripts'/name for name in SOURCES})
    return {name:digest(path) for name,path in paths.items()}


def expected_outputs():
    from hermes_born_ratio_plots import REQUIRED_OUTPUTS as born
    from hermes_projection_ratio_plots import REQUIRED_OUTPUTS as projections
    return born | projections | {'index.html', 'ratios.json', 'ratios.csv'}


def _valid_cache(directory, key):
    try:
        if directory.is_symlink() or not directory.is_dir():return False
        record=json.loads((directory/MANIFEST).read_text())
        expected=record['outputs']
        if record['cache_key']!=key or not isinstance(expected,dict) or set(expected)!=expected_outputs():return False
        if 'index.html' not in expected or 'ratios.json' not in expected or 'ratios.csv' not in expected:return False
        if hashlib.sha256(json.dumps({'inputs':record['inputs'],'settings':record['settings']},sort_keys=True).encode()).hexdigest()!=key:return False
        for name,checksum in expected.items():
            if Path(name).name!=name:return False
            path=directory/name
            if path.is_symlink() or not path.is_file() or not path.stat().st_size or digest(path)!=checksum:return False
        return set(p.name for p in directory.iterdir())==set(expected)|{MANIFEST}
    except (OSError,ValueError,KeyError,TypeError):return False


def _info(campaign, directory, key, created):
    return {'directory':Path(os.path.relpath(directory,campaign)).as_posix(),
            'index':Path(os.path.relpath(directory/'index.html',campaign)).as_posix(),
            'cache_key':key,'created':created,'ratio':'MC/data','theory_family':'nominal'}


def _write_index(output, comparisons):
    def figure(stem,title,full=None):
        extra=f' · <a href="{full}.pdf">Full-range PDF</a>' if full else ''
        return (f'<article><h3>{html.escape(title)}</h3><p><a href="{stem}.pdf">PDF</a> · '
                f'<a href="{stem}.png">PNG</a>{extra}</p><a href="{stem}.pdf">'
                f'<img src="{stem}.png" alt="{html.escape(title)}"></a></article>')
    blocks=[]
    for target,label in (('P','Proton'),('D','Deuteron')):
        stem=f'AParallel_{target}_BornCells_Ratio_5x4'
        blocks.append(figure(stem,label+': 45 Born cells, MC/data',stem+'_FullRange'))
    blocks.append('<h2>Absolute comparisons with ratio panels</h2><div class="pairs">')
    for comp in comparisons:
        title = comp['title']
        if comp.get('rivet_id'):
            target = 'proton' if comp['id'] == 'A1_P_Q2GT1' else 'deuteron'
            title = f"{comp['rivet_id']} — HERMES {target} A1 (Q² > 1 GeV²)"
        blocks.append(figure(comp['id']+'_ratio',title,comp['id']+'_ratio_fullrange'))
    blocks.append('</div><h2>Individual fixed-x Born comparisons</h2>')
    for target,label in (('P','Proton'),('D','Deuteron')):
        blocks.append(f'<details><summary>{label}: all 19 fixed-x slices</summary><div class="pairs">')
        for index in range(1,20):
            stem=f'AParallel_{target}_X{index:02d}_WithRatio'
            blocks.append(figure(stem,f'{label}, x slice {index}'))
        blocks.append('</div></details>')
    page='''<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1"><title>HERMES MC/data ratios</title>
<style>body{font:16px system-ui,sans-serif;max-width:1150px;margin:2rem auto;padding:0 1rem;color:#222}p{line-height:1.5}img{width:100%;height:auto}article{margin:1.5rem 0}a{color:#075b9c}.pairs{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:1.4rem}summary{cursor:pointer;padding:1rem;background:#f3f3f3}@media(max-width:700px){.pairs{grid-template-columns:1fr}}</style>
<h1>HERMES: additional MC/data comparisons</h1>
<p>Original absolute plots are retained. Orange ratios divide the nominal Herwig prediction by the experimental central value. Orange bands show MC statistical error divided by the absolute data central value. Dark/light gray bands about unity show relative experimental statistical/outer errors separately. These are display uncertainties, not Gaussian confidence intervals for an uncertain quotient.</p>
<p>Open orange markers mean the experimental outer interval contains zero. Signed nonzero denominators are retained; exact-zero denominators are undefined. Triangles mark truncated bands or ratios; truncated central values are labelled. Full-range companions show every plotted uncertainty without clipping. Their panel scales may differ.</p>
<p>Born-cell ratios use direct LL/UU. A1 uses the existing inverse-D proxy and published A1 values. Reconstructed projections compare the same GD11-weighted cells on both sides; their outer data errors retain the conservative systematic bound, and MC errors retain the independent-cell approximation. No new fit, depolarization, target correction, or uncertainty is added. The original experimental covariance remains in the reconstructed-data snapshot; no independent-bin goodness-of-fit is inferred here.</p>
<p><a href="ratios.json">Numerical JSON and provenance</a> · <a href="ratios.csv">Numerical CSV</a></p>
'''+''.join(blocks)+'</html>\n'
    (output/'index.html').write_text(page,encoding='utf-8')


def ensure_campaign_plots(campaign_dir, output_dir, integrated_path=None):
    import hermes_born_ratio_plots as born
    import hermes_projection_ratio_plots as projections
    campaign,output=Path(campaign_dir).resolve(),Path(output_dir).resolve()
    summary_path=campaign/'postprocess/summary.json'
    if integrated_path is None:
        record=json.loads((campaign/'manifest.json').read_text())
        integrated_path=campaign/record['plots']['reconstructed_born_projections']['directory']/'integrated.json'
    integrated_path=Path(integrated_path).resolve()
    inputs=input_hashes(summary_path,integrated_path)
    summary=json.loads(summary_path.read_text())
    integrated=json.loads(integrated_path.read_text())
    if summary.get('measurement')!=MEASUREMENT or integrated.get('measurement')!=MEASUREMENT:
        raise ValueError('Ratio plots require matching HERMES inputs')
    if summary.get('tag')!=integrated.get('campaign_tag'):
        raise ValueError('Ratio inputs belong to different campaigns')
    if integrated['input_files']['summary']['sha256']!=inputs['summary.json']:
        raise ValueError('Reconstructed ratios use a stale theory summary')
    if integrated['input_files']['reference']['sha256']!=inputs['born-apar-reference.json']:
        raise ValueError('Reconstructed ratios use a different Born reference')
    provenance={'inputs':inputs,'settings':SETTINGS}
    key=hashlib.sha256(json.dumps(provenance,sort_keys=True).encode()).hexdigest()
    cache=output/MEASUREMENT/'data-ratios'
    for directory in sorted(cache.glob(key+'*')):
        if _valid_cache(directory,key):
            if inputs!=input_hashes(summary_path,integrated_path):raise ValueError('Ratio inputs changed; retry')
            return _info(campaign,directory,key,False)
    direct=born.build_snapshot(summary,json.loads(BORN_REFERENCE.read_text()))
    comparisons=projections.build_comparisons(summary,json.loads(REFERENCE.read_text()),integrated)
    cache.mkdir(parents=True,exist_ok=True)
    attempt=0
    while True:
        directory=cache/(key if not attempt else f'{key}-{attempt:03d}')
        try:directory.mkdir();break
        except FileExistsError:attempt+=1
    paths=born.render(direct,directory)+projections.render(comparisons,directory)
    _write_index(directory,comparisons)
    snapshot={'measurement':MEASUREMENT,'tag':summary.get('tag'),'settings':SETTINGS,
              'provenance':provenance,'born_cells':direct,'comparisons':comparisons}
    (directory/'ratios.json').write_text(json.dumps(snapshot,indent=2,allow_nan=False)+'\n')
    columns=['comparison','bin','low','high','marker','data','stat','total','mc','mc_stat','supported',
             'ratio','ratio_mc_stat','data_relative_stat','data_relative_total',
             'denominator_zero_compatible','ratio_missing_reason']
    with (directory/'ratios.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=columns,extrasaction='ignore');writer.writeheader()
        for target,block in direct['targets'].items():
            for panel in block['panels']:
                for row in panel['ratio_rows']:
                    writer.writerow({'comparison':panel['selection'],**row})
        for comparison in comparisons:
            for row in comparison['rows']:writer.writerow({'comparison':comparison['id'],**row})
    paths += [directory/name for name in ('index.html','ratios.json','ratios.csv')]
    outputs={p.name:digest(p) for p in paths}
    if len(outputs)!=len(paths) or set(outputs)!=expected_outputs() or any(not p.stat().st_size for p in paths):raise ValueError('Invalid ratio output inventory')
    if inputs!=input_hashes(summary_path,integrated_path):raise ValueError('Ratio inputs changed while rendering; retry')
    record={'schema_version':1,'cache_key':key,**provenance,'outputs':outputs}
    (directory/MANIFEST).write_text(json.dumps(record,indent=2)+'\n')
    if not _valid_cache(directory,key):raise ValueError('Incomplete ratio gallery')
    return _info(campaign,directory,key,True)


def gallery_section(campaign, output, info, page_directory):
    campaign,output,page_directory=map(lambda p:Path(p).resolve(),(campaign,output,page_directory))
    gallery=(campaign/info['directory']).resolve()
    gallery.relative_to(output/MEASUREMENT/'data-ratios')
    if not _valid_cache(gallery,info['cache_key']):raise ValueError('Ratio gallery failed checksum validation')
    def link(name):return html.escape(Path(os.path.relpath(gallery/name,page_directory)).as_posix())
    items=['<h3 id="a1-ratio-comparisons">Published A<sub>1</sub> comparisons with ratio panels</h3>',
           '<div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,360px),1fr));gap:1.5rem">']
    for target,label,rivet_id in (('P','Proton','d14-x01-y01'),('D','Deuteron','d14-x01-y02')):
        stem=f'A1_{target}_Q2GT1_ratio'
        original=output/MEASUREMENT/f'{rivet_id}.pdf'
        original_link=(' · <a href="'+html.escape(Path(os.path.relpath(original,page_directory)).as_posix())
                       +'">Original absolute plot</a>') if original.is_file() else ''
        items.append(f'<article id="{rivet_id}-ratio"><h4>{rivet_id}: {label} A<sub>1</sub></h4>'
                     f'<a href="{link(stem+".pdf")}"><img src="{link(stem+".png")}" alt="{rivet_id}: {label} A1 with MC/data panel"></a>'
                     f'<p><a href="{link(stem+".pdf")}">Paper PDF</a> · '
                     f'<a href="{link(stem+".png")}">PNG</a> · '
                     f'<a href="{link(stem+"_fullrange.pdf")}">Full range</a>{original_link}</p></article>')
    items.append('</div><h3>Born-cell ratio panels</h3>')
    for target,label in (('P','Proton'),('D','Deuteron')):
        stem=f'AParallel_{target}_BornCells_Ratio_5x4'
        items.append(f'<article><h3>{label}: MC/data in the 45 Born cells</h3>'
                     f'<a href="{link(stem+".pdf")}"><img src="{link(stem+".png")}" alt="{label} Born MC/data ratios"></a>'
                     f'<p><a href="{link(stem+".pdf")}">Paper PDF</a> · '
                     f'<a href="{link(stem+"_FullRange.pdf")}">Full range</a></p></article>')
    return (BEGIN+'<section class="plot"><h2>Additional MC/data ratio figures</h2>'
            '<p>Separate absolute-plus-ratio figures and 5×4 ratio panels. Orange bands are MC statistical; '
            'gray reference bands show relative experimental uncertainties. Open markers flag data compatible with zero. '
            'Triangles mark truncation; labelled central ratios and full-range companions retain large values.</p>'
            f'<p><a href="{link("index.html")}">All ratio figures: A1, Born cells and reconstructed projections</a> · '
            f'<a href="{link("ratios.csv")}">CSV</a></p>'+''.join(items)+'</section>'+END)


def append_to_existing_indexes(campaign, output, info):
    """Change only the bounded additive HTML section; never rerender old plots."""
    campaign,output=Path(campaign).resolve(),Path(output).resolve()
    for page in (output/'index.html',output/MEASUREMENT/'index.html'):
        original=page.read_text(encoding='utf-8')
        section=gallery_section(campaign,output,info,page.parent)
        if BEGIN in original:
            if original.count(BEGIN)!=1 or original.count(END)!=1:raise ValueError('Ambiguous ratio HTML markers')
            updated=re.sub(re.escape(BEGIN)+'.*?'+re.escape(END),lambda m:section,original,flags=re.S)
        else:
            if '</body>' not in original:raise ValueError('Cannot locate existing gallery body')
            updated=original.replace('</body>',section+'\n</body>',1)
        if updated!=original:
            backup=campaign/'work/ratio-gallery-index-backups'/digest(page)
            backup.mkdir(parents=True,exist_ok=True)
            (backup/('root.html' if page.parent==output else 'analysis.html')).write_text(original,encoding='utf-8')
            temp=page.with_name('.index-ratios.tmp');temp.write_text(updated,encoding='utf-8');temp.replace(page)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign',required=True,type=Path)
    parser.add_argument('--integrated',type=Path)
    parser.add_argument('--publish',action='store_true',help='Append to existing indexes and record the additive gallery')
    args=parser.parse_args()
    campaign=args.campaign.resolve();output=campaign/'plots'
    info=ensure_campaign_plots(campaign,output,args.integrated)
    if args.publish:
        from run_experimental_campaign import atomic_write_json,utc_now
        append_to_existing_indexes(campaign,output,info)
        path=campaign/'manifest.json';record=json.loads(path.read_text())
        record.setdefault('plots',{})['data_ratios']=info
        record.setdefault('history',[]).append({'at':utc_now(),'action':'add-data-ratio-plots'})
        atomic_write_json(path,record)
    print(json.dumps(info,indent=2))


if __name__=='__main__':main()
