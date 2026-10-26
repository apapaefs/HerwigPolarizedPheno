#!/usr/bin/env python3
"""Compile and exercise PHENIX with controlled HepMC events (not production MC)."""
import argparse
import json
import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess

ROOT=Path(__file__).resolve().parents[1]
ID='PHENIX_2023_I2033856'
EDGES=[6,6.5,7,7.5,8,8.5,9,9.5,10,12,14,16,18,20,22,24,26,28,30]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--variant',choices=['expanded','MEGammaJet-Only'],default='expanded')
    args=parser.parse_args()
    output=args.output.resolve();output.mkdir(parents=True,exist_ok=True)
    legacy=args.variant=='MEGammaJet-Only'
    source_dir=ROOT/'analyses/rivet/pp'/('MEGammaJet-Only' if legacy else '')
    env=dict(os.environ)
    env['RIVET_ANALYSIS_PATH']=str(output)+os.pathsep+env.get('RIVET_ANALYSIS_PATH','')
    env['RIVET_DATA_PATH']=str(source_dir)+os.pathsep+env.get('RIVET_DATA_PATH','')
    def run(command,name):
        with (output/(name+'.log')).open('w') as log:
            subprocess.run(command,cwd=output,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    flags=lambda option:shlex.split(subprocess.check_output(['rivet-config',option],text=True))
    cxx=shlex.split(os.environ.get('RIVET_CXX',shutil.which('g++') or 'c++'))
    run(cxx+['-std=c++17','-shared','-fPIC']+flags('--cppflags')+
        [str(source_dir/f'{ID}.cc')]+flags('--libs')+
        ['-o',str(output/f'Rivet{ID}.so')],'build-plugin')
    run(cxx+['-std=c++17']+flags('--cppflags')+
        [str(ROOT/'scripts/tests/phenix2023_fidelity_events.cc')]+flags('--libs')+
        ['-o',str(output/'fixture-events')],'build-fixtures')
    run([str(output/'fixture-events'),str(output)],'create-fixtures')
    run(['rivet','--show-analysis',ID],'analysis-info')
    info=(output/'analysis-info.log').read_text()
    assert 'Parse error' not in info and 'UNVALIDATED' in info and '255' in info,info
    import yoda
    report={}
    rejected={'below-pt','above-pt','outside-eta','hadron-decay','zero-photons'}
    nonisolated={'cone-fail','charged-momentum','neutral-above','energy-not-et'}
    if legacy: nonisolated.remove('charged-momentum')
    for file in sorted(output.glob('*.hepmc')):
        name=file.stem
        run(['rivet','-a',ID,'-o',str(output/(name+'.yoda')),str(file)],name)
        histos=yoda.read(str(output/(name+'.yoda')))
        def values(label):
            return [float(b.val())*(float(b.xMax())-float(b.xMin())) for b in histos[f'/{ID}/{label}'].bins()]
        inc=values('InclusivePhotonCrossSection');iso=values('IsolatedPhotonCrossSection');yy=values('IsolatedPhotonYield')
        def variances(label):
            out=[]
            for b in histos[f'/{ID}/{label}'].bins():
                try: error=float(b.errAvg('stats'))
                except Exception:
                    assert float(b.val())==0., (name,label,b)
                    error=0.
                out.append(error**2*(float(b.xMax())-float(b.xMin()))**2)
            return out
        if name in {'two-same-bin','weighted-multiplicity'}:
            weighted=name=='weighted-multiplicity'
            factor=2. if weighted else 1.  # crossSection/sum(all event weights)
            a=2. if weighted else 1.; b=3. if weighted else 0.
            photon_sum=2*a+b
            photon_sum2=(2*a*a if legacy else (2*a)**2)+b*b
            inv_sum=a*(1/7.1+1/7.4)+b/7.2
            inv_sum2=(a*a*(1/7.1**2+1/7.4**2) if legacy
                      else (a*(1/7.1+1/7.4))**2)+(b/7.2)**2
            for label,vals in [('InclusivePhotonCrossSection',inc),('IsolatedPhotonCrossSection',iso)]:
                assert math.isclose(sum(vals),factor*inv_sum/math.pi,rel_tol=2e-6),(name,vals)
                vv=variances(label)
                assert math.isclose(sum(vv),(factor/math.pi)**2*inv_sum2,rel_tol=3e-6),(name,vv,inv_sum2)
            assert math.isclose(sum(yy),factor*photon_sum,rel_tol=2e-6),(name,yy)
            yyvar=variances('IsolatedPhotonYield')
            assert math.isclose(sum(yyvar),factor**2*photon_sum2,rel_tol=3e-6),(name,yyvar,photon_sum2)
        elif name.startswith('bin-'):
            i=int(name[4:]);pt=(EDGES[i]+EDGES[i+1])/2
            expected=1/(math.pi*pt)
            assert math.isclose(inc[i],expected,rel_tol=1e-6),(name,inc,expected)
            assert math.isclose(sum(inc),expected,rel_tol=1e-6)
            assert inc==iso
            assert math.isclose(sum(yy),1. if pt<20 else 0.,abs_tol=1e-6),(name,yy)
        else:
            pt=8. if name=='quark-radiation' else 6.2
            expected=0. if name in rejected else (1/7.+1/8.)/math.pi if name=='two-photons' else 1/(math.pi*pt)
            assert math.isclose(sum(inc),expected,rel_tol=1e-6,abs_tol=1e-8),(name,inc,expected)
            expected_iso=0. if name in nonisolated else expected
            assert math.isclose(sum(iso),expected_iso,rel_tol=1e-6,abs_tol=1e-8),(name,iso,expected_iso)
            assert math.isclose(sum(yy),0. if name in rejected|nonisolated else 2. if name=='two-photons' else 1.,abs_tol=1e-6),(name,yy)
        report[name]={'inclusive_integral_pb_per_GeV':sum(inc),'isolated_integral_pb_per_GeV':sum(iso),'isolated_yield_pb':sum(yy),'sum_diagonal_yield_variances_pb2':sum(variances('IsolatedPhotonYield'))}
    (output/'fixture-report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'Passed {len(report)} controlled events: binning, prompt ancestry, isolation and normalization, and event-level variances ('+args.variant+').')

if __name__=='__main__':main()
