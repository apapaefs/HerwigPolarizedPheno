"""Paired aperture errors checked against independent event-count derivatives."""
import math
from pathlib import Path
import sys
import unittest
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
import run_experimental_campaign as c

class AcceptanceCovarianceTests(unittest.TestCase):
    def calculate(self, outside):
        labels=("PP","PM","MP","MM")
        n=[70.,20.,30.,80.]
        kp=[1.3,1.5,1.7,1.9]
        ko=[1.1,1.2,1.4,1.6]
        uu={s:.25 for s in labels};ll=dict(zip(labels,[.25,-.25,-.25,.25]))
        def result(extra):
            ordinary={s:c.BinSeries([0.,1.],[n[i]+extra[i]],[n[i]+extra[i]]) for i,s in enumerate(labels)}
            weighted={s:c.BinSeries([0.,1.],[kp[i]*n[i]+ko[i]*extra[i]],[kp[i]**2*n[i]+ko[i]**2*extra[i]]) for i,s in enumerate(labels)}
            cov={s:c.BinSeries([0.,1.],[0.],[kp[i]*n[i]+ko[i]*extra[i]]) for i,s in enumerate(labels)}
            return c._asymmetry_outputs(ordinary,weighted,cov,uu,ll)
        p=result([0.]*4);r=result(outside)
        c._paired_acceptance_difference(p,r)
        for obs,k in (("a_parallel",[1.]*4),("a1",kp)):
            nk="sigma_ll" if obs=="a_parallel" else "sigma_ll_over_d"
            up=p["sigma_uu"].values[0];ur=r["sigma_uu"].values[0]
            lp=p[nk].values[0];lr=r[nk].values[0]
            variance=0.
            for i,s in enumerate(labels):
                dp=ll[s]*k[i]/up-lp*uu[s]/up**2
                dr=ll[s]*k[i]/ur-lr*uu[s]/ur**2
                outer_k=1. if obs=="a_parallel" else ko[i]
                do=ll[s]*outer_k/ur-lr*uu[s]/ur**2
                variance += n[i]*(dr-dp)**2+outside[i]*do**2
            self.assertAlmostEqual(r["paired_difference"][obs]["errors"][0]**2,variance,places=14)
        return r
    def test_identical_acceptances_have_zero_difference_uncertainty(self):
        r=self.calculate([0.]*4)
        for value in r["paired_difference"].values():
            self.assertAlmostEqual(value["errors"][0],0.,places=8)
            self.assertEqual(value["values"][0],0.)
    def test_subset_error_matches_independent_poisson_gradients(self):
        self.calculate([40.,50.,10.,20.])
    def test_control_descriptor_keeps_primary_paths_and_suppresses_doubled_galleries(self):
        m=c.get_measurement("HERMES_2007_I726689")
        for sel,out in m["outputs"].items():
            if not sel.startswith("RingControl_"):
                self.assertEqual(out["acceptance"],"rectangle_intersect_polar_ring")
                ring=m["outputs"]["RingControl_"+sel]
                self.assertTrue(ring["summary_only"])
                self.assertEqual(ring["paired_primary"],sel)

if __name__ == "__main__": unittest.main()
