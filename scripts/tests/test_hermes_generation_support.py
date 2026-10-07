#!/usr/bin/env python3
"""Check the generation envelope from on-shell kinematics, not cut equality.

This regression complements (and cannot replace) the live Herwig low-W2 smoke.
"""
from __future__ import annotations

import math
from pathlib import Path
import re
import unittest

ROOT = Path(__file__).resolve().parents[2]
CARD = ROOT / "cards/experimental/HERMES_2007_I726689/HERMES_2007_I726689-Common.in"


def generator_variables(mass, x, q2):
    """Reconstruct q and an on-shell collinear Born parton in the target rest frame."""
    energy = 27.6
    nu = q2 / (2 * mass * x)
    scattered = energy - nu
    cosine = 1 - q2 / (2 * energy * scattered)
    if not -1 <= cosine <= 1:
        return None
    theta = math.acos(cosine)
    q0, qz = nu, energy - scattered * cosine
    qx = -scattered * math.sin(theta)
    # p=(epsilon,0,0,-epsilon): (p+q)^2=0 fixes epsilon independently.
    parton_energy = q2 / (2 * (q0 + qz))
    outgoing_mass2 = (parton_energy + q0)**2 - qx*qx - (qz-parton_energy)**2
    xi = 2 * parton_energy / mass  # incoming target minus light-cone component is M
    invariant_s = mass*mass + 2*mass*energy
    return {
        "theta": theta,
        "y": nu / energy,
        "W2": (mass+q0)**2-qx*qx-qz*qz,
        "parton_mass2": outgoing_mass2,
        "W2_gen": q2 * (1-xi) / xi,
        "y_gen": q2 / (invariant_s*xi),
    }


class HERMESGenerationEnvelopeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        card = CARD.read_text()
        cls.cuts = {name: float(re.search(rf":{name} ([0-9.e+-]+)", card)[1])
                    for name in ("MinQ2", "MaxQ2", "MinW2", "Miny", "Maxy")}

    def test_all_physical_boundary_samples_fit_both_target_generation_envelopes(self):
        # Concentrate on the physical W2 and y lower boundaries as well as a
        # broad grid. Compute the variables using independent four-vectors.
        counts = {}
        for mass in (.9382720813, .9395654133):
            checked = 0
            for ix in range(1, 601):
                x = .0212 + (.9-.0212)*ix/601
                q2_candidates = [1+(20-1)*i/150 for i in range(151)]
                q2_candidates += [(3.24+epsilon-mass*mass)*x/(1-x)
                                  for epsilon in (1e-8, 1e-4, .01, .1)]
                q2_candidates += [2*mass*27.6*x*y for y in (.10000001, .1001, .101, .91)]
                for q2 in q2_candidates:
                    if not 1 <= q2 <= 20:
                        continue
                    dis = generator_variables(mass, x, q2)
                    if dis is None or not (.1 < dis["y"] <= .91 and dis["W2"] > 3.24
                                           and .04 <= dis["theta"] <= .22):
                        continue
                    checked += 1
                    self.assertAlmostEqual(dis["parton_mass2"], 0., delta=1e-9)
                    self.assertGreaterEqual(dis["W2_gen"], self.cuts["MinW2"])
                    self.assertGreaterEqual(dis["y_gen"], self.cuts["Miny"])
                    self.assertLessEqual(dis["y_gen"], self.cuts["Maxy"])
            counts[mass] = checked
        self.assertTrue(all(count > 10000 for count in counts.values()), counts)

    def test_pre_v5_cuts_reject_valid_high_x_and_low_y_events(self):
        for mass in (.9382720813, .9395654133):
            high_x = generator_variables(mass, .7, 6.)
            self.assertGreater(high_x["W2"], 3.24)
            self.assertLess(high_x["W2_gen"], 3.24)
            self.assertGreater(high_x["W2_gen"], self.cuts["MinW2"])
            low_y = generator_variables(mass, .3, 2*mass*27.6*.3*.1005)
            self.assertGreater(low_y["W2_gen"], 3.24)
            self.assertGreater(low_y["y"], .1)
            self.assertLess(low_y["y_gen"], .1)
            self.assertGreater(low_y["y_gen"], self.cuts["Miny"])

    def test_envelope_keeps_margin_and_the_original_Q2_support(self):
        self.assertEqual(self.cuts, {"MinQ2": 1., "MaxQ2": 20.,
                                    "MinW2": 2., "Miny": .095, "Maxy": .91})
        for mass in (.9382720813, .9395654133):
            self.assertGreater(.1/(1+mass/(2*27.6)), self.cuts["Miny"] + .003)
            self.assertGreater(3.24-mass*mass+mass/(2*27.6), self.cuts["MinW2"] + .37)


if __name__ == "__main__":
    unittest.main()
