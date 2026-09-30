#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import hermes_unpolarized_fit as fit
import run_experimental_campaign as campaign


class HermesUnpolarizedFitTests(unittest.TestCase):
    def test_primary_fortran_fit_values_for_both_targets(self):
        # Independent fixture: original Lara De Nardo Fortran compiled in /tmp,
        # not values recomputed by this module or from a rounded paper table.
        self.assertEqual(len(fit.SNAPSHOT["crosscheck"]["points"]), 10)
        for point in fit.SNAPSHOT["crosscheck"]["points"]:
            value = fit.f2(point["target"], point["x"], point["q2"])
            self.assertAlmostEqual(value/point["f2"], 1., delta=3e-14)
        self.assertEqual(fit.SNAPSHOT["ALLM97_P_parameters"][6], -.44812)

    def test_vector_R_fits_agree_with_existing_audited_scalar_helpers(self):
        x = np.array([.025, .065, .13, .36, .7])
        q2 = np.array([1.1, 2., 4., 8., 14.])
        for model, function in (("R1998", campaign.r1998), ("R1990", campaign.r1990)):
            np.testing.assert_allclose(fit.longitudinal_ratio(x, q2, model),
                                       [function(a, b) for a, b in zip(x, q2)], rtol=2e-15)

    def test_UU_units_and_equivalent_F1_mass_kinematics(self):
        # Eq37 of the2007 paper provides an independent F1/F2 form of Eq2.6.
        x, q2 = .1, 2.5
        for target in ("P", "D"):
            f2 = fit.f2(target, x, q2)
            mass = fit.DEFAULT_CUTS["nucleon_mass_GeV"]
            energy = fit.DEFAULT_CUTS["beam_energy_GeV"]
            y = q2/(2*mass*energy*x)
            gamma2 = 4*mass*mass*x*x/q2
            f1 = f2*(1+gamma2)/(2*x*(1+campaign.r1998(x, q2)))
            expected = 4*math.pi*fit.SNAPSHOT["constants"]["alpha_em"]**2/q2**2
            expected *= f1*y*y + f2/x*(1-y-y*y*gamma2/4)
            expected *= 389379323.
            self.assertAlmostEqual(fit.differential_cross_section(target, x, q2)/expected, 1., delta=1e-14)

    def test_HERMES_geometry_boundaries_match_direct_kinematics(self):
        energy = fit.DEFAULT_CUTS["beam_energy_GeV"]
        mass = fit.DEFAULT_CUTS["nucleon_mass_GeV"]
        for x in np.linspace(.0212, .9, 51):
            low, high = fit.q2_limits(x)
            for q2 in np.linspace(1, 20, 151):
                y = q2/(2*mass*energy*x)
                w2 = mass*mass + q2*(1/x-1)
                direct = False
                if 0 < y < 1:
                    theta = 2*math.asin(math.sqrt(q2/(4*energy*energy*(1-y)))) if q2/(4*energy*energy*(1-y)) <= 1 else math.inf
                    direct = .1 < y <= .91 and w2 > 3.24 and .04 <= theta <= .22
                analytic = low < q2 < high
                if min(abs(q2-low), abs(q2-high)) > 1e-12:
                    self.assertEqual(direct, analytic)

    def test_integral_Jacobian_and_disjoint_cell_additivity(self):
        with patch.object(fit, "differential_cross_section", side_effect=lambda target, x, q2, *args: np.ones_like(q2)):
            self.assertAlmostEqual(fit.integrate_cross_section("P", .1, .105, 2., 2.1), .005*.1, delta=1e-16)
        total = fit.integrate_cross_section("D", .0212, .9, 1, 20)
        pieces = sum(fit.integrate_cross_section("D", a, b, c, d)
                     for a, b in ((.0212, .1), (.1, .3), (.3, .9))
                     for c, d in ((1., 2.), (2., 4.), (4., 20.)))
        self.assertAlmostEqual(pieces/total, 1., delta=2e-12)
        self.assertEqual(fit.integrate_cross_section("P", .5, .7, 1., 1.1), 0.)
        self.assertEqual(fit.integrate_cross_section("P", .001, .01, 1., 20.), 0.)

    def test_all_supported_cells_and_splits_converge_and_track_fit_domain(self):
        snapshot = json.loads((ROOT/"data/experimental/HERMES_2007_I726689/born-apar-reference.json").read_text())
        for dataset in snapshot["datasets"]:
            for point in dataset["points"]:
                if not point["supported_prediction"]:
                    continue
                for cut in (1., 4.):
                    q2low = max(point["q2_low"], cut)
                    q2high = point["q2_high"]
                    if q2high <= q2low:
                        continue
                    args = (dataset["target"], dataset["x_low"], dataset["x_high"], q2low, q2high)
                    a = fit.integrate_cross_section(*args, quadrature_order=24)
                    b = fit.integrate_cross_section(*args, quadrature_order=48)
                    if a == 0.:
                        self.assertEqual(cut, 4.)
                        self.assertEqual(b, 0.)
                        continue
                    self.assertGreater(a, 0.)
                    self.assertAlmostEqual(a/b, 1., delta=3e-11)
                    control = fit.integrate_cross_section(*args, w2_min=4.)
                    self.assertGreaterEqual(control, 0.)
                    self.assertLessEqual(control, a*(1+1e-12))

    def test_controls_change_shapes_and_raw_covariance_is_not_promoted(self):
        for target in ("P", "D"):
            audit = fit.SNAPSHOT["models"][target]["covariance_audit"]
            self.assertFalse(audit["used_for_errors"])
            self.assertFalse(audit["positive_semidefinite"])
            cov = np.asarray(fit.SNAPSHOT["models"][target]["covariance"])
            self.assertEqual(cov.shape, (23, 23))
            np.testing.assert_array_equal(cov, cov.T)
            self.assertLess(np.linalg.eigvalsh(cov)[0], -1e-8)
            a = fit.integrate_cross_section(target, .06, .08, 1., 4.)
            b = fit.integrate_cross_section(target, .06, .08, 1., 4., model="ALLM97_HYBRID")
            c = fit.integrate_cross_section(target, .06, .08, 1., 4., r_model="R1990")
            self.assertNotAlmostEqual(a, b, delta=a*1e-4)
            self.assertNotAlmostEqual(a, c, delta=a*1e-7)
        hybrid_ratio = fit.f2("D", .08, 2., "ALLM97_HYBRID")/fit.f2("P", .08, 2., "ALLM97_HYBRID")
        self.assertAlmostEqual(hybrid_ratio, fit.f2("D", .08, 2.)/fit.f2("P", .08, 2.), delta=1e-14)

    def test_invalid_inputs_fail_without_silent_model_or_cut_substitution(self):
        for call in (lambda: fit.f2("N", .1, 2.), lambda: fit.f2("P", 0., 2.),
                     lambda: fit.integrate_cross_section("P", .1, .2, 1., 2., model="made_up"),
                     lambda: fit.integrate_cross_section("P", .2, .1, 1., 2.),
                     lambda: fit.integrate_cross_section("P", .1, .2, 1., 2., w2_min=3.),
                     lambda: fit.integrate_cross_section("P", .1, .2, 1., 2., quadrature_order=2)):
            with self.assertRaises(ValueError):
                call()


if __name__ == "__main__":
    unittest.main()
