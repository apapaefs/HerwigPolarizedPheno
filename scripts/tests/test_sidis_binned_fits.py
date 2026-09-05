"""Numerical closure and failure cases for the published COMPASS estimators."""
import json
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import compass_multiplicity_fiducial as fiducial
import run_experimental_campaign as experimental
import sidis_binned_fits as fits


def series(values, variances):
    return experimental.BinSeries(list(range(len(values)+1)),
                                  list(values), list(variances))


def packed_covariance(covariance):
    upper = covariance[np.triu_indices(len(covariance), 1)]
    # Deliberately unrelated central values: the second moment is in sumW2.
    return series(np.full(len(upper), 12345.), upper)


def phi_sample(coefficients=(.16, -.08, .05), epsilon=.7, normalization=1000.):
    edges = np.linspace(0., 2.*np.pi, 17)
    centers = .5*(edges[1:]+edges[:-1])
    # Independent analytic integral: sinc(n*width/2) times the center value.
    # Do not generate closure data using the fitter's design-matrix function.
    counts = normalization*(1. + coefficients[0]*np.sinc(1./16.)*np.cos(centers)
                            + coefficients[1]*np.sinc(2./16.)*np.cos(2.*centers)
                            + coefficients[2]*np.sinc(1./16.)*np.sin(centers))
    inputs = np.r_[counts, epsilon*sum(counts)]
    # Independent Poisson phi counts, with the same events in the epsilon sum.
    covariance = np.zeros((17, 17))
    covariance[:16, :16] = np.diag(counts)
    covariance[:16, 16] = covariance[16, :16] = epsilon*counts
    covariance[16, 16] = epsilon**2*sum(counts)
    return inputs, covariance


POLICY = {"phi_edges": np.linspace(0., 2.*np.pi, 17).tolist(),
          "excluded_phi_bins": [0, 15]}


class AzimuthalFitTests(unittest.TestCase):
    def fit(self, inputs, covariance, harmonic=1):
        return fits.azimuthal_target_fit(
            {"P": series(inputs, np.diag(covariance))},
            {"P": packed_covariance(covariance)}, {"P": 1.},
            harmonic, 1, POLICY)

    def test_bin_integrated_cosines_and_sine_nuisance_close(self):
        inputs, covariance = phi_sample()
        for harmonic, expected in ((1, .16/.7), (2, -.08/.7)):
            result = self.fit(inputs, covariance, harmonic)
            self.assertAlmostEqual(result["values"][0], expected, places=12)
            detail = result["fit_diagnostics"][0]
            self.assertEqual(detail["retained_phi_bins"], list(range(1, 15)))
            self.assertEqual(detail["ndof"], 10)
            self.assertAlmostEqual(detail["sin_coefficient_nuisance"], .05)
            self.assertGreater(result["errors"][0], 0.)

    def test_excluded_phi_contamination_does_not_enter_harmonic_fit(self):
        inputs, covariance = phi_sample()
        original = self.fit(inputs, covariance)["values"][0]
        inputs[[0, 15]] += 5000.
        inputs[16] = .7*sum(inputs[:16])
        self.assertAlmostEqual(self.fit(inputs, covariance)["values"][0], original)

    def test_higher_harmonic_uses_truncated_fit_not_full_phi_moment(self):
        edges = np.asarray(POLICY["phi_edges"])
        inputs, covariance = phi_sample((0., 0., 0.))
        cos3 = np.diff(np.sin(3.*edges))/(3.*np.diff(edges))
        inputs[:16] += 400.*cos3
        inputs[16] = .7*sum(inputs[:16])
        design = fits.azimuthal_design(edges)[1:15]
        beta = np.linalg.lstsq(design, inputs[1:15], rcond=None)[0]
        expected = beta[1]/beta[0]/.7
        self.assertGreater(abs(expected), .05)
        self.assertAlmostEqual(self.fit(inputs, covariance)["values"][0], expected)

    def test_shared_epsilon_covariance_matches_finite_difference_gradient(self):
        inputs, covariance = phi_sample()
        # Add independent epsilon fluctuations, preserving positive covariance.
        covariance[16, 16] += 400.
        design = fits.azimuthal_design(POLICY["phi_edges"])
        retained = np.arange(1, 15)
        value, error, _, gradient = fits._azimuthal_fit(
            inputs, covariance, 1, design, retained)
        numerical = []
        for index in range(17):
            plus, minus = inputs.copy(), inputs.copy()
            step = .001
            plus[index] += step
            minus[index] -= step
            numerical.append((fits._azimuthal_fit(plus, covariance, 1, design, retained)[0]
                              - fits._azimuthal_fit(minus, covariance, 1, design, retained)[0])/(2.*step))
        np.testing.assert_allclose(gradient, numerical, atol=1.e-12, rtol=1.e-7)
        self.assertAlmostEqual(error**2, gradient @ covariance @ gradient)
        self.assertTrue(np.isfinite(value))

    def test_target_yields_are_combined_before_fit_and_depolarization(self):
        proton, cp = phi_sample((.1, 0., 0.), .6, 1000.)
        neutron, cn = phi_sample((-.2, 0., 0.), .8, 500.)
        result = fits.azimuthal_target_fit(
            {"P": series(proton, np.diag(cp)), "N": series(neutron, np.diag(cn))},
            {"P": packed_covariance(cp), "N": packed_covariance(cn)},
            {"P": .5, "N": .5}, 1, 1, POLICY)
        self.assertAlmostEqual(result["values"][0], 0., places=12)
        self.assertAlmostEqual(result["fit_diagnostics"][0]["epsilon_mean"], 2./3.)

    def test_missing_variance_in_retained_phi_bin_masks_cell(self):
        inputs, covariance = phi_sample()
        covariance[4, 4] = 0.
        result = self.fit(inputs, covariance)
        self.assertIsNone(result["values"][0])
        self.assertEqual(result["fit_diagnostics"][0]["status"], "missing_variance_in_fit_range")

    def test_target_second_moments_use_squared_coefficients(self):
        samples = {"P": series([3., 5.], [4., 9.]),
                   "N": series([1., 2.], [1., 4.])}
        proxies = {"P": series([999.], [2.]), "N": series([999.], [1.])}
        means, covariance = fits._moments(samples, {"P": .5, "N": 2.}, 2, 1, proxies)
        np.testing.assert_allclose(means, [[3.5, 6.5]])
        np.testing.assert_allclose(covariance, [[[5., 4.5], [4.5, 18.25]]])


class SlopeFitTests(unittest.TestCase):
    edges = np.r_[.01+.02*np.arange(36), .7225]

    def fit(self, yields, covariance):
        return fits.pt2_slope_target_fit(
            {"P": series(yields, np.diag(covariance))}, {"P": 1.}, self.edges, 1,
            {"P": packed_covariance(covariance)})

    def test_integrated_exponential_closes_with_correlated_bins(self):
        yields = 1.e5*fits.exponential_integrals(self.edges, .27)
        covariance = np.diag(yields) + .0001*np.outer(yields, yields)
        result = self.fit(yields, covariance)
        self.assertAlmostEqual(result["values"][0], .27, places=9)
        self.assertEqual(result["retained_pt2_bins"][0], list(range(36)))
        self.assertEqual(result["fit_diagnostics"][0]["ndof"], 34)
        self.assertGreater(result["errors"][0], 0.)

    def test_three_positive_bins_cannot_define_full_range_slope(self):
        yields = np.r_[100., 90., 80., -np.ones(33)]
        result = self.fit(yields, np.eye(36))
        self.assertIsNone(result["values"][0])
        self.assertEqual(result["retained_pt2_bins"][0], [])
        self.assertEqual(result["fit_diagnostics"][0]["status"], "insufficient_signal_across_fit_range")

    def test_signed_bin_is_retained_and_bad_shape_is_reported(self):
        yields = 1.e5*fits.exponential_integrals(self.edges, .27)
        covariance = np.diag(yields)
        yields[17] = -1.
        result = self.fit(yields, covariance)
        self.assertIsNotNone(result["values"][0])
        self.assertEqual(result["fit_diagnostics"][0]["negative_bins_retained"], 1)
        self.assertTrue(result["fit_diagnostics"][0]["poor_exponential_shape"])
        self.assertEqual(result["retained_pt2_bins"][0], list(range(36)))

    def test_empty_and_singular_covariance_are_masked(self):
        yields = 1.e5*fits.exponential_integrals(self.edges, .27)
        empty = np.diag(yields)
        empty[-1, -1] = 0.
        for covariance, status in ((empty, "missing_variance_in_fit_range"),
                                   (np.ones((36, 36)), "singular_covariance")):
            result = self.fit(yields, covariance)
            self.assertIsNone(result["values"][0])
            self.assertEqual(result["fit_diagnostics"][0]["status"], status)


class FixedBeamSupportTests(unittest.TestCase):
    def test_published_cells_without_nominal_beam_support_are_explicit(self):
        expected = {"COMPASS_2017_I1444985": 5, "COMPASS_2017_I1483098": 4,
                    "COMPASS_2025_I2840545": 0, "COMPASS_2026_I3096394": 0}
        for measurement, count in expected.items():
            snapshot = json.loads((ROOT / "data/phenomenology" / measurement / "reference.json").read_text())
            for species, dataset in snapshot["datasets"].items():
                with self.subTest(measurement=measurement, species=species):
                    masked = fiducial.unsupported_cells(snapshot, dataset, species, {"P": .5, "N": .5})
                    self.assertEqual(len(masked), count)

    def test_unrelated_measurements_are_not_masked(self):
        self.assertEqual(fiducial.unsupported_cells({}, {}, "hplus", {"P": 1.}), [])


if __name__ == "__main__":
    unittest.main()
