#!/usr/bin/env python3
"""Independent HERMES closure for one eventwise D correction and one D-state factor.

The histogram fixtures originate from UU and virtual-photon A1 strata with
different D, rather than preweighted input histograms. Expected ratio errors
come from derivatives with respect to independent Poisson event counts.
"""
from __future__ import annotations

import copy
from fractions import Fraction as F
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import types
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import run_experimental_campaign as campaign

MEASUREMENT = "HERMES_2007_I726689"
SOURCE = ROOT / "analyses" / "rivet" / "dis" / f"{MEASUREMENT}.cc"
SNAPSHOT = ROOT / "data" / "experimental" / MEASUREMENT / "reference.json"
SIGNS = {"PP": 1, "PM": -1, "MP": -1, "MM": 1}
STRATA = {
    # (unpolarized cross section, photon A1, longitudinal depolarization D).
    "P": [(F(8), F(3, 10), F(1, 4)), (F(2), F(3, 5), F(3, 4))],
    "N": [(F(3), F(-1, 5), F(2, 5)), (F(17), F(1, 10), F(4, 5))],
}
EVENTS = {
    (component, helicity): [
        (uu + sign * uu * d * a1, d) for uu, a1, d in strata
    ]
    for component, strata in STRATA.items() for helicity, sign in SIGNS.items()
}
KINEMATICS = ((.0264, .7, 1.12), (.173, .4, 4.31), (.4583, .3, 8.53),
              (.5819, .5, 10.16), (.7248, .32, 12.21))


def epsilon_depolarization(x, y, q2):
    """HERMES epsilon definition and D=[1-(1-y)epsilon]/[1+epsilon R]."""
    # Independent coefficient inputs from Whitlow R1990 fit B. Do not call
    # either production R or production D while constructing the expectation.
    fit_b = (.0635, .5747, -.3534)
    theta = 1. + 12. * q2 / (q2 + 1.) * .015625 / (.015625 + x**2)
    r = fit_b[0] * theta / math.log(q2 / .04) + fit_b[1] / q2 + fit_b[2] / (q2**2 + .09)
    gamma_squared = 4. * .9382720813**2 * x**2 / q2
    epsilon = ((1. - y - gamma_squared * y**2 / 4.)
               / (1. - y + y**2 / 2. + gamma_squared * y**2 / 4.))
    return r, (1. - (1. - y) * epsilon) / (1. + epsilon * r)


def event_moments():
    """Raw Poisson moments, before target mixing or asymmetry construction."""
    return {
        label: (sum(w for w, d in rows), sum(w*w for w, d in rows),
                sum(w/d for w, d in rows), sum(w*w/(d*d) for w, d in rows),
                sum(w*w/d for w, d in rows))
        for label, rows in EVENTS.items()
    }


def independent_prediction(target_weights):
    """Compute means from physical strata, and errors from count derivatives."""
    denominator = sum(u * sum(row[0] for row in STRATA[c])
                      for c, (u, l) in target_weights.items())
    longitudinal = sum(l * sum(uu*d*a1 for uu, a1, d in STRATA[c])
                       for c, (u, l) in target_weights.items())
    photon_longitudinal = sum(l * sum(uu*a1 for uu, a1, d in STRATA[c])
                              for c, (u, l) in target_weights.items())
    apar, a1 = longitudinal/denominator, photon_longitudinal/denominator
    variances = {"apar": F(0), "a1": F(0)}
    covariances = {"apar": F(0), "a1": F(0)}
    for (component, helicity), rows in EVENTS.items():
        u, l = target_weights[component]
        a, b = u/4, l*SIGNS[helicity]/4
        for w, d in rows:
            # d(N/M)/dn_i = (dN/dn_i - (N/M)dM/dn_i)/M.
            variances["apar"] += ((b*w-apar*a*w)/denominator)**2
            variances["a1"] += ((b*w/d-a1*a*w)/denominator)**2
            covariances["apar"] += a*b*w*w
            covariances["a1"] += a*b*w*w/d
    return {"uu": denominator, "ll": longitudinal, "ll_over_d": photon_longitudinal,
            "apar": apar, "a1": a1, "variances": variances, "covariances": covariances}


EXPECTED = {
    "Q2GT1": independent_prediction({"P": (F(1), F(1)), "N": (F(0), F(0))}),
    "D_Q2GT1": independent_prediction({"P": (F(1, 2), F(37, 80)), "N": (F(1, 2), F(37, 80))}),
}


def postprocess_moments(moments, measurement=None):
    measurement = copy.deepcopy(measurement or campaign.get_measurement(MEASUREMENT))
    measurement["diagnostics"] = []
    family = campaign.campaign_family_specs(measurement, False)["nominal"]

    def load(groups, directory, analysis, name, family_id, component, helicities, orders, *variation):
        result = {}
        for helicity in helicities:
            ordinary, ordinary_variance, weighted, weighted_variance, covariance = moments[component, helicity]
            if name.startswith("SigmaOverD_"):
                value, variance = weighted, weighted_variance
            elif name.startswith("CovarianceProxy_"):
                value, variance = 0., covariance
            else:
                value, variance = ordinary, ordinary_variance
            result[helicity] = campaign.BinSeries([.1, .3], [float(value)], [float(variance)])
        return result

    with mock.patch.object(campaign, "_load_component_series", side_effect=load):
        results, _ = campaign._load_asymmetry_family_products({}, ROOT, measurement, "nominal", family)
    return results


def assert_closure(testcase, results):
    # Known values are independent of the histogram/postprocessor implementation.
    testcase.assertEqual(EXPECTED["Q2GT1"]["apar"], F(3, 20))
    testcase.assertEqual(EXPECTED["Q2GT1"]["a1"], F(9, 25))
    testcase.assertEqual(EXPECTED["D_Q2GT1"]["apar"], F(4847, 60000))
    testcase.assertEqual(EXPECTED["D_Q2GT1"]["a1"], F(1739, 12000))
    for selection, expected in EXPECTED.items():
        actual = results[selection]
        comparisons = {
            "uu": actual["sigma_uu"].values[0],
            "ll": actual["sigma_ll"].values[0],
            "ll_over_d": actual["sigma_ll_over_d"].values[0],
            "apar": actual["apar_values"][0], "a1": actual["a1_values"][0],
        }
        for key, value in comparisons.items():
            testcase.assertAlmostEqual(value, float(expected[key]), places=13, msg=f"{selection} {key}")
        for observable in ("apar", "a1"):
            testcase.assertAlmostEqual(actual[f"{observable}_errors"][0]**2,
                                       float(expected["variances"][observable]), places=13)
    # Identical test strata feed both nested selections; the target correction
    # must also be identical when the Q2>4 selection is postprocessed.
    for primary, conservative in (("Q2GT1", "Q2GT4"), ("D_Q2GT1", "D_Q2GT4")):
        testcase.assertEqual(results[primary]["a1_values"], results[conservative]["a1_values"])


class DepolarizationArithmeticTests(unittest.TestCase):
    def test_python_D_agrees_with_independent_epsilon_definition(self):
        for x, y, q2 in KINEMATICS:
            with self.subTest(x=x, y=y, q2=q2):
                expected_r, expected_d = epsilon_depolarization(x, y, q2)
                self.assertAlmostEqual(campaign.r1990(x, q2), expected_r, places=14)
                self.assertAlmostEqual(campaign.depolarization(x, y, q2), expected_d, places=14)

    def test_varying_D_physical_strata_and_poisson_ratio_variances(self):
        assert_closure(self, postprocess_moments(event_moments()))

    def test_covariance_uses_one_D_inverse_and_UU_times_LL_coefficients(self):
        moments = event_moments()
        measurement = campaign.get_measurement(MEASUREMENT)
        family = campaign.campaign_family_specs(measurement, False)["nominal"]
        for selection, expected in EXPECTED.items():
            uu_target = campaign._target_component_coefficients(measurement, family, "unpolarized", selection)
            ll_target = campaign._target_component_coefficients(measurement, family, "longitudinal", selection)
            uu = {f"{c}:{h}": uu_target[c]/4. for c, h in moments}
            ll = {f"{c}:{h}": ll_target[c]*SIGNS[h]/4. for c, h in moments}
            for observable, variance_index in (("apar", 1), ("a1", 4)):
                proxies = {
                    f"{c}:{h}": campaign.BinSeries([.1, .3], [0.], [float(values[variance_index])])
                    for (c, h), values in moments.items()
                }
                covariance = campaign._covariance_array(proxies, ll, uu)[0]
                self.assertAlmostEqual(covariance, float(expected["covariances"][observable]), places=13)
        self.assertEqual(EXPECTED["D_Q2GT1"]["covariances"]["a1"], F(18019, 3200))

    def test_double_D_mean_D_and_double_Dstate_are_rejected(self):
        expected = EXPECTED["D_Q2GT1"]
        double_d = F(37, 80) * sum(uu*a1/d for rows in STRATA.values() for uu, a1, d in rows) / expected["uu"]
        average_d = sum(uu*d for rows in STRATA.values() for uu, a1, d in rows) / F(30)
        wrong = (expected["apar"], double_d, expected["apar"]/average_d, F(37, 40)*expected["a1"])
        actual = postprocess_moments(event_moments())["D_Q2GT1"]["a1_values"][0]
        for value in wrong:
            self.assertFalse(math.isclose(actual, float(value), rel_tol=1e-10, abs_tol=1e-10))
        # d_state_factor is explanatory metadata. The numerical correction
        # resides in the .4625 component coefficients, rather than both fields.
        measurement = campaign.get_measurement(MEASUREMENT)
        measurement["outputs"]["D_Q2GT1"]["d_state_factor"] = .31
        self.assertAlmostEqual(postprocess_moments(event_moments(), measurement)["D_Q2GT1"]["a1_values"][0], actual)


def cpp_function(source, name):
    match = re.search(rf"    static (?:double|void) {re.escape(name)}\(", source)
    if match is None:
        raise AssertionError(f"Canonical C++ function {name} was not found")
    opening = source.index("{", match.start())
    depth = 1
    end = opening + 1
    while depth:
        depth += (source[end] == "{") - (source[end] == "}")
        end += 1
    return source[match.start():end]


class CanonicalCppDepolarizationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compiler = shutil.which("c++") or shutil.which("g++") or shutil.which("clang++")
        if compiler is None:
            raise unittest.SkipTest("A C++ compiler is unavailable for the canonical-fill audit")
        temporary = tempfile.TemporaryDirectory(prefix="hermes-D-audit-")
        cls.addClassCleanup(temporary.cleanup)
        directory = Path(temporary.name)
        source = SOURCE.read_text(encoding="utf-8")
        code = '''#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
struct DISKinematicsView { double x=.2; };
struct Hist {
  double eventWeight=0, sumW=0, sumW2=0;
  void fill(double, double factor=1) {
    const double w=eventWeight*factor; sumW+=w; sumW2+=w*w;
  }
};
using Histo1DPtr=std::shared_ptr<Hist>;
'''
        code += "\n".join(cpp_function(source, name) for name in ("r1990", "depolarization", "fillMeasurement"))
        code += '\nint main() {\nstd::cout << std::setprecision(17);\n'
        for index, (x, y, q2) in enumerate(KINEMATICS):
            code += (f'std::cout << "D {index} " << r1990({x:.17g},{q2:.17g}) << " " '
                     f'<< depolarization({x:.17g},{y:.17g},{q2:.17g}) << "\\n";\n')
        for (component, helicity), rows in EVENTS.items():
            code += '{ auto o=std::make_shared<Hist>(); auto w=std::make_shared<Hist>(); auto c=std::make_shared<Hist>();\n'
            for weight, d in rows:
                code += f'o->eventWeight=w->eventWeight=c->eventWeight={float(weight):.17g};\n'
                code += f'fillMeasurement(DISKinematicsView{{}}, {float(d):.17g}, o, w, c);\n'
            code += (f'std::cout << "M {component} {helicity} " << o->sumW << " " << o->sumW2 '
                     '<< " " << w->sumW << " " << w->sumW2 << " " << c->sumW2 << "\\n"; }\n')
        code += '}\n'
        cpp = directory / "canonical-hermes-fill.cc"
        executable = directory / "canonical-hermes-fill"
        cpp.write_text(code, encoding="utf-8")
        compilation = subprocess.run([compiler, "-std=c++17", str(cpp), "-o", str(executable)],
                                     capture_output=True, text=True, timeout=60)
        if compilation.returncode:
            raise AssertionError(f"Canonical-fill shim failed to compile:\n{compilation.stderr}")
        execution = subprocess.run([str(executable)], capture_output=True, text=True, check=True, timeout=30)
        cls.depolarizations, cls.moments = {}, {}
        for line in execution.stdout.splitlines():
            fields = line.split()
            if fields[0] == "D":
                cls.depolarizations[int(fields[1])] = tuple(map(float, fields[2:]))
            elif fields[0] == "M":
                cls.moments[fields[1], fields[2]] = tuple(map(float, fields[3:]))

    def test_actual_cpp_D_agrees_with_independent_epsilon_definition(self):
        for index, kinematics in enumerate(KINEMATICS):
            expected_r, expected_d = epsilon_depolarization(*kinematics)
            actual_r, actual_d = self.depolarizations[index]
            self.assertAlmostEqual(actual_r, expected_r, places=14)
            self.assertAlmostEqual(actual_d, expected_d, places=14)

    def test_actual_cpp_fill_has_one_D_inverse_and_closes_through_postprocessing(self):
        # The shim reproduces histogram event-weight multiplication, rather
        # than running Rivet projections. Production bodies are copied verbatim.
        for label, expected in event_moments().items():
            for actual, value in zip(self.moments[label], expected):
                self.assertTrue(math.isclose(actual, float(value), rel_tol=2e-14, abs_tol=2e-14))
        assert_closure(self, postprocess_moments(self.moments))


# Independent transcription of all five columns from the official archive:
# source SHA256 2567ee83422f8da7c3c602e2aa716c2501cc539820ba0f4f9657515ed2785997,
# A1p_a15 SHA256 f9681449a1b478c89986b72fd17a2598899b703fa06c259b0a8781b8a956f700,
# A1d_a15 SHA256 51acfe83a21d0e9950477090882e59fcf7a5c752a1ed7adcd17de4901bc23176.
# These expectations must never be constructed from the reference snapshot.
PUBLISHED_X = [.0264, .0329, .0403, .0506, .0648, .0829, .1059, .1354, .173,
               .2209, .282, .3598, .4583, .5819, .7248]
PUBLISHED_Q2 = [1.12, 1.25, 1.38, 1.54, 2.01, 2.45, 2.97, 3.59, 4.31,
                5.13, 6.11, 7.24, 8.53, 10.16, 12.21]
PUBLISHED = {
    "P": {
        "value": [.1113, .0975, .0897, .1061, .109, .1997, .2011, .2681, .3246,
                  .3361, .4094, .5169, .6573, .6647, 1.1976],
        "stat": [.0351, .0268, .0236, .0201, .0202, .0206, .0218, .0236, .0263,
                 .0305, .037, .0486, .0714, .1302, .2763],
        "systematic_combined": [.0095, .0087, .0082, .0100, .0093, .0147, .0150, .0190,
                                .0223, .0218, .0263, .0333, .0416, .0469, .0827],
    },
    "D": {
        "value": [.0275, -.0123, .0039, .0210, .0254, .0645, .0503, .0907, .1318,
                  .1656, .1810, .3255, .3815, .4403, .8641],
        "stat": [.0167, .0133, .0120, .0104, .0103, .0107, .0115, .0128, .0146,
                 .0173, .0216, .0291, .0440, .0813, .1790],
        "systematic_combined": [.0020, .0007, .0011, .0020, .0023, .0039, .0031, .0049,
                                .0066, .0079, .0085, .0152, .0200, .0281, .0759],
    },
}


class CapturingBin:
    def __init__(self):
        self.value, self.errors = None, {}

    def setVal(self, value):
        self.value = value

    def setErr(self, low, high, label):
        self.errors[label] = (low, high)


class CapturingEstimate:
    def __init__(self, edges, path):
        self.path, self.annotations = path, {}
        self.bins = [CapturingBin() for _ in edges[1:]]

    def setTitle(self, value):
        self.annotations["Title"] = value

    def setAnnotation(self, key, value):
        self.annotations[key] = value

    def bin(self, index):
        return self.bins[index-1]


class PublishedReferenceNormalizationTests(unittest.TestCase):
    def test_all30_snapshot_values_and_errors_equal_published_A1_without_correction(self):
        snapshot = campaign.load_json(SNAPSHOT)
        datasets = {item["id"]: item for item in snapshot["datasets"]}
        self.assertEqual(set(datasets), {"P", "D"})
        self.assertEqual(sum(len(item["points"]) for item in datasets.values()), 30)
        for target, dataset in datasets.items():
            self.assertEqual([point["x_mean"] for point in dataset["points"]], PUBLISHED_X)
            self.assertEqual([point["q2_mean"] for point in dataset["points"]], PUBLISHED_Q2)
            for key, expected in PUBLISHED[target].items():
                self.assertEqual([point[key] for point in dataset["points"]], expected)

    def test_multi_reference_writer_copies_all30_A1_values_and_errors(self):
        snapshot = campaign.load_json(SNAPSHOT)
        yoda = types.SimpleNamespace(BinnedEstimate1D=CapturingEstimate)
        with mock.patch.object(campaign, "_import_yoda", return_value=yoda), \
                mock.patch.object(campaign, "_write_yoda_objects") as write:
            campaign.write_reference_yoda(snapshot, Path("unused-reference.yoda"))
        objects = {item.path: item for item in write.call_args.args[1]}
        for target, token in (("P", "y01"), ("D", "y02")):
            obj = objects[f"/REF/{MEASUREMENT}/d14-x01-{token}"]
            self.assertEqual([item.value for item in obj.bins], PUBLISHED[target]["value"])
            self.assertEqual([item.errors["stat"][1] for item in obj.bins], PUBLISHED[target]["stat"])
            self.assertNotIn("DeuteronDStateFactor", obj.annotations)

    def test_actual_reference_YODA_contains_all30_unchanged_published_A1_values(self):
        try:
            import yoda
        except (ImportError, OSError):
            self.skipTest("YODA bindings are unavailable for actual reference-file inspection")
        objects = yoda.read(str(ROOT / "analyses" / "rivet" / "dis" / f"{MEASUREMENT}.yoda.gz"))
        for target, token in (("P", "y01"), ("D", "y02")):
            obj = objects[f"/REF/{MEASUREMENT}/d14-x01-{token}"]
            self.assertEqual(obj.numBins(), 15)
            for index, (value, stat) in enumerate(zip(PUBLISHED[target]["value"], PUBLISHED[target]["stat"]), 1):
                self.assertAlmostEqual(obj.bin(index).val(), value, places=13)
                self.assertAlmostEqual(obj.bin(index).errAvg("stat"), stat, places=13)


if __name__ == "__main__":
    unittest.main()
