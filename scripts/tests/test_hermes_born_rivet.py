#!/usr/bin/env python3
"""Execute canonical Rivet fills on controlled weights and acceptance fixtures."""
from __future__ import annotations

import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "analyses/rivet/dis/HERMES_2007_I726689.cc"


def cpp_definition(source, name):
    match = re.search(rf"    (?:static )?(?:double|void|bool|size_t) {re.escape(name)}\(", source)
    if match is None:
        raise AssertionError(f"Canonical C++ function {name} was not found")
    opening = source.index("{", match.start())
    depth, end = 1, opening + 1
    while depth:
        depth += (source[end] == "{") - (source[end] == "}")
        end += 1
    return source[match.start():end]


class CanonicalBornFillTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compiler = shutil.which("c++") or shutil.which("g++") or shutil.which("clang++")
        if compiler is None:
            raise unittest.SkipTest("A C++ compiler is unavailable for the canonical Born-fill audit")
        temporary = tempfile.TemporaryDirectory(prefix="hermes-Born-audit-")
        cls.addClassCleanup(temporary.cleanup)
        directory = Path(temporary.name)
        source = SOURCE.read_text(encoding="utf-8")
        code = r'''#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>
struct DISKinematicsView { double x=.2; };
struct Hist {
  double eventWeight=0, sumW=0, sumW2=0;
  std::vector<double> coordinates, factors;
  void fill(double coordinate, double factor=1) {
    const double w=eventWeight*factor; sumW+=w; sumW2+=w*w;
    coordinates.push_back(coordinate); factors.push_back(factor);
  }
};
using Histo1DPtr=std::shared_ptr<Hist>;
'''
        code += cpp_definition(source, "fillMeasurement")
        code += r'''
int main() {
  std::cout << std::setprecision(17);
  const double Ds[] = {.2, 0., -1., std::numeric_limits<double>::quiet_NaN()};
  for (double d : Ds) {
    auto ordinary=std::make_shared<Hist>();
    auto weighted=std::make_shared<Hist>();
    auto covariance=std::make_shared<Hist>();
    for (double w : {7.,3.,4.,8.}) {
      ordinary->eventWeight=weighted->eventWeight=covariance->eventWeight=w;
      fillMeasurement(DISKinematicsView{}, d, ordinary, weighted, covariance);
    }
    std::cout << ordinary->sumW << " " << ordinary->sumW2 << " "
              << ordinary->coordinates.size() << " " << weighted->sumW << " "
              << weighted->sumW2 << " " << covariance->sumW2 << "\n";
  }
}
'''
        cpp = directory / "canonical-hermes-born.cc"
        executable = directory / "canonical-hermes-born"
        cpp.write_text(code, encoding="utf-8")
        compilation = subprocess.run([compiler, "-std=c++17", str(cpp), "-o", str(executable)],
                                     capture_output=True, text=True, timeout=60)
        if compilation.returncode:
            raise AssertionError(f"Canonical Born-fill shim failed to compile:\n{compilation.stderr}")
        execution = subprocess.run([str(executable)], capture_output=True, text=True, check=True, timeout=30)
        cls.rows = [tuple(map(float, row.split())) for row in execution.stdout.splitlines()]

    def test_direct_helicity_counts_survive_zero_negative_and_nonfinite_D(self):
        # Independently supplied PP, PM, MP, MM ordinary event weights.
        for row in self.rows:
            with self.subTest(row=row):
                self.assertEqual(row[:3], (22., 138., 4.))
        expected_apar = (7. - 3. - 4. + 8.) / (7. + 3. + 4. + 8.)
        self.assertAlmostEqual(expected_apar, 4. / 11.)

    def test_only_the_A1_conversion_rejects_unusable_D(self):
        valid = self.rows[0]
        self.assertAlmostEqual(valid[3], 22. / .2)
        self.assertAlmostEqual(valid[4], 138. / .2**2)
        self.assertAlmostEqual(valid[5], 138. / .2)
        for row in self.rows[1:]:
            self.assertEqual(row[3:], (0., 0., 0.))
            self.assertTrue(all(math.isfinite(value) for value in row))


class CanonicalBornAcceptanceTests(unittest.TestCase):
    """Run the actual init/analyze helpers with controlled DIS kinematics."""
    @classmethod
    def setUpClass(cls):
        compiler = shutil.which("c++") or shutil.which("g++") or shutil.which("clang++")
        if compiler is None:
            raise unittest.SkipTest("A C++ compiler is unavailable for the canonical acceptance audit")
        temporary = tempfile.TemporaryDirectory(prefix="hermes-Born-acceptance-")
        cls.addClassCleanup(temporary.cleanup)
        directory = Path(temporary.name)
        source = SOURCE.read_text(encoding="utf-8")
        initialization = cpp_definition(source, "init")
        initialization = initialization.replace(
            '      declare(PromptFinalState(Cuts::pid == -11), "PromptPositrons");', '')
        fields = source[source.index("    std::vector<double> _xEdges,"):]
        fields = fields[:fields.index("  };")]
        code = r'''#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>
#define vetoEvent return
struct DISKinematicsView {
  bool valid=true;
  double Q2=2., x=.063, y=.4, W2=10., theta=.1;
};
struct Event { DISKinematicsView dis; };
double fixtureWeight=2.;
struct Hist {
  std::vector<double> edges, sumW, sumW2;
  size_t fills=0;
  Hist(std::vector<double> bins): edges(bins), sumW(bins.size()-1), sumW2(bins.size()-1) {}
  void fill(double coordinate, double factor=1) {
    ++fills;
    if (coordinate<edges.front() || coordinate>=edges.back()) return;
    const size_t bin=std::upper_bound(edges.begin(),edges.end(),coordinate)-edges.begin()-1;
    const double w=fixtureWeight*factor; sumW[bin]+=w; sumW2[bin]+=w*w;
  }
};
using Histo1DPtr=std::shared_ptr<Hist>;
std::string toString(size_t value) {return std::to_string(value);}
std::vector<double> logspace(int count,double low,double high) {
  std::vector<double> result;
  for(int i=0;i<=count;++i) result.push_back(low*std::pow(high/low,double(i)/count));
  return result;
}
class CanonicalAnalysis {
public:
  std::map<std::string,Histo1DPtr> booked;
  double fixtureD=.2;
  void book(Histo1DPtr& hist,const std::string& name,const std::vector<double>& edges) {
    hist=std::make_shared<Hist>(edges); booked[name]=hist;
  }
  void book(Histo1DPtr& hist,const std::string& name,int count,double low,double high) {
    std::vector<double> edges;
    for(int i=0;i<=count;++i) edges.push_back(low+(high-low)*i/count);
    book(hist,name,edges);
  }
  DISKinematicsView disKinematics(const Event& event) {return event.dis;}
  double depolarization(double,double,double) {return fixtureD;}
'''
        code += initialization
        code += cpp_definition(source, "analyze")
        for method in ("fillBornCell", "fillMeasurement", "fillDiagnostics"):
            code += cpp_definition(source, method)
        code += fields + '\n};\n'
        code += r'''
int main() {
  std::cout << std::setprecision(17);
  CanonicalAnalysis init;
  init.init();
  std::cout << "X";
  for(double edge:init._bornXEdges) std::cout << " " << edge;
  std::cout << "\n";
  for(const auto& item:init.booked) {
    if(item.first.find("BornSigmaQ2_")!=0 && item.first.find("SigmaQ2_")!=0) continue;
    std::cout << "E " << item.first;
    for(double edge:item.second->edges) std::cout << " " << edge;
    std::cout << "\n";
  }
  std::set<Hist*> scaled;
  for(const auto& hist:init._scaled) scaled.insert(hist.get());
  std::cout << "S " << init._scaled.size() << " " << scaled.size() << "\n";
  for(int scenario=0;scenario<14;++scenario) {
    CanonicalAnalysis analysis;
    analysis.init();
    Event event;
    if(scenario==1) analysis.fixtureD=std::numeric_limits<double>::quiet_NaN();
    if(scenario==2) analysis.fixtureD=0.;
    if(scenario==3) analysis.fixtureD=-1.;
    if(scenario==4) event.dis.Q2=.999999;
    if(scenario==5) event.dis.Q2=4.;
    if(scenario==6) event.dis.Q2=4.000001;
    if(scenario==7) event.dis.y=.10;
    if(scenario==8) event.dis.W2=3.24;
    if(scenario==9) event.dis.theta=.039999;
    if(scenario==10) event.dis.x=.021239;
    if(scenario==11) event.dis.x=.02124;
    if(scenario==12) event.dis.x=.056809;
    if(scenario==13) event.dis.x=.05681;
    analysis.analyze(event);
    std::cout << "F " << scenario;
    for(const std::string name:{"SigmaX_Q2GT1","SigmaOverD_X_Q2GT1",
      "SigmaQ2_Q2GT1","SigmaQ2_Q2GT4","Accepted_X_Q2GT1",
      "BornSigmaQ2_X04","BornSigmaQ2_X05","BornSigmaQ2_X08","BornSigmaQ2_X09"}) {
      const auto& hist=analysis.booked[name];
      double total=0.; for(double value:hist->sumW) total+=value;
      std::cout << " " << hist->fills << " " << total;
    }
    std::cout << "\n";
  }
}
'''
        cpp = directory / "canonical-hermes-acceptance.cc"
        executable = directory / "canonical-hermes-acceptance"
        cpp.write_text(code, encoding="utf-8")
        compilation = subprocess.run([compiler, "-std=c++17", str(cpp), "-o", str(executable)],
                                     capture_output=True, text=True, timeout=60)
        if compilation.returncode:
            raise AssertionError(f"Canonical acceptance shim failed to compile:\n{compilation.stderr}")
        execution = subprocess.run([str(executable)], capture_output=True, text=True, check=True, timeout=30)
        cls.edges, cls.rows = {}, {}
        for row in execution.stdout.splitlines():
            fields = row.split()
            if fields[0] == "E":
                cls.edges[fields[1]] = list(map(float, fields[2:]))
            elif fields[0] == "F":
                cls.rows[int(fields[1])] = list(map(float, fields[2:]))
            elif fields[0] == "S":
                cls.scaled = tuple(map(int, fields[1:]))
            elif fields[0] == "X":
                cls.xedges = list(map(float, fields[1:]))

    def test_ordinary_cells_and_both_projections_are_independent_of_D(self):
        for scenario in (0, 1, 2, 3):
            row = self.rows[scenario]
            self.assertEqual(row[0:2], [1., 2.])  # x projection
            self.assertEqual(row[4:6], [1., 2.])  # Q2 projection
            self.assertEqual(row[8:10], [1., 2.])  # diagnostic
            self.assertEqual(row[16:18], [1., 2.])  # complete Born cell
            if scenario:
                self.assertEqual(row[2:4], [0., 0.])

    def test_acceptance_and_nested_cut_are_applied_to_direct_observable(self):
        for scenario in (4, 7, 8, 9):
            self.assertEqual(self.rows[scenario], [0.] * 18)
        self.assertEqual(self.rows[5][6:8], [0., 0.])
        self.assertEqual(self.rows[6][6:8], [1., 2.])

    def test_physical_x_cell_edges_are_not_rounded_projection_edges(self):
        # The x=.0212 projection extends below the physical .02124 Born cell.
        self.assertEqual(self.rows[10][0:2], [1., 2.])
        self.assertEqual(self.rows[10][10:18], [0.] * 8)
        self.assertEqual(self.rows[11][12:14], [1., 2.])
        self.assertEqual(self.rows[12][14:16], [1., 2.])
        self.assertEqual(self.rows[12][16:18], [0., 0.])
        self.assertEqual(self.rows[13][14:16], [0., 0.])
        self.assertEqual(self.rows[13][16:18], [1., 2.])

    def test_primary_Q2_cell_cuts_and_explicit_projection_binning(self):
        # Independent primary-analysis inputs, Ehrenfried thesis Appendix C.
        cut_pairs = [(1.505, 2.265), (1.620, 2.623), (1.740, 3.026),
                     (1.882, 3.491), (2.061, 4.032), (2.237, 4.614),
                     (2.657, 5.491), (3.305, 6.645), (4.093, 7.967),
                     (5.043, 9.458), (7.655, 12.528)]
        for index, pair in enumerate(cut_pairs, 9):
            self.assertEqual(self.edges[f"BornSigmaQ2_X{index:02d}"], [1., *pair, 20.])
        for index in range(1, 9):
            expected = [.18, 1.] if index < 5 else [.18, 1., 20.]
            self.assertEqual(self.edges[f"BornSigmaQ2_X{index:02d}"], expected)
        self.assertEqual(self.edges["SigmaQ2_Q2GT1"], [1., 1.5, 2., 3., 4., 6., 8., 12., 20.])
        self.assertEqual(self.edges["SigmaQ2_Q2GT4"], [4., 6., 8., 12., 20.])
        full_cells = sum(sum(low >= 1. for low in edges[:-1]) for name, edges in self.edges.items()
                         if name.startswith("BornSigmaQ2_"))
        self.assertEqual(full_cells, 37)

    def test_runtime_cell_geometry_matches_the_primary_reference_snapshot(self):
        snapshot = json.loads((ROOT / "data/experimental/HERMES_2007_I726689/born-apar-reference.json").read_text())
        geometry = snapshot["binning"]
        self.assertEqual(self.xedges, geometry["physical_x_edges"])
        for index, expected in enumerate(geometry["q2_edges_by_x"], 1):
            self.assertEqual(self.edges[f"BornSigmaQ2_X{index:02d}"], expected)

    def test_every_raw_histogram_is_normalized_exactly_once(self):
        self.assertEqual(self.scaled, (37, 37))


if __name__ == "__main__":
    unittest.main()
