// -*- C++ -*-
#include "COMPASSSIDIS.hh"
#include "SIDISAzimuthal.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include <array>
#include <map>
#include <string>
#include <vector>

namespace Rivet {

  /// Low-pT inverse slopes of unidentified charged-hadron multiplicities.
  class COMPASS_2013_I1236358 : public Analysis {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(COMPASS_2013_I1236358);

    struct DISCell {
      double xLow, xHigh, q2Low, q2High;
    };

    void init() {
      declare(FinalState(), "LabFS");
      declare(PromptFinalState(Cuts::pid == -13), "PromptMuons");
      _zEdges = {.20,.25,.30,.35,.40,.50,.60,.70,.80};
      _pt2Edges.reserve(37);
      for (size_t index = 0; index < 36; ++index)
        _pt2Edges.push_back(.01 + .02*double(index));
      _pt2Edges.push_back(.7225);
      const size_t cells = disCells().size()*8*36;
      for (const std::string charge : {"hplus", "hminus"}) {
        book(_yield[charge], "HadronYield_" + charge + "_pt2_cells",
             cells, 0.0, double(cells));
        _scaled.push_back(_yield[charge]);
      }
      book(_acceptedX, "Accepted_X", logspace(40, .0045, .12));
      book(_acceptedQ2, "Accepted_Q2", logspace(40, 1.0, 10.0));
      book(_acceptedY, "Accepted_Y", 32, .1, .9);
      book(_acceptedZ, "Accepted_Z", 30, .2, .8);
      book(_acceptedPt2, "Accepted_HadronPt2", 36, .01, .7225);
      _scaled.insert(_scaled.end(), {
        _acceptedX, _acceptedQ2, _acceptedY, _acceptedZ, _acceptedPt2});
    }

    void analyze(const Event& event) {
      const Particles prompt =
        apply<PromptFinalState>(event, "PromptMuons").particles();
      const COMPASSSIDIS::DISKinematics dis =
        COMPASSInclusiveDIS::reconstruct(event, prompt);
      if (!dis.valid || !COMPASSSIDIS::fixed160GeVBeam(dis)) vetoEvent;
      if (!(dis.Q2 > 1.0 && dis.W2 > 25.0 &&
            dis.y > .1 && dis.y < .9)) vetoEvent;
      const int disCell = findDISCell(dis.x, dis.Q2);
      if (disCell < 0) vetoEvent;

      _acceptedX->fill(dis.x);
      _acceptedQ2->fill(dis.Q2);
      _acceptedY->fill(dis.y);
      std::map<std::string, std::vector<double>> eventYield = {
        {"hplus", std::vector<double>(disCells().size()*8*36, 0.0)},
        {"hminus", std::vector<double>(disCells().size()*8*36, 0.0)},
      };
      for (const Particle& particle :
           apply<FinalState>(event, "LabFS").particles()) {
        if (!PID::isHadron(particle.pid()) || particle.charge3() == 0) continue;
        const auto hadron = COMPASSSIDIS::hadronKinematics(dis, particle, true);
        if (!hadron.valid) continue;
        const int zBin = SIDISAzimuthal::binIndex(hadron.z, _zEdges);
        const int pt2Bin =
          SIDISAzimuthal::binIndex(hadron.transverseMomentum2, _pt2Edges);
        if (zBin < 0 || pt2Bin < 0) continue;
        const size_t flat = (size_t(disCell)*8 + size_t(zBin))*36
                          + size_t(pt2Bin);
        const std::string charge = particle.charge3() > 0 ? "hplus" : "hminus";
        eventYield[charge][flat] += 1.0;
        _acceptedZ->fill(hadron.z);
        _acceptedPt2->fill(hadron.transverseMomentum2);
      }
      for (const auto& entry : eventYield)
        SIDISAzimuthal::fillEventBins(_yield.at(entry.first), entry.second);
    }

    void finalize() {
      if (sumW() == 0.0) return;
      const double factor = crossSection()/picobarn/sumW();
      for (Histo1DPtr& histogram : _scaled) scale(histogram, factor);
    }

  private:
    static const std::array<DISCell, 23>& disCells() {
      static const std::array<DISCell, 23> cells = {{
        {.0045,.0060,1.0,1.25}, {.0060,.0080,1.0,1.30},
        {.0060,.0080,1.3,1.70}, {.0080,.0120,1.0,1.50},
        {.0080,.0120,1.5,2.10}, {.0120,.0180,1.0,1.50},
        {.0120,.0180,1.5,2.50}, {.0120,.0180,2.5,3.50},
        {.0180,.0250,1.0,1.50}, {.0180,.0250,1.5,2.50},
        {.0180,.0250,2.5,3.50}, {.0180,.0250,3.5,5.00},
        {.0250,.0350,1.0,1.20}, {.0250,.0400,1.2,1.50},
        {.0250,.0400,1.5,2.50}, {.0250,.0400,2.5,3.50},
        {.0250,.0400,3.5,6.00}, {.0400,.0500,1.5,2.50},
        {.0400,.0700,2.5,3.50}, {.0400,.0700,3.5,6.00},
        {.0400,.0700,6.0,10.0}, {.0700,.1200,3.5,6.00},
        {.0700,.1200,6.0,10.0},
      }};
      return cells;
    }

    static int findDISCell(double x, double q2) {
      for (size_t index = 0; index < disCells().size(); ++index) {
        const DISCell& cell = disCells()[index];
        if (x >= cell.xLow && x < cell.xHigh &&
            q2 >= cell.q2Low && q2 < cell.q2High)
          return int(index);
      }
      return -1;
    }

    std::vector<double> _zEdges, _pt2Edges;
    std::map<std::string, Histo1DPtr> _yield;
    std::vector<Histo1DPtr> _scaled;
    Histo1DPtr _acceptedX, _acceptedQ2, _acceptedY;
    Histo1DPtr _acceptedZ, _acceptedPt2;
  };

  RIVET_DECLARE_PLUGIN(COMPASS_2013_I1236358);
}
