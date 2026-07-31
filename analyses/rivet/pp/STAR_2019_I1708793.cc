// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include <cmath>
#include <vector>

namespace Rivet {

  /// STAR W+/W- and Z/gamma* longitudinal-spin fiducial yields at 510 GeV.
  class STAR_2019_I1708793 : public Analysis {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(STAR_2019_I1708793);

    void init() {
      declare(FinalState(), "FS");
      const Cut wcut = Cuts::abspid == PID::ELECTRON && Cuts::pT > 25*GeV &&
                       Cuts::pT < 50*GeV && Cuts::abseta < 1.5;
      const Cut zcut = Cuts::abspid == PID::ELECTRON && Cuts::pT > 14*GeV &&
                       Cuts::abseta < 1.1;
      declare(LeptonFinder(0.1, wcut, LeptonOrigin::NODECAY,
                           PhotonOrigin::NODECAY), "WLeptons");
      declare(LeptonFinder(0.1, zcut, LeptonOrigin::NODECAY,
                           PhotonOrigin::NODECAY), "ZLeptons");
      const std::vector<double> etaEdges = {-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5};
      book(_wplus, "Yield_Wplus_Eta", etaEdges);
      book(_wminus, "Yield_Wminus_Eta", etaEdges);
      book(_z, "Yield_Zgamma", std::vector<double>{-0.5,0.5});
    }

    void analyze(const Event& event) {
      const Particles wleptons = apply<LeptonFinder>(event, "WLeptons").particles();
      if (wleptons.size() == 1) {
        if (wleptons[0].pid() == -11) _wplus->fill(wleptons[0].eta());
        if (wleptons[0].pid() ==  11) _wminus->fill(wleptons[0].eta());
      }
      const Particles zleptons = apply<LeptonFinder>(event, "ZLeptons").particles();
      const Particles finalState = apply<FinalState>(event, "FS").particles();
      for (size_t first = 0; first < zleptons.size(); ++first) {
        for (size_t second = first+1; second < zleptons.size(); ++second) {
          if (zleptons[first].pid()*zleptons[second].pid() >= 0) continue;
          const double mass = (zleptons[first].momentum()+zleptons[second].momentum()).mass()/GeV;
          if (mass <= 70.0 || mass >= 110.0) continue;
          if (!isolated(zleptons[first], finalState) ||
              !isolated(zleptons[second], finalState)) continue;
          _z->fill(0.0);
        }
      }
    }

    void finalize() {
      if (sumW() == 0.0) return;
      const double sf = crossSection()/picobarn/sumW();
      scale(_wplus, sf); scale(_wminus, sf); scale(_z, sf);
    }

  private:
    static bool constituent(const Particle& particle, const Particle& dressed) {
      for (const Particle& item : dressed.constituents())
        if (particle.genParticle() == item.genParticle()) return true;
      return false;
    }

    static bool isolated(const Particle& lepton, const Particles& finalState) {
      double coneEt = lepton.Et();
      for (const Particle& particle : finalState) {
        if (PID::isNeutrino(particle.pid()) || constituent(particle, lepton)) continue;
        if (deltaR(lepton, particle) < 0.7) coneEt += particle.Et();
      }
      return coneEt > 0.0 && lepton.Et()/coneEt > 0.88;
    }

    Histo1DPtr _wplus, _wminus, _z;
  };

  RIVET_DECLARE_PLUGIN(STAR_2019_I1708793);
}
