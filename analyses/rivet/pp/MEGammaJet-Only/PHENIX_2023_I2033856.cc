// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include <cmath>
#include <vector>

namespace Rivet {

  /// PHENIX inclusive and isolated prompt photons at 510 GeV.
  class PHENIX_2023_I2033856 : public Analysis {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(PHENIX_2023_I2033856);

    void init() {
      declare(FinalState(), "FS");
      declare(PromptFinalState(Cuts::abspid == PID::PHOTON &&
                               Cuts::abseta < 0.25 && Cuts::pT > 6*GeV &&
                               Cuts::pT < 30*GeV), "PromptPhotons");
      const std::vector<double> crossSectionEdges =
        {6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,12.0,14.0,16.0,
         18.0,20.0,22.0,24.0,26.0,28.0,30.0};
      const std::vector<double> asymmetryEdges =
        {6.0,7.0,8.0,9.0,10.0,12.0,15.0,20.0};
      book(_inclusive, "InclusivePhotonCrossSection", crossSectionEdges);
      book(_isolated, "IsolatedPhotonCrossSection", crossSectionEdges);
      book(_isolatedYield, "IsolatedPhotonYield", asymmetryEdges);
    }

    void analyze(const Event& event) {
      const Particles photons = apply<PromptFinalState>(event, "PromptPhotons").particles();
      const Particles finalState = apply<FinalState>(event, "FS").particles();
      for (const Particle& photon : photons) {
        const double pt = photon.pT()/GeV;
        _inclusive->fill(pt, 1.0/pt);
        if (!isolated(photon, finalState)) continue;
        _isolated->fill(pt, 1.0/pt);
        _isolatedYield->fill(pt);
      }
    }

    void finalize() {
      if (sumW() == 0.0) return;
      const double sf = crossSection()/picobarn/sumW();
      // E d3sigma/dp3 = (2 pi pT Delta eta)^-1 d sigma/dpT.
      scale(_inclusive, sf/(2.0*M_PI*0.5));
      scale(_isolated, sf/(2.0*M_PI*0.5));
      scale(_isolatedYield, sf);
    }

  private:
    static bool isolated(const Particle& photon, const Particles& finalState) {
      double coneEnergy = 0.0;
      for (const Particle& particle : finalState) {
        if (particle.genParticle() == photon.genParticle()) continue;
        if (PID::isNeutrino(particle.pid()) || deltaR(photon, particle) >= 0.5) continue;
        if (PID::isCharged(particle.pid())) {
          if (particle.pT()/GeV < 0.2) continue;
        } else if (particle.E()/GeV < 0.3) continue;
        coneEnergy += particle.E();
      }
      return coneEnergy < 0.10*photon.E();
    }

    Histo1DPtr _inclusive, _isolated, _isolatedYield;
  };

  RIVET_DECLARE_PLUGIN(PHENIX_2023_I2033856);
}
