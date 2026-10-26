// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include <algorithm>
#include <cmath>
#include <vector>

namespace Rivet {

  /// PHENIX inclusive and isolated prompt photons at 510 GeV.
  class PHENIX_2023_I2033856 : public Analysis {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(PHENIX_2023_I2033856);

    void init() {
      // Data are corrected for limited detector acceptance. Use a full cone,
      // and truth prompt photons instead of repeating the decay subtraction.
      declare(FinalState(), "FS");
      declare(PromptFinalState(Cuts::abspid == PID::PHOTON &&
                               Cuts::abseta < 0.25 && Cuts::pT > 6*GeV &&
                               Cuts::pT < 30*GeV), "PromptPhotons");
      _crossSectionEdges =
        {6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,12.0,14.0,16.0,
         18.0,20.0,22.0,24.0,26.0,28.0,30.0};
      _asymmetryEdges =
        {6.0,7.0,8.0,9.0,10.0,12.0,15.0,20.0};
      book(_inclusive, "InclusivePhotonCrossSection", _crossSectionEdges);
      book(_isolated, "IsolatedPhotonCrossSection", _crossSectionEdges);
      book(_isolatedYield, "IsolatedPhotonYield", _asymmetryEdges);
    }

    void analyze(const Event& event) {
      const Particles photons = apply<PromptFinalState>(event, "PromptPhotons").particles();
      const Particles finalState = apply<FinalState>(event, "FS").particles();
      // An event can supply several photons to a bin. Accumulate its full
      // bin contribution before filling: sumW2 must contain (w * sum f)^2,
      // not sum (w*f)^2 for correlated photons in the same collision.
      // Rivet applies the event weight once, so do not multiply by it here.
      std::vector<double> inclusive(_crossSectionEdges.size(), 0.0);
      std::vector<double> isolatedSpectrum(_crossSectionEdges.size(), 0.0);
      std::vector<double> isolatedYield(_asymmetryEdges.size(), 0.0);
      for (const Particle& photon : photons) {
        const double pt = photon.pT()/GeV;
        const size_t spectralBin = binIndex(_crossSectionEdges, pt);
        inclusive[spectralBin] += 1.0/pt;
        if (!isolated(photon, finalState)) continue;
        isolatedSpectrum[spectralBin] += 1.0/pt;
        isolatedYield[binIndex(_asymmetryEdges, pt)] += 1.0;
      }
      fillEventBins(_inclusive, _crossSectionEdges, inclusive);
      fillEventBins(_isolated, _crossSectionEdges, isolatedSpectrum);
      fillEventBins(_isolatedYield, _asymmetryEdges, isolatedYield);
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
    // PromptPhotons already enforces the common lower edge, pT > 6 GeV.
    // The final vector element is the overflow (needed for ALL at pT >=20).
    static size_t binIndex(const std::vector<double>& edges, double pt) {
      return std::upper_bound(edges.begin(), edges.end(), pt)-edges.begin()-1;
    }

    static void fillEventBins(Histo1DPtr histogram,
                              const std::vector<double>& edges,
                              const std::vector<double>& contributions) {
      for (size_t i=0; i<contributions.size(); ++i) {
        if (contributions[i] == 0.0) continue;
        // Use the bin center only as its address; the actual photon 1/pT
        // factors have already been summed, without a bin-center approximation.
        const double x = i+1<edges.size() ? 0.5*(edges[i]+edges[i+1]) : edges.back();
        histogram->fill(x, contributions[i]);
      }
    }

    static bool isolated(const Particle& photon, const Particles& finalState) {
      double coneEnergy = 0.0;
      for (const Particle& particle : finalState) {
        if (particle.genParticle() == photon.genParticle()) continue;
        if (PID::isNeutrino(particle.pid()) || deltaR(photon, particle) >= 0.5) continue;
        if (PID::isCharged(particle.pid())) {
          if (particle.p3().mod()/GeV < 0.2) continue;
        } else if (particle.E()/GeV < 0.3) continue;
        coneEnergy += particle.E();
      }
      return coneEnergy < 0.10*photon.E();
    }

    std::vector<double> _crossSectionEdges, _asymmetryEdges;
    Histo1DPtr _inclusive, _isolated, _isolatedYield;
  };

  RIVET_DECLARE_PLUGIN(PHENIX_2023_I2033856);
}
