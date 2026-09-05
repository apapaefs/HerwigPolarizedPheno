// -*- C++ -*-
#include "STARPolarizedJets.hh"

namespace Rivet {

  /// STAR inclusive-jet and dijet A_LL at sqrt(s)=200 GeV.
  class STAR_2021_I1850855 : public Analysis {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(STAR_2021_I1850855);

    void init() {
      declare(VisibleFinalState(), "VisibleFS");
      _hadronLevel = toUpper(getOption("LEVEL", "HARDPARTON")) == "HADRON";
      book(_inclusiveForward, "inclusive_forward_Yield",
           {6.0,7.1,8.4,9.9,11.7,13.8,16.3,19.2,22.7,26.8,31.6,37.3});
      book(_inclusiveCentral, "inclusive_central_Yield",
           {6.0,7.1,8.4,9.9,11.7,13.8,16.3,19.2,22.7,26.8,31.6,37.3});
      book(_inclusiveCombined, "inclusive_combined_Yield",
           {6.0,7.1,8.4,9.9,11.7,13.8,16.3,19.2,22.7,26.8,31.6,37.3});
      book(_dijetSame, "dijet_same_sign_Yield",
           {17.0,19.0,23.0,28.0,34.0,41.0,58.0,82.0});
      book(_dijetOpposite, "dijet_opposite_sign_Yield",
           {17.0,19.0,23.0,28.0,34.0,41.0,58.0,82.0});
      book(_provenance, "ProvenanceStatus", {0.0,1.0,2.0});
    }

    void analyze(const Event& event) {
      using namespace STARPolarizedJets;
      bool provenance = true;
      const std::vector<ClusteredJet> jets = _hadronLevel
        ? stableParticleJets(
            apply<VisibleFinalState>(event, "VisibleFS").particles(), 0.6)
        : hardPartonJets(event, 0.6, provenance);
      _provenance->fill(provenance ? 0.5 : 1.5);
      if (!provenance) vetoEvent;

      size_t selectedInclusive = 0;
      for (const ClusteredJet& jet : jets) {
        const double absEta = std::abs(jet.eta());
        const double pt = jet.pT()/GeV;
        if (absEta >= 1.0 || pt < 6.0) continue;
        if (selectedInclusive++ == 2) break;
        _inclusiveCombined->fill(pt);
        if (absEta < 0.5) _inclusiveCentral->fill(pt);
        else _inclusiveForward->fill(pt);
      }

      if (jets.size() < 2) return;
      const ClusteredJet& first = jets[0];
      const ClusteredJet& second = jets[1];
      if (first.pT()/GeV <= 8.0 || second.pT()/GeV <= 6.0) return;
      if (std::abs(first.eta()) >= 0.8 || std::abs(second.eta()) >= 0.8) return;
      if (deltaPhi(first, second) <= 2.0*M_PI/3.0) return;
      const double mass = (first.momentum+second.momentum).mass()/GeV;
      const double etaProduct = first.eta()*second.eta();
      if (etaProduct > 0.0) _dijetSame->fill(mass);
      else if (etaProduct < 0.0) _dijetOpposite->fill(mass);
    }

    void finalize() {
      if (sumW() == 0.0) return;
      const double factor = crossSection()/picobarn/sumW();
      scale(_inclusiveForward, factor);
      scale(_inclusiveCentral, factor);
      scale(_inclusiveCombined, factor);
      scale(_dijetSame, factor);
      scale(_dijetOpposite, factor);
    }

  private:
    bool _hadronLevel = false;
    Histo1DPtr _inclusiveForward, _inclusiveCentral, _inclusiveCombined;
    Histo1DPtr _dijetSame, _dijetOpposite, _provenance;
  };

  RIVET_DECLARE_PLUGIN(STAR_2021_I1850855);
}
