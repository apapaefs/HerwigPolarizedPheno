// -*- C++ -*-
#include "STARPolarizedJets.hh"

namespace Rivet {

  /// STAR inclusive-jet and topology-resolved dijet A_LL at 510 GeV.
  class STAR_2022_I1949588 : public Analysis {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(STAR_2022_I1949588);

    void init() {
      declare(VisibleFinalState(), "VisibleFS");
      _hadronLevel = toUpper(getOption("LEVEL", "HARDPARTON")) == "HADRON";
      book(_inclusive, "inclusive_Yield",
           {7.0,8.2,9.6,11.2,13.1,15.3,17.9,20.9,24.5,28.7,
            33.6,39.3,46.0,53.8,62.8});
      book(_dijetA, "dijet_A_Yield",
           {12.0,14.0,17.0,20.0,24.0,29.0,34.0,41.0,49.0,
            59.0,70.0,84.0,101.0});
      book(_dijetB, "dijet_B_Yield",
           {12.0,14.0,17.0,20.0,24.0,29.0,34.0,41.0,49.0,
            59.0,70.0,84.0,101.0,121.0});
      book(_dijetC, "dijet_C_Yield",
           {12.0,14.0,17.0,20.0,24.0,29.0,34.0,41.0,49.0,
            59.0,70.0,84.0,101.0});
      book(_dijetD, "dijet_D_Yield",
           {14.0,17.0,20.0,24.0,29.0,34.0,41.0,49.0,59.0,
            70.0,84.0,101.0,121.0});
      book(_provenance, "ProvenanceStatus", {0.0,1.0,2.0});
    }

    void analyze(const Event& event) {
      using namespace STARPolarizedJets;
      bool provenance = true;
      const std::vector<ClusteredJet> jets = _hadronLevel
        ? stableParticleJets(
            apply<VisibleFinalState>(event, "VisibleFS").particles(), 0.5)
        : hardPartonJets(event, 0.5, provenance);
      _provenance->fill(provenance ? 0.5 : 1.5);
      if (!provenance) vetoEvent;

      for (const ClusteredJet& jet : jets) {
        if (std::abs(jet.eta()) < 0.9) _inclusive->fill(jet.pT()/GeV);
      }
      if (jets.size() < 2) return;
      const ClusteredJet& first = jets[0];
      const ClusteredJet& second = jets[1];
      if (first.pT()/GeV <= 7.0 || second.pT()/GeV <= 5.0) return;
      if (deltaPhi(first, second) <= 2.0*M_PI/3.0) return;
      if (std::abs(first.eta()-second.eta()) >= 1.6) return;

      const double absFirst = std::abs(first.eta());
      const double absSecond = std::abs(second.eta());
      const bool firstForward = absFirst > 0.3 && absFirst < 0.9;
      const bool secondForward = absSecond > 0.3 && absSecond < 0.9;
      const bool firstCentral = absFirst < 0.3;
      const bool secondCentral = absSecond < 0.3;
      const bool sameSign = first.eta()*second.eta() > 0.0;
      const double mass = (first.momentum+second.momentum).mass()/GeV;

      if (firstForward && secondForward && sameSign) _dijetA->fill(mass);
      else if ((firstForward && secondCentral) ||
               (secondForward && firstCentral)) _dijetB->fill(mass);
      else if (firstCentral && secondCentral) _dijetC->fill(mass);
      else if (firstForward && secondForward && !sameSign) _dijetD->fill(mass);
    }

    void finalize() {
      if (sumW() == 0.0) return;
      const double factor = crossSection()/picobarn/sumW();
      scale(_inclusive, factor);
      scale(_dijetA, factor);
      scale(_dijetB, factor);
      scale(_dijetC, factor);
      scale(_dijetD, factor);
    }

  private:
    bool _hadronLevel = false;
    Histo1DPtr _inclusive, _dijetA, _dijetB, _dijetC, _dijetD;
    Histo1DPtr _provenance;
  };

  RIVET_DECLARE_PLUGIN(STAR_2022_I1949588);
}
