// -*- C++ -*-
#include "STARPolarizedJets.hh"

#include <algorithm>
#include <cmath>
#include <vector>

namespace Rivet {

  /// Loose inclusive dijet observables for polarized-shower comparisons.
  class MC_POLDIJETS : public Analysis {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(MC_POLDIJETS);

    void init() {
      declare(VisibleFinalState(), "VisibleFS");
      _hadronLevel = toUpper(getOption("LEVEL", "HARDPARTON")) == "HADRON";
      _radius = getOption<double>("R", 0.5);
      _leadingPtMin = getOption<double>("PTJ1MIN", 5.0)*GeV;
      _subleadingPtMin = getOption<double>("PTJ2MIN", 4.0)*GeV;
      _thirdPtMin = getOption<double>("PTJ3MIN", 2.0)*GeV;
      _etaMax = getOption<double>("ETAMAX", 1.5);

      book(_dijetRate, "dijet_rate_Yield", {0.0, 1.0});
      book(_jet1Pt, "jet1_pt_Yield",
           {5.0,6.0,7.0,8.0,10.0,12.0,15.0,20.0,30.0,45.0,
            70.0,110.0,170.0,255.0});
      book(_jet2Pt, "jet2_pt_Yield",
           {4.0,5.0,6.0,7.0,8.0,10.0,12.0,15.0,20.0,30.0,
            45.0,70.0,110.0,170.0});
      book(_jet3Pt, "jet3_pt_Yield",
           {2.0,3.0,4.0,5.0,6.0,8.0,10.0,15.0,20.0,30.0,
            50.0,80.0});
      book(_ptAverage, "pt_average_Yield",
           {4.5,5.5,6.5,8.0,10.0,12.0,15.0,20.0,30.0,45.0,
            70.0,110.0,170.0,255.0});
      book(_ptRatio21, "pt_ratio_21_Yield", 20, 0.0, 1.0);
      book(_ptRatio31, "pt_ratio_31_Yield", 20, 0.0, 1.0);
      book(_dijetMass, "dijet_mass_Yield",
           {8.0,10.0,12.0,15.0,20.0,25.0,30.0,40.0,55.0,
            75.0,100.0,140.0,200.0,300.0,500.0});
      book(_dijetPt, "dijet_pt_Yield",
           {0.0,1.0,2.0,3.0,4.0,5.0,7.0,10.0,15.0,20.0,
            30.0,50.0,80.0,130.0});
      book(_dijetHt, "dijet_ht_Yield",
           {9.0,11.0,13.0,16.0,20.0,25.0,32.0,45.0,65.0,
            90.0,130.0,190.0,280.0,500.0});
      book(_deltaPhi, "delta_phi_Yield", 18, 0.0, M_PI);
      book(_absDeltaEta, "abs_delta_eta_Yield", 12, 0.0, 3.0);
      book(_deltaR, "delta_r_Yield", 18, 0.0, 4.5);
      book(_etaBoost, "eta_boost_Yield", 12, -1.5, 1.5);
      book(_cos2DeltaPhi, "cos2_delta_phi_Yield", 20, -1.0, 1.0);
      book(_provenance, "ProvenanceStatus", {0.0, 1.0, 2.0});
    }

    void analyze(const Event& event) {
      using namespace STARPolarizedJets;
      bool provenance = true;
      const std::vector<ClusteredJet> allJets = _hadronLevel
        ? stableParticleJets(
            apply<VisibleFinalState>(event, "VisibleFS").particles(), _radius)
        : hardPartonJets(event, _radius, provenance);
      _provenance->fill(provenance ? 0.5 : 1.5);
      if (!provenance) vetoEvent;

      std::vector<ClusteredJet> jets;
      std::copy_if(
        allJets.begin(), allJets.end(), std::back_inserter(jets),
        [this](const ClusteredJet& jet) {
          return std::abs(jet.eta()) < _etaMax;
        });
      if (jets.size() < 2) return;

      const ClusteredJet& first = jets[0];
      const ClusteredJet& second = jets[1];
      if (first.pT() <= _leadingPtMin || second.pT() <= _subleadingPtMin) {
        return;
      }

      const double pt1 = first.pT()/GeV;
      const double pt2 = second.pT()/GeV;
      const FourMomentum pair = first.momentum + second.momentum;
      const double dphi = deltaPhi(first, second);
      const double deta = std::abs(first.eta()-second.eta());

      _dijetRate->fill(0.5);
      _jet1Pt->fill(pt1);
      _jet2Pt->fill(pt2);
      _ptAverage->fill(0.5*(pt1+pt2));
      _ptRatio21->fill(pt2/pt1);
      _dijetMass->fill(pair.mass()/GeV);
      _dijetPt->fill(pair.pT()/GeV);
      _dijetHt->fill(pt1+pt2);
      _deltaPhi->fill(dphi);
      _absDeltaEta->fill(deta);
      _deltaR->fill(std::hypot(deta, dphi));
      _etaBoost->fill(0.5*(first.eta()+second.eta()));
      _cos2DeltaPhi->fill(std::cos(2.0*dphi));

      if (jets.size() >= 3 && jets[2].pT() > _thirdPtMin) {
        const double pt3 = jets[2].pT()/GeV;
        _jet3Pt->fill(pt3);
        _ptRatio31->fill(pt3/pt1);
      }
    }

    void finalize() {
      if (sumW() == 0.0) return;
      const double factor = crossSection()/picobarn/sumW();
      for (Histo1DPtr histogram : {
             _dijetRate, _jet1Pt, _jet2Pt, _jet3Pt, _ptAverage,
             _ptRatio21, _ptRatio31, _dijetMass, _dijetPt, _dijetHt,
             _deltaPhi, _absDeltaEta, _deltaR, _etaBoost,
             _cos2DeltaPhi}) {
        scale(histogram, factor);
      }
    }

  private:
    bool _hadronLevel = false;
    double _radius = 0.5;
    double _leadingPtMin = 5.0*GeV;
    double _subleadingPtMin = 4.0*GeV;
    double _thirdPtMin = 2.0*GeV;
    double _etaMax = 1.5;

    Histo1DPtr _dijetRate, _jet1Pt, _jet2Pt, _jet3Pt, _ptAverage;
    Histo1DPtr _ptRatio21, _ptRatio31, _dijetMass, _dijetPt, _dijetHt;
    Histo1DPtr _deltaPhi, _absDeltaEta, _deltaR, _etaBoost;
    Histo1DPtr _cos2DeltaPhi, _provenance;
  };

  RIVET_DECLARE_PLUGIN(MC_POLDIJETS);
}
