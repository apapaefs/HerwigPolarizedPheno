// -*- C++ -*-
// Harder-jet spin-transfer pilot. No generator history or flavour labels.
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "PolJetShapesHard.hh"
#include <map>
#include <string>

namespace Rivet {
  class MC_POLJETSHAPES_HARD : public Analysis {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(MC_POLJETSHAPES_HARD);

    void init() {
      _radius = getOption<double>("R", 0.5);
      _pt1 = getOption<double>("PTJ1MIN", 5.0);
      _pt2 = getOption<double>("PTJ2MIN", 4.0);
      _eta = getOption<double>("ETAMAX", 1.5);
      declare(FastJets(VisibleFinalState(), JetAlg::ANTIKT, _radius,
                      JetMuons::ALL, JetInvisibles::NONE), "Jets");
      book(_cutflow, "CutflowFraction", 5, 0, 5);
      for (const std::string& window : _windows) {
        for (int jet = 1; jet <= 2; ++jet) {
          for (int kt = 1; kt <= 2; ++kt) {
            const std::string suffix = "_j"+std::to_string(jet)+"_kt"+std::to_string(kt);
            bookAngle(window+"_dpsi12"+suffix);
            bookAngle(window+"_hardplane_primary"+suffix);
            bookAngle(window+"_dpsi12_j"+std::to_string(jet)+"_hardshare_kt"+std::to_string(kt));
          }
        }
        for (int kt = 1; kt <= 2; ++kt)
          bookAngle(window+"_interjet_dpsi11_kt"+std::to_string(kt));
        bookAngle(window+"_resolved_dphi31");
        bookAngle(window+"_resolved_dpsi34");
        bookSpectrum(window+"_dijet_rate", {0,1});
        for (int jet = 3; jet <= 4; ++jet)
          bookSpectrum(window+"_jet"+std::to_string(jet)+"_pt",
                       {2,3,4,5,6,8,10,12,15,20,30,45,70,110,170,255});
        for (const std::string name : {"dijet_threshold_denominator", "ge3_threshold", "ge4_threshold"})
          bookSpectrum(window+"_"+name, {1.5,2.5,3.5,4.5,5.5,7,9,11});
        for (const std::string name : {"pt31_cumulative_tail", "pt41_cumulative_tail"})
          bookSpectrum(window+"_"+name, {0.05,0.15,0.25,0.35,0.45,0.55});
      }
      book(_accepted, "AcceptedEntriesPerEvent", _angleNames.size(), 0, _angleNames.size());
    }

    void analyze(const Event& event) {
      _cutflow->fill(0.5);
      const Jets jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::abseta < _eta);
      if (jets.size() < 2) vetoEvent;
      _cutflow->fill(1.5);
      if (jets[0].pT() <= _pt1*GeV || jets[1].pT() <= _pt2*GeV) vetoEvent;
      _cutflow->fill(2.5);
      // Individual-jet windows for intra-jet angles; leading-jet window for
      // event-level observables. Keep the original loose dijet denominator.
      const int leadingWindow = PolJetShapesHard::ptWindow(jets[0].pT()/GeV);
      const int secondWindow = PolJetShapesHard::ptWindow(jets[1].pT()/GeV);
      if (leadingWindow < 0 && secondWindow < 0) return;
      std::array<PolJetShapesHard::Angles, 2> result;
      for (size_t j = 0; j < 2; ++j) {
        const int window = PolJetShapesHard::ptWindow(jets[j].pT()/GeV);
        // Jet 2 is also needed for inter-jet angles when only jet 1 is in a window.
        result[j] = PolJetShapesHard::angles(jets[j], _radius);
        if (window < 0) continue;
        for (size_t k = 0; k < 2; ++k) {
          const std::string suffix = "_j"+std::to_string(j+1)+"_kt"+std::to_string(k+1);
          if (result[j].primary[k])
            fillAngle(_windows[window]+"_hardplane_primary"+suffix, result[j].primaryPsi[k]);
          if (result[j].pair[k])
            fillAngle(_windows[window]+"_dpsi12"+suffix, result[j].pairPsi[k]);
          if (result[j].hardshare[k])
            fillAngle(_windows[window]+"_dpsi12_j"+std::to_string(j+1)+"_hardshare_kt"+std::to_string(k+1),
                      result[j].hardsharePsi[k]);
        }
      }
      if (leadingWindow < 0) return;
      _cutflow->fill(3.5+leadingWindow);
      const std::string& window = _windows[leadingWindow];
      fill(window+"_dijet_rate", 0.5);
      const double pt1 = jets[0].pT()/GeV;
      const double pt3 = jets.size() > 2 ? jets[2].pT()/GeV : 0;
      const double pt4 = jets.size() > 3 ? jets[3].pT()/GeV : 0;
      if (pt3 > 2) fill(window+"_jet3_pt", pt3);
      if (pt4 > 2) fill(window+"_jet4_pt", pt4);
      for (double threshold : {2.,3.,4.,5.,6.,8.,10.}) {
        fill(window+"_dijet_threshold_denominator", threshold);
        if (pt3 > threshold) fill(window+"_ge3_threshold", threshold);
        if (pt4 > threshold) fill(window+"_ge4_threshold", threshold);
      }
      for (double threshold : {0.1,0.2,0.3,0.4,0.5}) {
        if (pt3/pt1 > threshold) fill(window+"_pt31_cumulative_tail", threshold);
        if (pt4/pt1 > threshold) fill(window+"_pt41_cumulative_tail", threshold);
      }
      for (size_t k = 0; k < 2; ++k)
        if (result[0].primary[k] && result[1].primary[k])
          fillAngle(window+"_interjet_dpsi11_kt"+std::to_string(k+1),
                    result[0].primaryPsi[k]-result[1].primaryPsi[k]);
      Vector3 beamPlane, thirdPlane, fourthPlane;
      double angle = 0;
      if (pt3 > 2 && PolJetShapesHard::planeNormal(jets[0].p3(), jets[2].p3(), thirdPlane)) {
        if (PolJetShapesHard::planeNormal(Vector3(0,0,1), jets[0].p3(), beamPlane) &&
            PolJetShapesHard::signedPlaneAngle(beamPlane, thirdPlane, jets[0].p3(), angle))
          fillAngle(window+"_resolved_dphi31", angle);
        if (pt4 > 2 && PolJetShapesHard::planeNormal(jets[0].p3(), jets[3].p3(), fourthPlane) &&
            PolJetShapesHard::signedPlaneAngle(thirdPlane, fourthPlane, jets[0].p3(), angle))
          fillAngle(window+"_resolved_dpsi34", angle);
      }
    }

    void finalize() {
      if (sumW() == 0) return;
      // Rivet applies the compensated Herwig event weight automatically.
      // Counts here are weighted fractions, NOT literal sampled-event counts.
      for (auto& item : _histograms)
        scale(item.second, crossSection()/picobarn/sumW());
      scale(_cutflow, 1.0/sumW());
      scale(_accepted, 1.0/sumW());
    }

  private:
    void bookSpectrum(const std::string& name, const std::vector<double>& edges) {
      book(_histograms[name], name+"_Yield", edges);
    }
    void bookAngle(const std::string& name) {
      book(_histograms[name], name+"_Yield", PolJetShapesHard::NANGLE, -M_PI, M_PI);
      _angleNames.push_back(name);
    }
    void fill(const std::string& name, double value) {
      _histograms.at(name)->fill(value);
    }
    void fillAngle(const std::string& name, double value) {
      fill(name, PolJetShapesHard::wrapToPi(value));
      const auto found = std::find(_angleNames.begin(), _angleNames.end(), name);
      _accepted->fill(std::distance(_angleNames.begin(), found)+0.5);
    }
    double _radius, _pt1, _pt2, _eta;
    const std::array<std::string, 2> _windows{{"pt20_30", "pt30_45"}};
    std::map<std::string, Histo1DPtr> _histograms;
    std::vector<std::string> _angleNames;
    Histo1DPtr _cutflow, _accepted;
  };
  RIVET_DECLARE_PLUGIN(MC_POLJETSHAPES_HARD);
}
