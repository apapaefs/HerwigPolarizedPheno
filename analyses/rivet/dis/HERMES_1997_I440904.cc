// -*- C++ -*-
#include "COMPASSInclusiveDIS.hh"
#include "HERMES1997DIS.hh"

namespace Rivet {

  /// Inclusive positron DIS inputs for HERMES1997 A1(neutron).
  /// The campaign combines independently normalized neutron helicity samples.
  class HERMES_1997_I440904 : public Analysis {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(HERMES_1997_I440904);

    void init() {
      declare(PromptFinalState(Cuts::pid == -11), "PromptPositrons");
      const std::vector<double> edges(HERMES1997DIS::xEdges.begin(), HERMES1997DIS::xEdges.end());
      for (size_t i = 0; i < _windows.size(); ++i) {
        Window& w = _windows[i];
        const std::string suffix = HERMES1997DIS::windowNames[i];
        book(w.ordinary, "SigmaX_"+suffix, edges);
        book(w.weighted, "SigmaOverD_X_"+suffix, edges);
        book(w.covariance, "CovarianceProxy_X_"+suffix, edges);
        book(w.x, "Accepted_X_"+suffix, edges);
        book(w.q2, "Accepted_Q2_"+suffix, 56, 1.0, 15.0);
        book(w.y, "Accepted_Y_"+suffix, 50, 0.0, 1.0);
        book(w.w2, "Accepted_W2_"+suffix, 92, 4.0, 50.0);
        book(w.theta, "Accepted_Theta_"+suffix, 50, 0.0, 0.25);
        book(w.eprime, "Accepted_EPrime_"+suffix, 55, 0.0, 27.5);
        // First moments per x bin expose the actual cross-section-weighted
        // kinematics, which need not equal the experimental published means.
        book(w.xMoment, "XMoment_"+suffix, edges);
        book(w.q2Moment, "Q2Moment_"+suffix, edges);
        for (const Histo1DPtr& h : {w.ordinary, w.weighted, w.covariance,
             w.x, w.q2, w.y, w.w2, w.theta, w.eprime, w.xMoment, w.q2Moment})
          _scaled.push_back(h);
      }
    }

    void analyze(const Event& event) {
      const auto dis = COMPASSInclusiveDIS::reconstructForLeptons(
        event, apply<PromptFinalState>(event, "PromptPositrons").particles(), {-11});
      if (!dis.valid || dis.targetPid != 2112 || HERMES1997DIS::xBin(dis.x) < 0) vetoEvent;
      if (!HERMES1997DIS::passesFiducial(dis.Q2, dis.W2, dis.y, dis.theta)) vetoEvent;

      const double mass = dis.target.mass()/GeV;
      const double D = HERMES1997DIS::depolarization(dis.x, dis.y, dis.Q2, mass);
      if (!std::isfinite(D) || D <= 0.0) vetoEvent;
      for (size_t i = 0; i < _windows.size(); ++i) {
        if (!HERMES1997DIS::inWindow(dis.x, dis.Q2, HERMES1997DIS::windowFractions[i])) continue;
        Window& w = _windows[i];
        COMPASSInclusiveDIS::fillMeasurement(dis, D, w.ordinary, w.weighted, w.covariance);
        w.x->fill(dis.x);
        w.q2->fill(dis.Q2);
        w.y->fill(dis.y);
        w.w2->fill(dis.W2);
        w.theta->fill(dis.theta);
        w.eprime->fill(dis.outgoingEnergy);
        w.xMoment->fill(dis.x, dis.x);
        w.q2Moment->fill(dis.x, dis.Q2);
      }
    }

    void finalize() {
      if (sumW() == 0.0) return;
      const double sf = crossSection()/picobarn/sumW();
      for (Histo1DPtr& hist : _scaled) scale(hist, sf);
    }

  private:
    struct Window {
      Histo1DPtr ordinary, weighted, covariance;
      Histo1DPtr x, q2, y, w2, theta, eprime, xMoment, q2Moment;
    };
    std::array<Window, 3> _windows;
    std::vector<Histo1DPtr> _scaled;
  };

  RIVET_DECLARE_PLUGIN(HERMES_1997_I440904);
}
