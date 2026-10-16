// -*- C++ -*-
#include "COMPASSInclusiveDIS.hh"
#include "E142DIS.hh"

namespace Rivet {

  /// Inclusive electron DIS inputs for E142 A1(helium-3) and A1(neutron).
  /// The campaign combines independently normalized p/n and helicity samples.
  class E142_1996_I424108 : public Analysis {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(E142_1996_I424108);

    void init() {
      declare(PromptFinalState(Cuts::pid == 11), "PromptElectrons");
      const std::vector<double> edges(E142DIS::xEdges.begin(), E142DIS::xEdges.end());
      for (size_t i = 0; i < _windows.size(); ++i) {
        Window& w = _windows[i];
        const std::string suffix = E142DIS::windowNames[i];
        book(w.ordinary, "SigmaX_"+suffix, edges);
        book(w.weighted, "SigmaOverD_X_"+suffix, edges);
        book(w.covariance, "CovarianceProxy_X_"+suffix, edges);
        book(w.x, "Accepted_X_"+suffix, edges);
        book(w.q2, "Accepted_Q2_"+suffix, 51, 1.0, 6.1);
        book(w.y, "Accepted_Y_"+suffix, 50, 0.0, 1.0);
        book(w.w2, "Accepted_W2_"+suffix, 50, 4.0, 44.0);
        book(w.theta, "Accepted_Theta_"+suffix, 50, 0.0, 0.25);
        book(w.eprime, "Accepted_EPrime_"+suffix, 52, 0.0, 26.0);
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
        event, apply<PromptFinalState>(event, "PromptElectrons").particles(), {11});
      if (!dis.valid || E142DIS::xBin(dis.x) < 0) vetoEvent;
      if (dis.Q2 <= 1.0 || dis.W2 <= 4.0 || dis.y <= 0.0 || dis.y >= 1.0) vetoEvent;

      const double mass = dis.target.mass()/GeV;
      const double D = E142DIS::depolarization(dis.x, dis.y, dis.Q2, mass);
      if (!std::isfinite(D) || D <= 0.0) vetoEvent;
      for (size_t i = 0; i < _windows.size(); ++i) {
        if (!E142DIS::inWindow(dis.x, dis.Q2, E142DIS::windowFractions[i])) continue;
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

  RIVET_DECLARE_PLUGIN(E142_1996_I424108);
}
