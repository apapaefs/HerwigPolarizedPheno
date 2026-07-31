// -*- C++ -*-
#include "COMPASSInclusiveDIS.hh"
#include <vector>

namespace Rivet {

  /// COMPASS 2011 inclusive longitudinally polarized muon-proton DIS.
  class COMPASS_2016_I1357198 : public Analysis {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(COMPASS_2016_I1357198);

    void init() {
      declare(PromptFinalState(Cuts::pid == -13), "PromptMuons");
      _xEdges = {0.0025, 0.004, 0.005, 0.006, 0.008, 0.010, 0.014,
                 0.020, 0.030, 0.040, 0.060, 0.100, 0.150, 0.200,
                 0.250, 0.350, 0.500, 0.700};
      bookSelection("Q2GT1", 1.0);
      bookSelection("Q2GT4", 4.0);
    }

    void analyze(const Event& event) {
      const Particles finalState =
        apply<PromptFinalState>(event, "PromptMuons").particles();
      const COMPASSInclusiveDIS::Kinematics dis =
        COMPASSInclusiveDIS::reconstruct(event, finalState);
      if (!dis.valid || dis.targetPid != 2212) vetoEvent;
      if (dis.Q2 <= 1.0 || dis.Q2 > 190.0) vetoEvent;
      if (dis.y <= 0.10 || dis.y >= 0.90) vetoEvent;
      if (dis.W2 <= 12.0) vetoEvent;
      if (dis.x <= _xEdges.front() || dis.x >= _xEdges.back()) vetoEvent;
      const double D = COMPASSInclusiveDIS::depolarization(dis.x, dis.y, dis.Q2);
      if (!std::isfinite(D) || D <= 0.0) vetoEvent;
      fillSelection(dis, D, _q2gt1);
      if (dis.Q2 > 4.0) fillSelection(dis, D, _q2gt4);
    }

    void finalize() {
      if (sumW() == 0.0) return;
      const double sf = crossSection()/picobarn/sumW();
      for (Histo1DPtr& hist : _scaled) scale(hist, sf);
    }

  private:
    struct SelectionHistograms {
      Histo1DPtr ordinary, weighted, covariance;
      Histo1DPtr x, q2, y, w2, theta, outgoingEnergy;
    };

    void bookSelection(const std::string& suffix, double q2Minimum) {
      SelectionHistograms& h = suffix == "Q2GT1" ? _q2gt1 : _q2gt4;
      book(h.ordinary, "SigmaX_" + suffix, _xEdges);
      book(h.weighted, "SigmaOverD_X_" + suffix, _xEdges);
      book(h.covariance, "CovarianceProxy_X_" + suffix, _xEdges);
      book(h.x, "Accepted_X_" + suffix, _xEdges);
      book(h.q2, "Accepted_Q2_" + suffix, logspace(40, q2Minimum, 190.0));
      book(h.y, "Accepted_Y_" + suffix, 32, 0.10, 0.90);
      book(h.w2, "Accepted_W2_" + suffix, logspace(40, 12.0, 380.0));
      book(h.theta, "Accepted_Theta_" + suffix, 40, 0.0, 0.20);
      book(h.outgoingEnergy, "Accepted_EPrime_" + suffix, 40, 0.0, 200.0);
      _scaled.insert(_scaled.end(), {
        h.ordinary, h.weighted, h.covariance, h.x, h.q2, h.y,
        h.w2, h.theta, h.outgoingEnergy});
    }

    static void fillSelection(const COMPASSInclusiveDIS::Kinematics& dis,
                              double D, SelectionHistograms& h) {
      COMPASSInclusiveDIS::fillMeasurement(
        dis, D, h.ordinary, h.weighted, h.covariance);
      COMPASSInclusiveDIS::fillDiagnostics(
        dis, h.x, h.q2, h.y, h.w2, h.theta, h.outgoingEnergy);
    }

    std::vector<double> _xEdges;
    std::vector<Histo1DPtr> _scaled;
    SelectionHistograms _q2gt1, _q2gt4;
  };

  RIVET_DECLARE_PLUGIN(COMPASS_2016_I1357198);
}
