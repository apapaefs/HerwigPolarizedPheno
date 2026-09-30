// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Tools/Beams.hh"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace Rivet {

  /// HERMES 27.6 GeV DIS with independent proton/neutron target components.
  class HERMES_2007_I726689 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(HERMES_2007_I726689);

    struct DISKinematicsView {
      bool valid = false;
      double Q2 = -1.0;
      double x = -1.0;
      double y = -1.0;
      double W2 = -1.0;
      double theta = -1.0;
    };

    void init() {
      declare(PromptFinalState(Cuts::pid == -11), "PromptPositrons");

      _xEdges = {0.0212, 0.0295, 0.0362, 0.0444, 0.0568, 0.0727,
                 0.0929, 0.119, 0.152, 0.194, 0.249, 0.318, 0.406,
                 0.520, 0.665, 0.900};

      book(_hSigmaQ2GT1, "SigmaX_Q2GT1", _xEdges);
      book(_hSigmaOverDQ2GT1, "SigmaOverD_X_Q2GT1", _xEdges);
      book(_hCovarianceQ2GT1, "CovarianceProxy_X_Q2GT1", _xEdges);
      book(_hSigmaQ2GT4, "SigmaX_Q2GT4", _xEdges);
      book(_hSigmaOverDQ2GT4, "SigmaOverD_X_Q2GT4", _xEdges);
      book(_hCovarianceQ2GT4, "CovarianceProxy_X_Q2GT4", _xEdges);

      book(_hAcceptedXQ2GT1, "Accepted_X_Q2GT1", _xEdges);
      book(_hAcceptedQ2Q2GT1, "Accepted_Q2_Q2GT1", logspace(30, 1.0, 20.0));
      book(_hAcceptedYQ2GT1, "Accepted_Y_Q2GT1", 27, 0.10, 0.91);
      book(_hAcceptedW2Q2GT1, "Accepted_W2_Q2GT1", 30, 3.24, 60.0);
      book(_hAcceptedThetaQ2GT1, "Accepted_Theta_Q2GT1", 18, 0.04, 0.22);
      book(_hAcceptedXQ2GT4, "Accepted_X_Q2GT4", _xEdges);
      book(_hAcceptedQ2Q2GT4, "Accepted_Q2_Q2GT4", logspace(20, 4.0, 20.0));
      book(_hAcceptedYQ2GT4, "Accepted_Y_Q2GT4", 27, 0.10, 0.91);
      book(_hAcceptedW2Q2GT4, "Accepted_W2_Q2GT4", 30, 3.24, 60.0);
      book(_hAcceptedThetaQ2GT4, "Accepted_Theta_Q2GT4", 18, 0.04, 0.22);

      _scaled = {
        _hSigmaQ2GT1, _hSigmaOverDQ2GT1, _hCovarianceQ2GT1,
        _hSigmaQ2GT4, _hSigmaOverDQ2GT4, _hCovarianceQ2GT4,
        _hAcceptedXQ2GT1, _hAcceptedQ2Q2GT1, _hAcceptedYQ2GT1,
        _hAcceptedW2Q2GT1, _hAcceptedThetaQ2GT1,
        _hAcceptedXQ2GT4, _hAcceptedQ2Q2GT4, _hAcceptedYQ2GT4,
        _hAcceptedW2Q2GT4, _hAcceptedThetaQ2GT4
      };

      // The derived postprocessed objects carry the R1990 and neglected-A2
      // annotations. Avoid dereferencing multiplexed Rivet histogram pointers
      // here, before a concrete event-weight stream has become active.
    }

    void analyze(const Event& event) {
      const DISKinematicsView dis = disKinematics(event);
      if (!dis.valid) vetoEvent;

      if (dis.Q2 < 1.0 || dis.Q2 > 20.0) vetoEvent;
      if (dis.y <= 0.10 || dis.y > 0.91) vetoEvent;
      if (dis.W2 <= 3.24) vetoEvent;
      if (dis.theta < 0.04 || dis.theta > 0.22) vetoEvent;
      if (dis.x < _xEdges.front() || dis.x > _xEdges.back()) vetoEvent;

      const double D = depolarization(dis.x, dis.y, dis.Q2);
      if (!std::isfinite(D) || D <= 0.0) vetoEvent;

      fillMeasurement(dis, D, _hSigmaQ2GT1, _hSigmaOverDQ2GT1,
                      _hCovarianceQ2GT1);
      fillDiagnostics(dis, _hAcceptedXQ2GT1, _hAcceptedQ2Q2GT1,
                      _hAcceptedYQ2GT1, _hAcceptedW2Q2GT1,
                      _hAcceptedThetaQ2GT1);

      if (dis.Q2 > 4.0) {
        fillMeasurement(dis, D, _hSigmaQ2GT4, _hSigmaOverDQ2GT4,
                        _hCovarianceQ2GT4);
        fillDiagnostics(dis, _hAcceptedXQ2GT4, _hAcceptedQ2Q2GT4,
                        _hAcceptedYQ2GT4, _hAcceptedW2Q2GT4,
                        _hAcceptedThetaQ2GT4);
      }
    }

    void finalize() {
      if (sumW() == 0.0) return;
      const double sf = crossSection() / picobarn / sumW();
      for (Histo1DPtr& hist : _scaled) scale(hist, sf);
    }

  private:

    Particle scatteredPositron(const Event& event, const Particle& incoming) const {
      Particles candidates;
      for (const Particle& particle :
           apply<PromptFinalState>(event, "PromptPositrons").particles()) {
        if (particle.pid() == incoming.pid()) candidates.push_back(particle);
      }
      if (candidates.empty()) return Particle();
      return *std::max_element(candidates.begin(), candidates.end(),
        [](const Particle& a, const Particle& b) { return a.E() < b.E(); });
    }

    DISKinematicsView disKinematics(const Event& event) const {
      DISKinematicsView out;
      const ParticlePair incoming = Rivet::beams(event);
      Particle lepton;
      Particle hadron;
      if (PID::isLepton(incoming.first.pid()) && PID::isHadron(incoming.second.pid())) {
        lepton = incoming.first;
        hadron = incoming.second;
      } else if (PID::isLepton(incoming.second.pid()) && PID::isHadron(incoming.first.pid())) {
        lepton = incoming.second;
        hadron = incoming.first;
      } else {
        return out;
      }
      // Deuteron observables combine separately normalized proton/neutron samples
      // externally; preserve identical raw histograms and laboratory selections.
      if (lepton.pid() != -11 ||
          (hadron.pid() != 2212 && hadron.pid() != 2112)) return out;

      const Particle outgoing = scatteredPositron(event, lepton);
      if (outgoing.pid() == PID::ANY) return out;

      const FourMomentum k = lepton.momentum();
      const FourMomentum kp = outgoing.momentum();
      const FourMomentum P = hadron.momentum();
      const FourMomentum q = k - kp;
      const double Pdotq = P * q;
      const double Pdotk = P * k;
      if (Pdotq <= 0.0 || Pdotk <= 0.0) return out;

      out.Q2 = -q.mass2() / GeV2;
      out.x = (-q.mass2()) / (2.0 * Pdotq);
      out.y = Pdotq / Pdotk;
      out.W2 = (P + q).mass2() / GeV2;
      out.theta = k.angle(kp);
      out.valid = std::isfinite(out.Q2) && std::isfinite(out.x) &&
                  std::isfinite(out.y) && std::isfinite(out.W2) &&
                  std::isfinite(out.theta);
      return out;
    }

    static double r1990(double x, double Q2) {
      if (!(x > 0.0) || !(Q2 > 0.0)) return std::numeric_limits<double>::quiet_NaN();
      const double scale = 0.125 * 0.125;
      const double theta = 1.0 + 12.0 * Q2 / (Q2 + 1.0) * scale / (scale + x*x);
      return 0.0635 / std::log(Q2 / 0.04) * theta
           + 0.5747 / Q2
           - 0.3534 / (Q2*Q2 + 0.09);
    }

    static double depolarization(double x, double y, double Q2) {
      constexpr double protonMass = 0.9382720813;
      const double gamma2 = 4.0 * protonMass * protonMass * x * x / Q2;
      const double R = r1990(x, Q2);
      const double numerator = y * (2.0 - y) * (1.0 + 0.5 * gamma2 * y);
      const double denominator = y*y * (1.0 + gamma2)
        + 2.0 * (1.0 + R) * (1.0 - y - 0.25 * gamma2 * y*y);
      return numerator / denominator;
    }

    static void fillMeasurement(const DISKinematicsView& dis, double D,
                                const Histo1DPtr& ordinary,
                                const Histo1DPtr& weighted,
                                const Histo1DPtr& covariance) {
      ordinary->fill(dis.x);
      weighted->fill(dis.x, 1.0 / D);
      // sumW2 of this auxiliary fill is sum(eventWeight^2 / D), i.e. the
      // within-sample covariance of the ordinary and 1/D-weighted estimates.
      covariance->fill(dis.x, 1.0 / std::sqrt(D));
    }

    static void fillDiagnostics(const DISKinematicsView& dis,
                                const Histo1DPtr& x,
                                const Histo1DPtr& q2,
                                const Histo1DPtr& y,
                                const Histo1DPtr& w2,
                                const Histo1DPtr& theta) {
      x->fill(dis.x);
      q2->fill(dis.Q2);
      y->fill(dis.y);
      w2->fill(dis.W2);
      theta->fill(dis.theta);
    }

    std::vector<double> _xEdges;
    std::vector<Histo1DPtr> _scaled;
    Histo1DPtr _hSigmaQ2GT1, _hSigmaOverDQ2GT1, _hCovarianceQ2GT1;
    Histo1DPtr _hSigmaQ2GT4, _hSigmaOverDQ2GT4, _hCovarianceQ2GT4;
    Histo1DPtr _hAcceptedXQ2GT1, _hAcceptedQ2Q2GT1, _hAcceptedYQ2GT1;
    Histo1DPtr _hAcceptedW2Q2GT1, _hAcceptedThetaQ2GT1;
    Histo1DPtr _hAcceptedXQ2GT4, _hAcceptedQ2Q2GT4, _hAcceptedYQ2GT4;
    Histo1DPtr _hAcceptedW2Q2GT4, _hAcceptedThetaQ2GT4;
  };

  RIVET_DECLARE_PLUGIN(HERMES_2007_I726689);

}
