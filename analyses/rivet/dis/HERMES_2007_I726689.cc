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
      double thetaX = -1.0;
      double thetaY = -1.0;
    };

    struct HistogramSet {
      std::vector<Histo1DPtr> _hBornSigmaQ2;
      Histo1DPtr _hSigmaQ2ProjectionGT1, _hSigmaQ2ProjectionGT4;
      Histo1DPtr _hSigmaQ2GT1, _hSigmaOverDQ2GT1, _hCovarianceQ2GT1;
      Histo1DPtr _hSigmaQ2GT4, _hSigmaOverDQ2GT4, _hCovarianceQ2GT4;
      Histo1DPtr _hAcceptedXQ2GT1, _hAcceptedQ2Q2GT1, _hAcceptedYQ2GT1;
      Histo1DPtr _hAcceptedW2Q2GT1, _hAcceptedThetaQ2GT1;
      Histo1DPtr _hAcceptedXQ2GT4, _hAcceptedQ2Q2GT4, _hAcceptedYQ2GT4;
      Histo1DPtr _hAcceptedW2Q2GT4, _hAcceptedThetaQ2GT4;
      Histo1DPtr _hAcceptedW2FineQ2GT1;
    };

    void init() {
      declare(PromptFinalState(Cuts::pid == -11), "PromptPositrons");

      _xEdges = {0.0212, 0.0295, 0.0362, 0.0444, 0.0568, 0.0727,
                 0.0929, 0.119, 0.152, 0.194, 0.249, 0.318, 0.406,
                 0.520, 0.665, 0.900};

      // Born A_parallel cells use the analysis binning in Ehrenfried's
      // Appendix C, restricted to the final publication's 0.18<Q2<20 window.
      // Keep these physical x boundaries separate from rounded A1 x bins.
      _bornXEdges = {0.00406, 0.00733, 0.01177, 0.01677, 0.02124,
                    0.02948, 0.03619, 0.04442, 0.05681, 0.07265,
                    0.09291, 0.11882, 0.15195, 0.19432, 0.24850,
                    0.31780, 0.40641, 0.51974, 0.66466, 0.90000};
      _bornQ2Edges = {
        {0.18, 1.0}, {0.18, 1.0}, {0.18, 1.0}, {0.18, 1.0},
        {0.18, 1.0, 20.0}, {0.18, 1.0, 20.0},
        {0.18, 1.0, 20.0}, {0.18, 1.0, 20.0},
        {1.0, 1.505, 2.265, 20.0}, {1.0, 1.620, 2.623, 20.0},
        {1.0, 1.740, 3.026, 20.0}, {1.0, 1.882, 3.491, 20.0},
        {1.0, 2.061, 4.032, 20.0}, {1.0, 2.237, 4.614, 20.0},
        {1.0, 2.657, 5.491, 20.0}, {1.0, 3.305, 6.645, 20.0},
        {1.0, 4.093, 7.967, 20.0}, {1.0, 5.043, 9.458, 20.0},
        {1.0, 7.655, 12.528, 20.0}
      };
      bookSelection(_primary, "");
      bookSelection(_ringControl, "RingControl_");

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
      // Both selections see the same event and normalization. Every primary
      // contribution is also in the ring control; its squared-weight moments
      // therefore supply the shared-event covariance in postprocessing.
      fillSelection(dis, D, _ringControl);
      if (inRectangularAperture(dis.thetaX, dis.thetaY)) fillSelection(dis, D, _primary);
    }

    void finalize() {
      if (sumW() == 0.0) return;
      const double sf = crossSection() / picobarn / sumW();
      for (Histo1DPtr& hist : _scaled) scale(hist, sf);
    }

  private:

    void bookSelection(HistogramSet& hist, const std::string& prefix) {
      for (size_t index = 0; index < _bornQ2Edges.size(); ++index) {
        const std::string token = index+1 < 10 ? "0"+toString(index+1) : toString(index+1);
        Histo1DPtr histogram;
        book(histogram, prefix+"BornSigmaQ2_X"+token, _bornQ2Edges[index]);
        hist._hBornSigmaQ2.push_back(histogram);
        _scaled.push_back(histogram);
      }
      // The Q2 projections integrate the same rounded x range as the x
      // projections. They are MC ratios of integrated cross sections.
      book(hist._hSigmaQ2ProjectionGT1, prefix+"SigmaQ2_Q2GT1",
           std::vector<double>{1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 20.0});
      book(hist._hSigmaQ2ProjectionGT4, prefix+"SigmaQ2_Q2GT4",
           std::vector<double>{4.0, 6.0, 8.0, 12.0, 20.0});

      book(hist._hSigmaQ2GT1, prefix+"SigmaX_Q2GT1", _xEdges);
      book(hist._hSigmaOverDQ2GT1, prefix+"SigmaOverD_X_Q2GT1", _xEdges);
      book(hist._hCovarianceQ2GT1, prefix+"CovarianceProxy_X_Q2GT1", _xEdges);
      book(hist._hSigmaQ2GT4, prefix+"SigmaX_Q2GT4", _xEdges);
      book(hist._hSigmaOverDQ2GT4, prefix+"SigmaOverD_X_Q2GT4", _xEdges);
      book(hist._hCovarianceQ2GT4, prefix+"CovarianceProxy_X_Q2GT4", _xEdges);

      book(hist._hAcceptedXQ2GT1, prefix+"Accepted_X_Q2GT1", _xEdges);
      book(hist._hAcceptedQ2Q2GT1, prefix+"Accepted_Q2_Q2GT1", logspace(30, 1.0, 20.0));
      book(hist._hAcceptedYQ2GT1, prefix+"Accepted_Y_Q2GT1", 27, 0.10, 0.91);
      book(hist._hAcceptedW2Q2GT1, prefix+"Accepted_W2_Q2GT1", 30, 3.24, 60.0);
      book(hist._hAcceptedThetaQ2GT1, prefix+"Accepted_Theta_Q2GT1", 18, 0.04, 0.22);
      book(hist._hAcceptedXQ2GT4, prefix+"Accepted_X_Q2GT4", _xEdges);
      book(hist._hAcceptedQ2Q2GT4, prefix+"Accepted_Q2_Q2GT4", logspace(20, 4.0, 20.0));
      book(hist._hAcceptedYQ2GT4, prefix+"Accepted_Y_Q2GT4", 27, 0.10, 0.91);
      book(hist._hAcceptedW2Q2GT4, prefix+"Accepted_W2_Q2GT4", 30, 3.24, 60.0);
      book(hist._hAcceptedThetaQ2GT4, prefix+"Accepted_Theta_Q2GT4", 18, 0.04, 0.22);

      book(hist._hAcceptedW2FineQ2GT1, prefix+"Accepted_W2Fine_Q2GT1", 28, 3.24, 6.04);

      _scaled.insert(_scaled.end(), {
        hist._hAcceptedW2FineQ2GT1,
        hist._hSigmaQ2ProjectionGT1, hist._hSigmaQ2ProjectionGT4,
        hist._hSigmaQ2GT1, hist._hSigmaOverDQ2GT1, hist._hCovarianceQ2GT1,
        hist._hSigmaQ2GT4, hist._hSigmaOverDQ2GT4, hist._hCovarianceQ2GT4,
        hist._hAcceptedXQ2GT1, hist._hAcceptedQ2Q2GT1, hist._hAcceptedYQ2GT1,
        hist._hAcceptedW2Q2GT1, hist._hAcceptedThetaQ2GT1,
        hist._hAcceptedXQ2GT4, hist._hAcceptedQ2Q2GT4, hist._hAcceptedYQ2GT4,
        hist._hAcceptedW2Q2GT4, hist._hAcceptedThetaQ2GT4
      });

    }

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
      // Project the laboratory scattering direction onto orthogonal planes
      // containing the incoming beam. Longitudinal inclusive DIS is azimuthally
      // symmetric, so a rotation of which transverse axis is called vertical
      // does not change the predicted aperture average. The polar-angle veto
      // is retained separately from this rectangular spectrometer aperture.
      const Vector3 beam = k.p3().unit();
      const Vector3 reference = std::abs(beam.x()) < 0.9 ? Vector3(1,0,0) : Vector3(0,1,0);
      const Vector3 horizontal = (reference - beam*reference.dot(beam)).unit();
      const Vector3 vertical = beam.cross(horizontal);
      const double longitudinal = kp.p3().dot(beam);
      out.thetaX = std::atan2(kp.p3().dot(horizontal), longitudinal);
      out.thetaY = std::atan2(kp.p3().dot(vertical), longitudinal);
      out.valid = std::isfinite(out.Q2) && std::isfinite(out.x) &&
                  std::isfinite(out.y) && std::isfinite(out.W2) &&
                  std::isfinite(out.theta) && std::isfinite(out.thetaX) &&
                  std::isfinite(out.thetaY);
      return out;
    }

    static bool inRectangularAperture(double thetaX, double thetaY) {
      return std::abs(thetaX) < 0.17 && std::abs(thetaY) > 0.04 &&
             std::abs(thetaY) < 0.14;
    }

    void fillSelection(const DISKinematicsView& dis, double D, HistogramSet& hist) {
      hist._hAcceptedW2FineQ2GT1->fill(dis.W2);
      // Direct Born A_parallel needs no empirical D or R conversion.
      fillBornCell(dis, hist._hBornSigmaQ2);
      hist._hSigmaQ2ProjectionGT1->fill(dis.Q2);
      if (dis.Q2 > 4.0) hist._hSigmaQ2ProjectionGT4->fill(dis.Q2);

      fillMeasurement(dis, D, hist._hSigmaQ2GT1, hist._hSigmaOverDQ2GT1,
                      hist._hCovarianceQ2GT1);
      fillDiagnostics(dis, hist._hAcceptedXQ2GT1, hist._hAcceptedQ2Q2GT1,
                      hist._hAcceptedYQ2GT1, hist._hAcceptedW2Q2GT1,
                      hist._hAcceptedThetaQ2GT1);

      if (dis.Q2 > 4.0) {
        fillMeasurement(dis, D, hist._hSigmaQ2GT4, hist._hSigmaOverDQ2GT4,
                        hist._hCovarianceQ2GT4);
        fillDiagnostics(dis, hist._hAcceptedXQ2GT4, hist._hAcceptedQ2Q2GT4,
                        hist._hAcceptedYQ2GT4, hist._hAcceptedW2Q2GT4,
                        hist._hAcceptedThetaQ2GT4);
      }
    }

    void fillBornCell(const DISKinematicsView& dis,
                      const std::vector<Histo1DPtr>& bornHistograms) {
      if (dis.x < _bornXEdges.front() || dis.x >= _bornXEdges.back()) return;
      const size_t slice = std::upper_bound(_bornXEdges.begin(), _bornXEdges.end(), dis.x)
                         - _bornXEdges.begin() - 1;
      if (slice >= bornHistograms.size()) return;
      // The generator's Q2>=1 support excludes the eight lower cells. They
      // remain empty raw bins and are explicitly masked by the postprocessor.
      if (dis.Q2 < 1.0 || dis.Q2 < _bornQ2Edges[slice].front() ||
          dis.Q2 >= _bornQ2Edges[slice].back()) return;
      bornHistograms[slice]->fill(dis.Q2);
    }

    static double r1990(double x, double Q2) {
      if (!(x > 0.0) || !(Q2 > 0.0)) return std::numeric_limits<double>::quiet_NaN();
      const double scale = 0.125 * 0.125;
      const double theta = 1.0 + 12.0 * Q2 / (Q2 + 1.0) * scale / (scale + x*x);
      const double common = theta / std::log(Q2 / 0.04);
      // Whitlow thesis Eqs. (5.31)-(5.35), also used in the E143 R1990
      // prescription: the mean of fits A, B and C, not fit B alone.
      const double ra = 0.0672*common
        + 0.4671/std::pow(std::pow(Q2, 4) + std::pow(1.8979, 4), 0.25);
      const double rb = 0.0635*common + 0.5747/Q2 - 0.3534/(Q2*Q2 + 0.09);
      const double threshold = 5.0*std::pow(1.0-x, 5);
      const double rc = 0.0599*common
        + 0.5088/std::sqrt((Q2-threshold)*(Q2-threshold) + 2.1081*2.1081);
      return (ra + rb + rc)/3.0;
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
      // Direct A_parallel uses ordinary helicity cross sections irrespective
      // of the empirical A1 conversion. Only inverse-D estimators require D.
      if (!std::isfinite(D) || D <= 0.0) return;
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

    std::vector<double> _xEdges, _bornXEdges;
    std::vector<std::vector<double>> _bornQ2Edges;
    HistogramSet _primary, _ringControl;
    std::vector<Histo1DPtr> _scaled;
  };

  RIVET_DECLARE_PLUGIN(HERMES_2007_I726689);

}
