// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Tools/Beams.hh"
#include <algorithm>
#include <cmath>
#include <vector>

namespace Rivet {

  /// HERMES Table 7 Born A_parallel in 19 x slices and 45 (x,Q2) bins.
  class HERMES_2007_I726689_LEGACY : public Analysis {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(HERMES_2007_I726689_LEGACY);

    struct Kinematics {
      bool valid = false;
      double Q2 = -1.0, x = -1.0, y = -1.0, W2 = -1.0, theta = -1.0;
    };

    void init() {
      declare(PromptFinalState(Cuts::pid == -11), "PromptPositrons");
      _xEdges = {0.0041, 0.0073, 0.0118, 0.0168, 0.0212, 0.0295,
                 0.0362, 0.0444, 0.0568, 0.0727, 0.0929, 0.119,
                 0.152, 0.194, 0.249, 0.318, 0.406, 0.520, 0.665, 0.9};
      _q2Edges = {
        {0.18,20.0},{0.18,20.0},{0.18,20.0},{0.18,20.0},
        {0.18,1.0,20.0},{0.18,1.0,20.0},{0.18,1.0,20.0},{0.18,1.0,20.0},
        {0.18,1.505,2.27,20.0},{0.18,1.62,2.62,20.0},
        {0.18,1.74,3.02,20.0},{0.18,1.88,3.49,20.0},
        {0.18,2.06,4.02,20.0},{0.18,2.23,4.61,20.0},
        {0.18,2.66,5.5,20.0},{0.18,3.3,6.625,20.0},
        {0.18,4.09,7.98,20.0},{0.18,5.04,9.455,20.0},
        {0.18,7.645,12.5,20.0}
      };
      for (size_t index = 0; index < _q2Edges.size(); ++index) {
        const std::string token = index+1 < 10 ? "0"+toString(index+1) : toString(index+1);
        Slice histograms;
        book(histograms.ordinary, "SigmaQ2_X"+token, _q2Edges[index]);
        book(histograms.duplicate, "SigmaDuplicateQ2_X"+token, _q2Edges[index]);
        book(histograms.covariance, "CovarianceProxyQ2_X"+token, _q2Edges[index]);
        _scaled.insert(_scaled.end(), {histograms.ordinary, histograms.duplicate,
                                       histograms.covariance});
        _slices.push_back(histograms);
      }
      book(_acceptedX, "Accepted_X", _xEdges);
      book(_acceptedQ2, "Accepted_Q2", logspace(50, 0.18, 20.0));
      book(_acceptedY, "Accepted_Y", 32, 0.10, 0.91);
      book(_acceptedW2, "Accepted_W2", 40, 3.24, 60.0);
      book(_acceptedTheta, "Accepted_Theta", 36, 0.04, 0.22);
      book(_acceptedLowQ2, "Accepted_Q2_Below1", 32, 0.18, 1.0);
      _scaled.insert(_scaled.end(), {_acceptedX, _acceptedQ2, _acceptedY,
                                     _acceptedW2, _acceptedTheta, _acceptedLowQ2});
    }

    void analyze(const Event& event) {
      const Kinematics dis = reconstruct(event);
      if (!dis.valid) vetoEvent;
      if (dis.Q2 <= 0.18 || dis.Q2 >= 20.0) vetoEvent;
      if (dis.y <= 0.10 || dis.y >= 0.91) vetoEvent;
      if (dis.W2 <= 3.24) vetoEvent;
      if (dis.theta <= 0.04 || dis.theta >= 0.22) vetoEvent;
      if (dis.x < _xEdges.front() || dis.x >= _xEdges.back()) vetoEvent;
      const size_t slice = std::upper_bound(_xEdges.begin(), _xEdges.end(), dis.x)
                         - _xEdges.begin() - 1;
      if (slice >= _slices.size()) vetoEvent;
      _slices[slice].ordinary->fill(dis.Q2);
      // A_parallel has no depolarization reweighting.  The duplicate and
      // covariance proxy keep the generic normalized-bin combiner applicable.
      _slices[slice].duplicate->fill(dis.Q2);
      _slices[slice].covariance->fill(dis.Q2);
      _acceptedX->fill(dis.x); _acceptedQ2->fill(dis.Q2);
      _acceptedY->fill(dis.y); _acceptedW2->fill(dis.W2);
      _acceptedTheta->fill(dis.theta);
      if (dis.Q2 < 1.0) _acceptedLowQ2->fill(dis.Q2);
    }

    void finalize() {
      if (sumW() == 0.0) return;
      const double sf = crossSection()/picobarn/sumW();
      for (Histo1DPtr& hist : _scaled) scale(hist, sf);
    }

  private:
    struct Slice { Histo1DPtr ordinary, duplicate, covariance; };

    Kinematics reconstruct(const Event& event) const {
      Kinematics out;
      const ParticlePair incoming = Rivet::beams(event);
      Particle lepton, hadron;
      if (PID::isLepton(incoming.first.pid()) && PID::isHadron(incoming.second.pid())) {
        lepton = incoming.first; hadron = incoming.second;
      } else if (PID::isLepton(incoming.second.pid()) && PID::isHadron(incoming.first.pid())) {
        lepton = incoming.second; hadron = incoming.first;
      } else return out;
      if (lepton.pid() != -11 || hadron.pid() != 2212) return out;
      Particles candidates;
      for (const Particle& particle :
           apply<PromptFinalState>(event, "PromptPositrons").particles())
        if (particle.pid() == lepton.pid()) candidates.push_back(particle);
      if (candidates.empty()) return out;
      const Particle outgoing = *std::max_element(candidates.begin(), candidates.end(),
        [](const Particle& a, const Particle& b) { return a.E() < b.E(); });
      const FourMomentum k = lepton.momentum(), kp = outgoing.momentum();
      const FourMomentum P = hadron.momentum(), q = k-kp;
      const double pdotq = P*q, pdotk = P*k;
      if (pdotq <= 0.0 || pdotk <= 0.0) return out;
      out.Q2 = -q.mass2()/GeV2; out.x = (-q.mass2())/(2.0*pdotq);
      out.y = pdotq/pdotk; out.W2 = (P+q).mass2()/GeV2;
      out.theta = k.angle(kp);
      out.valid = std::isfinite(out.Q2) && std::isfinite(out.x) &&
                  std::isfinite(out.y) && std::isfinite(out.W2) &&
                  std::isfinite(out.theta);
      return out;
    }

    std::vector<double> _xEdges;
    std::vector<std::vector<double>> _q2Edges;
    std::vector<Slice> _slices;
    std::vector<Histo1DPtr> _scaled;
    Histo1DPtr _acceptedX, _acceptedQ2, _acceptedY, _acceptedW2;
    Histo1DPtr _acceptedTheta, _acceptedLowQ2;
  };

  RIVET_DECLARE_PLUGIN(HERMES_2007_I726689_LEGACY);
}
