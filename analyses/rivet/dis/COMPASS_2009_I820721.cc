// -*- C++ -*-
#include "COMPASSSIDIS.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include <array>
#include <cmath>
#include <map>
#include <string>
#include <vector>

namespace Rivet {

  /// COMPASS identified-pion/kaon longitudinal SIDIS asymmetries on deuterium.
  class COMPASS_2009_I820721 : public Analysis {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(COMPASS_2009_I820721);

    struct SpeciesHistograms {
      Histo1DPtr ordinary, inverseD, covariance;
    };

    void init() {
      declare(FinalState(), "LabFS");
      declare(PromptFinalState(Cuts::pid == -13), "PromptMuons");
      _xEdges = {0.004,0.006,0.010,0.020,0.030,0.040,
                 0.060,0.100,0.150,0.200,0.300};
      for (const std::string species :
           {"piplus", "piminus", "kplus", "kminus"}) {
        SpeciesHistograms& histograms = _species[species];
        book(histograms.ordinary, "Yield_" + species + "_x", _xEdges);
        book(histograms.inverseD, "YieldOverD_" + species + "_x", _xEdges);
        book(histograms.covariance,
             "CovarianceProxy_" + species + "_x", _xEdges);
        _scaled.insert(_scaled.end(), {
          histograms.ordinary, histograms.inverseD, histograms.covariance});
      }
      book(_acceptedX, "Accepted_X", _xEdges);
      book(_acceptedQ2, "Accepted_Q2", logspace(40, 1.0, 100.0));
      book(_acceptedY, "Accepted_Y", 40, 0.1, 0.9);
      book(_acceptedW, "Accepted_W", 40, 4.0, 20.0);
      book(_acceptedZ, "Accepted_Z", 40, 0.2, 0.85);
      book(_acceptedMomentum, "Accepted_HadronMomentum", 40, 10.0, 50.0);
      book(_acceptedTheta, "Accepted_HadronTheta", 50, 0.0, 0.25);
      _scaled.insert(_scaled.end(), {
        _acceptedX, _acceptedQ2, _acceptedY, _acceptedW, _acceptedZ,
        _acceptedMomentum, _acceptedTheta});
    }

    void analyze(const Event& event) {
      const Particles prompt =
        apply<PromptFinalState>(event, "PromptMuons").particles();
      const COMPASSSIDIS::DISKinematics dis =
        COMPASSInclusiveDIS::reconstruct(event, prompt);
      if (!dis.valid || !COMPASSSIDIS::fixed160GeVBeam(dis)) vetoEvent;
      if (!(dis.Q2 > 1.0 && dis.y > 0.1 && dis.y < 0.9 &&
            dis.x > 0.004 && dis.x < 0.3)) vetoEvent;
      const double D = COMPASSInclusiveDIS::depolarization(dis.x, dis.y, dis.Q2);
      if (!(D > 0.0) || !std::isfinite(D)) vetoEvent;

      _acceptedX->fill(dis.x);
      _acceptedQ2->fill(dis.Q2);
      _acceptedY->fill(dis.y);
      if (dis.W2 > 0.0) _acceptedW->fill(std::sqrt(dis.W2));

      std::map<std::string, double> counts;
      const Particles particles = apply<FinalState>(event, "LabFS").particles();
      for (const Particle& particle : particles) {
        std::string species;
        if (particle.pid() ==  211) species = "piplus";
        if (particle.pid() == -211) species = "piminus";
        if (particle.pid() ==  321) species = "kplus";
        if (particle.pid() == -321) species = "kminus";
        if (species.empty()) continue;
        const COMPASSSIDIS::HadronKinematics hadron =
          COMPASSSIDIS::hadronKinematics(dis, particle);
        if (!hadron.valid || !(hadron.z > 0.2 && hadron.z < 0.85) ||
            !(hadron.momentum > 10.0 && hadron.momentum < 50.0)) continue;
        counts[species] += 1.0;
        _acceptedZ->fill(hadron.z);
        _acceptedMomentum->fill(hadron.momentum);
        _acceptedTheta->fill(hadron.theta);
      }

      for (const auto& entry : counts) {
        if (!(entry.second > 0.0)) continue;
        SpeciesHistograms& histograms = _species.at(entry.first);
        histograms.ordinary->fill(dis.x, entry.second);
        histograms.inverseD->fill(dis.x, entry.second/D);
        // sumW2 = n^2/D, the ordinary--inverse-D event covariance.
        histograms.covariance->fill(dis.x, entry.second/std::sqrt(D));
      }
    }

    void finalize() {
      if (sumW() == 0.0) return;
      const double factor = crossSection()/picobarn/sumW();
      for (Histo1DPtr& histogram : _scaled) scale(histogram, factor);
    }

  private:
    std::vector<double> _xEdges;
    std::map<std::string, SpeciesHistograms> _species;
    std::vector<Histo1DPtr> _scaled;
    Histo1DPtr _acceptedX, _acceptedQ2, _acceptedY, _acceptedW;
    Histo1DPtr _acceptedZ, _acceptedMomentum, _acceptedTheta;
  };

  RIVET_DECLARE_PLUGIN(COMPASS_2009_I820721);

} // namespace Rivet
