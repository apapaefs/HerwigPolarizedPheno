// -*- C++ -*-
#include "COMPASSSIDIS.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include <cmath>
#include <map>
#include <string>
#include <vector>

namespace Rivet {

  /// COMPASS charged-kaon SIDIS multiplicities.
  class COMPASS_2017_I1483098 : public Analysis {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(COMPASS_2017_I1483098);

    struct SpeciesHistograms {
      Histo1DPtr numerator, denominator, covariance;
    };

    void init() {
      declare(FinalState(), "LabFS");
      declare(PromptFinalState(Cuts::pid == -13), "PromptMuons");
      _cells = &COMPASSSIDIS::kaonCells();
      for (const std::string species : {"kplus", "kminus"}) {
        SpeciesHistograms& histograms = _species[species];
        const size_t count = _cells->size();
        book(histograms.numerator,
             "HadronNumerator_" + species + "_cells", count, 0.0, count);
        book(histograms.denominator,
             "DISDenominator_" + species + "_cells", count, 0.0, count);
        book(histograms.covariance,
             "CovarianceProxy_" + species + "_cells", count, 0.0, count);
        _scaled.insert(_scaled.end(), {
          histograms.numerator, histograms.denominator, histograms.covariance});
      }
      book(_acceptedX, "Accepted_X", 45, 0.004, 0.4);
      book(_acceptedQ2, "Accepted_Q2", logspace(45, 1.0, 60.0));
      book(_acceptedY, "Accepted_Y", 40, 0.1, 0.7);
      book(_acceptedW, "Accepted_W", 40, 5.0, 20.0);
      book(_acceptedZ, "Accepted_Z", 40, 0.2, 0.85);
      book(_acceptedMomentum, "Accepted_HadronMomentum", 40, 12.0, 40.0);
      book(_acceptedTheta, "Accepted_HadronTheta", 44, 0.010, 0.120);
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
      if (!(dis.Q2 > 1.0 && dis.W2 > 25.0 &&
            dis.x > 0.004 && dis.x < 0.4 &&
            dis.y > 0.1 && dis.y < 0.7)) vetoEvent;

      std::vector<size_t> denominatorCells;
      for (size_t index = 0; index < _cells->size(); ++index) {
        if (COMPASSSIDIS::matchesDIS((*_cells)[index], dis.x, dis.y))
          denominatorCells.push_back(index);
      }
      if (denominatorCells.empty()) vetoEvent;

      _acceptedX->fill(dis.x);
      _acceptedQ2->fill(dis.Q2);
      _acceptedY->fill(dis.y);
      _acceptedW->fill(std::sqrt(dis.W2));

      std::map<std::string, std::vector<double>> counts;
      counts["kplus"] = std::vector<double>(_cells->size(), 0.0);
      counts["kminus"] = std::vector<double>(_cells->size(), 0.0);
      const Particles particles = apply<FinalState>(event, "LabFS").particles();
      for (const Particle& particle : particles) {
        if (std::abs(particle.pid()) != 321) continue;
        const COMPASSSIDIS::HadronKinematics hadron =
          COMPASSSIDIS::hadronKinematics(dis, particle);
        if (!hadron.valid || !(hadron.z >= 0.2 && hadron.z <= 0.85) ||
            !(hadron.momentum > 12.0 && hadron.momentum < 40.0) ||
            !(hadron.theta > 0.010 && hadron.theta < 0.120)) continue;
        const int cell = COMPASSSIDIS::findCell(*_cells, dis.x, dis.y, hadron.z);
        if (cell < 0) continue;
        counts[particle.pid() > 0 ? "kplus" : "kminus"][cell] += 1.0;
        _acceptedZ->fill(hadron.z);
        _acceptedMomentum->fill(hadron.momentum);
        _acceptedTheta->fill(hadron.theta);
      }

      for (auto& entry : _species) {
        const std::string& species = entry.first;
        SpeciesHistograms& histograms = entry.second;
        for (const size_t index : denominatorCells)
          histograms.denominator->fill(double(index) + 0.5);
        for (size_t index = 0; index < counts[species].size(); ++index) {
          const double count = counts[species][index];
          COMPASSSIDIS::fillEventCount(
            histograms.numerator, double(index) + 0.5, count);
          COMPASSSIDIS::fillMultiplicityCovariance(
            histograms.covariance, double(index) + 0.5, count);
        }
      }
    }

    void finalize() {
      if (sumW() == 0.0) return;
      const double factor = crossSection()/picobarn/sumW();
      for (Histo1DPtr& histogram : _scaled) scale(histogram, factor);
    }

  private:
    const std::vector<COMPASSSIDIS::MultiplicityCell>* _cells = nullptr;
    std::map<std::string, SpeciesHistograms> _species;
    std::vector<Histo1DPtr> _scaled;
    Histo1DPtr _acceptedX, _acceptedQ2, _acceptedY, _acceptedW;
    Histo1DPtr _acceptedZ, _acceptedMomentum, _acceptedTheta;
  };

  RIVET_DECLARE_PLUGIN(COMPASS_2017_I1483098);

} // namespace Rivet
