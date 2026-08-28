// -*- C++ -*-
#pragma once

#include "COMPASSSIDIS.hh"
#include "SIDISTrancheBinning.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include <cmath>
#include <map>
#include <string>
#include <vector>

namespace Rivet {

  /// Shared raw-yield implementation for the 2025 proton and corrected 2026
  /// isoscalar COMPASS multiplicity releases.  Published target combinations,
  /// density widths, and signed-NLO arithmetic are deliberately postprocessed.
  template <int ReleaseYear>
  class COMPASSModernMultiplicity : public Analysis {
  public:
    explicit COMPASSModernMultiplicity(const std::string& name)
      : Analysis(name) {}

    struct Histograms {
      Histo1DPtr numerator, denominator, covariance;
    };

    void init() override {
      declare(FinalState(), "LabFS");
      declare(PromptFinalState(Cuts::pid == -13), "PromptMuons");
      configureCells();
      for (const auto& entry : _cells) {
        const std::string& species = entry.first;
        const size_t count = entry.second->size();
        Histograms& histograms = _species[species];
        book(histograms.numerator,
             "HadronNumerator_" + species + "_cells", count, 0., count);
        book(histograms.denominator,
             "DISDenominator_" + species + "_cells", count, 0., count);
        book(histograms.covariance,
             "CovarianceProxy_" + species + "_cells", count, 0., count);
        _scaled.insert(_scaled.end(), {histograms.numerator,
          histograms.denominator, histograms.covariance});
      }
      book(_acceptedX, "Accepted_X", 45, .004, .4);
      book(_acceptedQ2, "Accepted_Q2", logspace(45, 1., 60.));
      book(_acceptedY, "Accepted_Y", 40, .1, .7);
      book(_acceptedW, "Accepted_W", 40, 5., 20.);
      book(_acceptedZ, "Accepted_Z", 40, .2, .85);
      book(_acceptedMomentum, "Accepted_HadronMomentum", 40, 12., 40.);
      book(_acceptedTheta, "Accepted_HadronTheta", 44, .010, .120);
      _scaled.insert(_scaled.end(), {_acceptedX, _acceptedQ2, _acceptedY,
        _acceptedW, _acceptedZ, _acceptedMomentum, _acceptedTheta});
    }

    void analyze(const Event& event) override {
      const Particles prompt =
        apply<PromptFinalState>(event, "PromptMuons").particles();
      const COMPASSSIDIS::DISKinematics dis =
        COMPASSInclusiveDIS::reconstruct(event, prompt);
      if (!dis.valid || !COMPASSSIDIS::fixed160GeVBeam(dis)) vetoEvent;
      if (!(dis.Q2 > 1. && dis.W2 > 25. && dis.x > .004 && dis.x < .4 &&
            dis.y > .1 && dis.y < .7)) vetoEvent;

      std::map<std::string, std::vector<size_t>> denominatorCells;
      bool accepted = false;
      for (const auto& entry : _cells) {
        for (size_t index = 0; index < entry.second->size(); ++index) {
          if (SIDISTrancheBinning::matchesDIS(
                (*entry.second)[index], dis.x, dis.y)) {
            denominatorCells[entry.first].push_back(index);
            accepted = true;
          }
        }
      }
      if (!accepted) vetoEvent;

      _acceptedX->fill(dis.x);
      _acceptedQ2->fill(dis.Q2);
      _acceptedY->fill(dis.y);
      _acceptedW->fill(std::sqrt(dis.W2));
      std::map<std::string, std::vector<double>> counts;
      for (const auto& entry : _cells)
        counts[entry.first] = std::vector<double>(entry.second->size(), 0.);

      const Particles particles = apply<FinalState>(event, "LabFS").particles();
      for (const Particle& particle : particles) {
        if (!PID::isHadron(particle.pid()) || particle.charge3() == 0) continue;
        const bool positive = particle.charge3() > 0;
        const COMPASSSIDIS::HadronKinematics unidentified =
          COMPASSSIDIS::hadronKinematics(dis, particle, true);
        fillHadron(positive ? "hplus" : "hminus", dis, unidentified, counts);
        if (std::abs(particle.pid()) == 211) {
          const auto pion = COMPASSSIDIS::hadronKinematics(dis, particle);
          fillHadron(particle.pid() > 0 ? "piplus" : "piminus", dis, pion, counts);
        } else if (std::abs(particle.pid()) == 321) {
          const auto kaon = COMPASSSIDIS::hadronKinematics(dis, particle);
          fillHadron(particle.pid() > 0 ? "kplus" : "kminus", dis, kaon, counts);
        }
      }

      for (auto& entry : _species) {
        const std::string& species = entry.first;
        Histograms& histograms = entry.second;
        for (const size_t index : denominatorCells[species])
          histograms.denominator->fill(double(index) + .5);
        for (size_t index = 0; index < counts[species].size(); ++index) {
          const double count = counts[species][index];
          COMPASSSIDIS::fillEventCount(histograms.numerator,
                                       double(index) + .5, count);
          COMPASSSIDIS::fillMultiplicityCovariance(histograms.covariance,
                                                   double(index) + .5, count);
        }
      }
    }

    void finalize() override {
      if (sumW() == 0.) return;
      const double factor = crossSection()/picobarn/sumW();
      for (Histo1DPtr& histogram : _scaled) scale(histogram, factor);
    }

  private:
    using Cells = std::vector<SIDISTrancheBinning::XYZCell>;

    void configureCells() {
      if (ReleaseYear == 2025) {
        _cells = {
          {"hplus", &SIDISTrancheBinning::compass2025_hplusCells()},
          {"hminus", &SIDISTrancheBinning::compass2025_hminusCells()},
          {"piplus", &SIDISTrancheBinning::compass2025_piplusCells()},
          {"piminus", &SIDISTrancheBinning::compass2025_piminusCells()},
          {"kplus", &SIDISTrancheBinning::compass2025_kplusCells()},
          {"kminus", &SIDISTrancheBinning::compass2025_kminusCells()},
        };
      } else {
        _cells = {
          {"hplus", &SIDISTrancheBinning::compass2026_hplusCells()},
          {"hminus", &SIDISTrancheBinning::compass2026_hminusCells()},
          {"piplus", &SIDISTrancheBinning::compass2026_piplusCells()},
          {"piminus", &SIDISTrancheBinning::compass2026_piminusCells()},
          {"kplus", &SIDISTrancheBinning::compass2026_kplusCells()},
          {"kminus", &SIDISTrancheBinning::compass2026_kminusCells()},
        };
      }
    }

    void fillHadron(const std::string& species,
                    const COMPASSSIDIS::DISKinematics& dis,
                    const COMPASSSIDIS::HadronKinematics& hadron,
                    std::map<std::string, std::vector<double>>& counts) {
      if (!hadron.valid || hadron.z < .2 || hadron.z > .85 ||
          hadron.momentum <= 12. || hadron.momentum >= 40. ||
          hadron.theta <= .010 || hadron.theta >= .120) return;
      const int cell = SIDISTrancheBinning::findCell(
        *_cells.at(species), dis.x, dis.y, hadron.z);
      if (cell < 0) return;
      counts[species][size_t(cell)] += 1.;
      _acceptedZ->fill(hadron.z);
      _acceptedMomentum->fill(hadron.momentum);
      _acceptedTheta->fill(hadron.theta);
    }

    std::map<std::string, const Cells*> _cells;
    std::map<std::string, Histograms> _species;
    std::vector<Histo1DPtr> _scaled;
    Histo1DPtr _acceptedX, _acceptedQ2, _acceptedY, _acceptedW;
    Histo1DPtr _acceptedZ, _acceptedMomentum, _acceptedTheta;
  };

} // namespace Rivet
