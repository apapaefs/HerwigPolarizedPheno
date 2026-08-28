// -*- C++ -*-
#include "COMPASSSIDIS.hh"
#include "SIDISTrancheBinning.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include <cmath>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace Rivet {

  /// COMPASS high-z antiproton/proton and K-/K+ multiplicity ratios.
  class COMPASS_2020_I1788430 : public Analysis {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(COMPASS_2020_I1788430);

    struct Histograms {
      Histo1DPtr negative, positive, covariance;
    };

    struct Dataset {
      const std::vector<SIDISTrancheBinning::XZMomentumCell>* cells;
      int absolutePid;
    };

    void init() override {
      declare(FinalState(), "LabFS");
      declare(PromptFinalState(Cuts::pid == -13), "PromptMuons");
      _datasets = {
        {"pbar_over_p_xz",
         {&SIDISTrancheBinning::compass2020_pbar_over_p_xzCells(), 2212}},
        {"pbar_over_p_lowx_zp",
         {&SIDISTrancheBinning::compass2020_pbar_over_p_lowx_zpCells(), 2212}},
        {"kminus_over_kplus_lowx_zp",
         {&SIDISTrancheBinning::compass2020_kminus_over_kplus_lowx_zpCells(), 321}},
      };
      for (const auto& entry : _datasets) {
        const std::string& dataset = entry.first;
        const size_t count = entry.second.cells->size();
        Histograms& histograms = _histograms[dataset];
        book(histograms.negative,
             "NegativeHadronYield_" + dataset + "_cells", count, 0., count);
        book(histograms.positive,
             "PositiveHadronYield_" + dataset + "_cells", count, 0., count);
        book(histograms.covariance,
             "ChargeCovarianceProxy_" + dataset + "_cells", count, 0., count);
        _scaled.insert(_scaled.end(), {
          histograms.negative, histograms.positive, histograms.covariance,
        });
      }
      book(_acceptedX, "Accepted_X", logspace(45, .01, .4));
      book(_acceptedQ2, "Accepted_Q2", logspace(45, 1., 60.));
      book(_acceptedY, "Accepted_Y", 45, .1, 1.);
      book(_acceptedW, "Accepted_W", 45, 5., 25.);
      book(_acceptedZ, "Accepted_Z", 48, .5, 1.1);
      book(_acceptedMomentum, "Accepted_HadronMomentum", 40, 20., 60.);
      book(_acceptedTheta, "Accepted_HadronTheta", 45, 0., .180);
      _scaled.insert(_scaled.end(), {
        _acceptedX, _acceptedQ2, _acceptedY, _acceptedW,
        _acceptedZ, _acceptedMomentum, _acceptedTheta,
      });
    }

    void analyze(const Event& event) override {
      const Particles prompt =
        apply<PromptFinalState>(event, "PromptMuons").particles();
      const COMPASSSIDIS::DISKinematics dis =
        COMPASSInclusiveDIS::reconstruct(event, prompt);
      if (!dis.valid || !COMPASSSIDIS::fixed160GeVBeam(dis)) vetoEvent;
      if (!(dis.Q2 > 1. && dis.W2 > 25. && dis.x > .01 && dis.x < .4 &&
            dis.y > .1 && dis.y < 1.)) vetoEvent;

      _acceptedX->fill(dis.x);
      _acceptedQ2->fill(dis.Q2);
      _acceptedY->fill(dis.y);
      _acceptedW->fill(std::sqrt(dis.W2));

      using ChargeCounts = std::pair<std::vector<double>, std::vector<double>>;
      std::map<std::string, ChargeCounts> counts;
      for (const auto& entry : _datasets) {
        const size_t size = entry.second.cells->size();
        counts[entry.first] = {
          std::vector<double>(size, 0.), std::vector<double>(size, 0.)
        };
      }

      const Particles particles = apply<FinalState>(event, "LabFS").particles();
      for (const Particle& particle : particles) {
        const int absolutePid = std::abs(particle.pid());
        if (absolutePid != 321 && absolutePid != 2212) continue;
        const auto hadron = COMPASSSIDIS::hadronKinematics(dis, particle);
        if (!hadron.valid || hadron.theta < 0. || hadron.theta >= .180) continue;
        bool accepted = false;
        for (const auto& entry : _datasets) {
          if (absolutePid != entry.second.absolutePid) continue;
          const int cell = SIDISTrancheBinning::findCell(
            *entry.second.cells, dis.x, hadron.z, hadron.momentum
          );
          if (cell < 0) continue;
          ChargeCounts& datasetCounts = counts[entry.first];
          std::vector<double>& chargeCounts =
            particle.pid() < 0 ? datasetCounts.first : datasetCounts.second;
          chargeCounts[size_t(cell)] += 1.;
          accepted = true;
        }
        if (accepted) {
          _acceptedZ->fill(hadron.z);
          _acceptedMomentum->fill(hadron.momentum);
          _acceptedTheta->fill(hadron.theta);
        }
      }

      for (auto& entry : _histograms) {
        const auto& negative = counts[entry.first].first;
        const auto& positive = counts[entry.first].second;
        for (size_t index = 0; index < negative.size(); ++index) {
          const double coordinate = double(index) + .5;
          COMPASSSIDIS::fillEventCount(
            entry.second.negative, coordinate, negative[index]
          );
          COMPASSSIDIS::fillEventCount(
            entry.second.positive, coordinate, positive[index]
          );
          COMPASSSIDIS::fillCountCovariance(
            entry.second.covariance, coordinate,
            negative[index], positive[index]
          );
        }
      }
    }

    void finalize() override {
      if (sumW() == 0.) return;
      const double factor = crossSection()/picobarn/sumW();
      for (Histo1DPtr& histogram : _scaled) scale(histogram, factor);
    }

  private:
    std::map<std::string, Dataset> _datasets;
    std::map<std::string, Histograms> _histograms;
    std::vector<Histo1DPtr> _scaled;
    Histo1DPtr _acceptedX, _acceptedQ2, _acceptedY, _acceptedW;
    Histo1DPtr _acceptedZ, _acceptedMomentum, _acceptedTheta;
  };

  RIVET_DECLARE_PLUGIN(COMPASS_2020_I1788430);
}
