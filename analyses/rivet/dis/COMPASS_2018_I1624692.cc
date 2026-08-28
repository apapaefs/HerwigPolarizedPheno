// -*- C++ -*-
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

  /// COMPASS isoscalar unidentified-hadron multiplicity versus z and PhT^2.
  class COMPASS_2018_I1624692 : public Analysis {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(COMPASS_2018_I1624692);

    struct Histograms { Histo1DPtr numerator, denominator, covariance; };

    void init() {
      declare(FinalState(), "LabFS");
      declare(PromptFinalState(Cuts::pid == -13), "PromptMuons");
      _cells = {
        {"hplus", &SIDISTrancheBinning::compass2018_hplusCells()},
        {"hminus", &SIDISTrancheBinning::compass2018_hminusCells()},
      };
      for (const auto& entry : _cells) {
        const size_t count = entry.second->size();
        Histograms& histograms = _histograms[entry.first];
        book(histograms.numerator,
             "HadronNumerator_" + entry.first + "_cells", count, 0., count);
        book(histograms.denominator,
             "DISDenominator_" + entry.first + "_cells", count, 0., count);
        book(histograms.covariance,
             "CovarianceProxy_" + entry.first + "_cells", count, 0., count);
        _scaled.insert(_scaled.end(), {histograms.numerator,
          histograms.denominator, histograms.covariance});
      }
      book(_acceptedX, "Accepted_X", logspace(48, .003, .4));
      book(_acceptedQ2, "Accepted_Q2", logspace(45, 1., 81.));
      book(_acceptedY, "Accepted_Y", 40, .1, .9);
      book(_acceptedW, "Accepted_W", 40, 5., 20.);
      book(_acceptedZ, "Accepted_Z", 40, .2, .8);
      book(_acceptedPt2, "Accepted_HadronPt2", logspace(50, .02, 3.));
      _scaled.insert(_scaled.end(), {_acceptedX, _acceptedQ2, _acceptedY,
        _acceptedW, _acceptedZ, _acceptedPt2});
    }

    void analyze(const Event& event) {
      const Particles prompt =
        apply<PromptFinalState>(event, "PromptMuons").particles();
      const COMPASSSIDIS::DISKinematics dis =
        COMPASSInclusiveDIS::reconstruct(event, prompt);
      if (!dis.valid || !COMPASSSIDIS::fixed160GeVBeam(dis)) vetoEvent;
      if (!(dis.Q2 > 1. && dis.W2 > 25. && dis.x > .003 && dis.x < .4 &&
            dis.y > .1 && dis.y < .9)) vetoEvent;

      std::map<std::string, std::vector<size_t>> denominatorCells;
      bool accepted = false;
      for (const auto& entry : _cells) {
        for (size_t index = 0; index < entry.second->size(); ++index) {
          if (SIDISTrancheBinning::matchesDIS(
                (*entry.second)[index], dis.x, dis.Q2)) {
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
        const std::string species = particle.charge3() > 0 ? "hplus" : "hminus";
        const auto hadron = COMPASSSIDIS::hadronKinematics(dis, particle, true);
        if (!hadron.valid || hadron.z < .2 || hadron.z > .8 ||
            hadron.transverseMomentum2 < .02 || hadron.transverseMomentum2 > 3.)
          continue;
        const int cell = SIDISTrancheBinning::findCell(
          *_cells.at(species), dis.x, dis.Q2, hadron.z,
          hadron.transverseMomentum2);
        if (cell < 0) continue;
        counts[species][size_t(cell)] += 1.;
        _acceptedZ->fill(hadron.z);
        _acceptedPt2->fill(hadron.transverseMomentum2);
      }

      for (auto& entry : _histograms) {
        for (const size_t index : denominatorCells[entry.first])
          entry.second.denominator->fill(double(index) + .5);
        for (size_t index = 0; index < counts[entry.first].size(); ++index) {
          COMPASSSIDIS::fillEventCount(entry.second.numerator,
                                       double(index) + .5,
                                       counts[entry.first][index]);
          COMPASSSIDIS::fillMultiplicityCovariance(entry.second.covariance,
                                                   double(index) + .5,
                                                   counts[entry.first][index]);
        }
      }
    }

    void finalize() {
      if (sumW() == 0.) return;
      const double factor = crossSection()/picobarn/sumW();
      for (Histo1DPtr& histogram : _scaled) scale(histogram, factor);
    }

  private:
    using Cells = std::vector<SIDISTrancheBinning::XQ2ZTransverseCell>;
    std::map<std::string, const Cells*> _cells;
    std::map<std::string, Histograms> _histograms;
    std::vector<Histo1DPtr> _scaled;
    Histo1DPtr _acceptedX, _acceptedQ2, _acceptedY, _acceptedW;
    Histo1DPtr _acceptedZ, _acceptedPt2;
  };

  RIVET_DECLARE_PLUGIN(COMPASS_2018_I1624692);
}
