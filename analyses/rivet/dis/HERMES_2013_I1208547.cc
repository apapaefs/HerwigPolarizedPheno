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

  /// HERMES VM-subtracted pion and kaon multiplicities in five independent
  /// three-dimensional binnings.  Only event-aggregated raw estimators are
  /// written here; target combinations, ratios, widths and projections are
  /// formed after signed-NLO and shard combination.
  class HERMES_2013_I1208547 : public Analysis {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(HERMES_2013_I1208547);

    struct SpeciesHistograms {
      Histo1DPtr numerator, covariance;
    };
    struct ConfigurationHistograms {
      Histo1DPtr denominator;
      std::map<std::string, SpeciesHistograms> species;
    };
    struct ProjectionHistograms {
      Histo1DPtr denominator;
      std::map<std::string, SpeciesHistograms> species;
      std::vector<std::vector<size_t>> cellGroups;
      std::vector<double> edges;
    };

    void init() {
      declare(FinalState(), "LabFS");
      declare(PromptFinalState(Cuts::abspid == PID::ELECTRON),
              "PromptLeptons");
      _cells = {
        {"z-3D", &SIDISTrancheBinning::hermes_z_3dCells()},
        {"zpt-3D", &SIDISTrancheBinning::hermes_zpt_3dCells()},
        {"zx-3D", &SIDISTrancheBinning::hermes_zx_3dCells()},
        {"zQ2-3D", &SIDISTrancheBinning::hermes_zq2_3dCells()},
        {"zxpt-3D", &SIDISTrancheBinning::hermes_zxpt_3dCells()},
      };
      for (const auto& configuration : _cells) {
        const std::string& name = configuration.first;
        const size_t count = configuration.second->size();
        ConfigurationHistograms& histograms = _histograms[name];
        book(histograms.denominator,
             "DISDenominator_" + name + "_cells", count, 0., count);
        _scaled.push_back(histograms.denominator);
        for (const std::string species :
             {"piplus", "piminus", "kplus", "kminus"}) {
          SpeciesHistograms& item = histograms.species[species];
          book(item.numerator,
               "HadronNumerator_" + species + "_" + name + "_cells",
               count, 0., count);
          book(item.covariance,
               "CovarianceProxy_" + species + "_" + name + "_cells",
               count, 0., count);
          _scaled.insert(_scaled.end(), {item.numerator, item.covariance});
        }
      }
      // These objects implement the released HERMES projections at the raw
      // yield level.  Numerators and inclusive-DIS denominators are integrated
      // independently on each event; postprocessing forms their ratio only
      // after shards, signed-NLO samples, and target components are combined.
      addProjection("z-3D", "z", {.2,.25,.3,.4,.5,.6,.7,.8,1.1},
                    "z", -1., -1.);
      const std::vector<std::pair<double, double>> projectionZ = {
        {.2,.3}, {.3,.4}, {.4,.6}, {.6,.8}
      };
      for (const auto& interval : projectionZ) {
        const std::string token = zToken(interval.first, interval.second);
        addProjection("zpt-3D", "pt_" + token,
                      {0.,.1,.2,.3,.4,.5,.6,.7,.8,1.2},
                      "transverse", interval.first, interval.second);
        addProjection("zx-3D", "x_" + token,
                      {.023,.04,.055,.075,.1,.14,.2,.3,.4,.6},
                      "x", interval.first, interval.second);
        addProjection("zQ2-3D", "q2_" + token,
                      {1.,1.25,1.5,1.75,2.,2.25,2.5,3.,5.,15.},
                      "q2", interval.first, interval.second);
      }
      book(_acceptedX, "Accepted_X", logspace(48, .023, .6));
      book(_acceptedQ2, "Accepted_Q2", logspace(45, 1., 15.));
      book(_acceptedY, "Accepted_Y", 40, .1, .85);
      book(_acceptedW2, "Accepted_W2", 40, 10., 50.);
      book(_acceptedZ, "Accepted_Z", 50, .2, 1.1);
      book(_acceptedMomentum, "Accepted_HadronMomentum", 45, 2., 15.);
      book(_acceptedPhPerp, "Accepted_HadronPhPerp", 48, 0., 1.2);
      _scaled.insert(_scaled.end(), {_acceptedX, _acceptedQ2, _acceptedY,
        _acceptedW2, _acceptedZ, _acceptedMomentum, _acceptedPhPerp});
    }

    void analyze(const Event& event) {
      const Particles prompt =
        apply<PromptFinalState>(event, "PromptLeptons").particles();
      const COMPASSSIDIS::DISKinematics dis =
        COMPASSInclusiveDIS::reconstructForLeptons(event, prompt, {11, -11});
      if (!dis.valid || !COMPASSSIDIS::fixedBeamEnergy(dis, 27.6)) vetoEvent;
      if (!(dis.Q2 > 1. && dis.W2 > 10. && dis.y > .1 && dis.y < .85))
        vetoEvent;

      std::map<std::string, std::vector<size_t>> denominatorCells;
      bool accepted = false;
      for (const auto& configuration : _cells) {
        for (size_t index = 0; index < configuration.second->size(); ++index) {
          if (SIDISTrancheBinning::matchesDIS(
                (*configuration.second)[index], dis.x, dis.Q2)) {
            denominatorCells[configuration.first].push_back(index);
            accepted = true;
          }
        }
      }
      if (!accepted) vetoEvent;

      _acceptedX->fill(dis.x);
      _acceptedQ2->fill(dis.Q2);
      _acceptedY->fill(dis.y);
      _acceptedW2->fill(dis.W2);
      std::map<std::string, std::map<std::string, std::vector<double>>> counts;
      for (const auto& configuration : _cells)
        for (const std::string species :
             {"piplus", "piminus", "kplus", "kminus"})
          counts[configuration.first][species] =
            std::vector<double>(configuration.second->size(), 0.);

      const Particles particles = apply<FinalState>(event, "LabFS").particles();
      for (const Particle& particle : particles) {
        std::string species;
        if (particle.pid() == 211) species = "piplus";
        if (particle.pid() == -211) species = "piminus";
        if (particle.pid() == 321) species = "kplus";
        if (particle.pid() == -321) species = "kminus";
        if (species.empty()) continue;
        const auto hadron = COMPASSSIDIS::hadronKinematics(dis, particle);
        if (!hadron.valid || hadron.momentum <= 2. || hadron.momentum >= 15.)
          continue;
        for (const auto& configuration : _cells) {
          const int cell = SIDISTrancheBinning::findCell(
            *configuration.second, dis.x, dis.Q2, hadron.z,
            hadron.transverseMomentum);
          if (cell >= 0)
            counts[configuration.first][species][size_t(cell)] += 1.;
        }
        _acceptedZ->fill(hadron.z);
        _acceptedMomentum->fill(hadron.momentum);
        _acceptedPhPerp->fill(hadron.transverseMomentum);
      }

      for (auto& configuration : _histograms) {
        for (const size_t index : denominatorCells[configuration.first])
          configuration.second.denominator->fill(double(index) + .5);
        for (auto& species : configuration.second.species) {
          const auto& eventCounts = counts[configuration.first][species.first];
          for (size_t index = 0; index < eventCounts.size(); ++index) {
            COMPASSSIDIS::fillEventCount(species.second.numerator,
                                         double(index) + .5,
                                         eventCounts[index]);
            COMPASSSIDIS::fillMultiplicityCovariance(species.second.covariance,
                                                     double(index) + .5,
                                                     eventCounts[index]);
          }
        }
      }
      for (auto& configuration : _projections) {
        const Cells& cells = *_cells.at(configuration.first);
        for (auto& projectionItem : configuration.second) {
          ProjectionHistograms& projection = projectionItem.second;
          for (size_t bin = 0; bin < projection.cellGroups.size(); ++bin) {
            const double coordinate = .5*(projection.edges[bin]
                                           + projection.edges[bin+1]);
            bool denominatorAccepted = false;
            for (const size_t cell : projection.cellGroups[bin]) {
              if (SIDISTrancheBinning::matchesDIS(cells[cell], dis.x, dis.Q2)) {
                denominatorAccepted = true;
                break;
              }
            }
            if (denominatorAccepted) projection.denominator->fill(coordinate);
            for (auto& species : projection.species) {
              double count = 0.;
              for (const size_t cell : projection.cellGroups[bin])
                count += counts[configuration.first][species.first][cell];
              COMPASSSIDIS::fillEventCount(
                species.second.numerator, coordinate, count);
              COMPASSSIDIS::fillMultiplicityCovariance(
                species.second.covariance, coordinate, count);
            }
          }
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

    static bool sameEdge(double left, double right) {
      return std::abs(left-right) < 1.e-10;
    }

    static int binIndex(double value, const std::vector<double>& edges) {
      for (size_t bin = 0; bin+1 < edges.size(); ++bin)
        if (value >= edges[bin] &&
            (value < edges[bin+1] ||
             (bin+2 == edges.size() && sameEdge(value, edges[bin+1]))))
          return int(bin);
      return -1;
    }

    static std::string zToken(double low, double high) {
      const int low100 = int(std::lround(100.*low));
      const int high100 = int(std::lround(100.*high));
      const std::string lowText = (low100 < 10 ? "00" : low100 < 100 ? "0" : "")
                                  + std::to_string(low100);
      const std::string highText = (high100 < 10 ? "00" : high100 < 100 ? "0" : "")
                                   + std::to_string(high100);
      return "z" + lowText + "_" + highText;
    }

    void addProjection(const std::string& configuration,
                       const std::string& key,
                       const std::vector<double>& edges,
                       const std::string& axis,
                       double zLow, double zHigh) {
      ProjectionHistograms& projection = _projections[configuration][key];
      projection.edges = edges;
      projection.cellGroups.resize(edges.size()-1);
      const Cells& cells = *_cells.at(configuration);
      for (size_t cellIndex = 0; cellIndex < cells.size(); ++cellIndex) {
        const auto& cell = cells[cellIndex];
        if (zLow >= 0. &&
            (!sameEdge(cell.zLow, zLow) || !sameEdge(cell.zHigh, zHigh)))
          continue;
        const double value = axis == "z" ? .5*(cell.zLow+cell.zHigh)
                           : axis == "x" ? .5*(cell.xLow+cell.xHigh)
                           : axis == "q2" ? .5*(cell.q2Low+cell.q2High)
                           : .5*(cell.transverseLow+cell.transverseHigh);
        const int bin = binIndex(value, edges);
        if (bin >= 0) projection.cellGroups[size_t(bin)].push_back(cellIndex);
      }
      for (const auto& group : projection.cellGroups)
        if (group.empty())
          throw Error("Empty HERMES readable-projection cell group");
      book(projection.denominator,
           "ProjectionDISDenominator_" + configuration + "_" + key, edges);
      _scaled.push_back(projection.denominator);
      for (const std::string species :
           {"piplus", "piminus", "kplus", "kminus"}) {
        SpeciesHistograms& item = projection.species[species];
        book(item.numerator,
             "ProjectionHadronNumerator_" + species + "_" +
             configuration + "_" + key, edges);
        book(item.covariance,
             "ProjectionCovarianceProxy_" + species + "_" +
             configuration + "_" + key, edges);
        _scaled.insert(_scaled.end(), {item.numerator, item.covariance});
      }
    }

    std::map<std::string, const Cells*> _cells;
    std::map<std::string, ConfigurationHistograms> _histograms;
    std::map<std::string, std::map<std::string, ProjectionHistograms>> _projections;
    std::vector<Histo1DPtr> _scaled;
    Histo1DPtr _acceptedX, _acceptedQ2, _acceptedY, _acceptedW2;
    Histo1DPtr _acceptedZ, _acceptedMomentum, _acceptedPhPerp;
  };

  RIVET_DECLARE_PLUGIN(HERMES_2013_I1208547);
}
