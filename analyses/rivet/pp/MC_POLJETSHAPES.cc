// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Math/LorentzTrans.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VisibleFinalState.hh"

#include "fastjet/ClusterSequence.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <string>
#include <vector>

namespace Rivet {

  namespace {

    constexpr double EPS = 1.0e-12;
    constexpr size_t NANGLE = 24;

    double clampUnit(double value) {
      return std::max(-1.0, std::min(1.0, value));
    }

    double wrapToPi(double angle) {
      while (angle <= -M_PI) angle += 2.0*M_PI;
      while (angle > M_PI) angle -= 2.0*M_PI;
      return angle;
    }

    Vector3 pseudoVector(const fastjet::PseudoJet& momentum) {
      return Vector3(momentum.px(), momentum.py(), momentum.pz());
    }

    bool planeNormal(const Vector3& first, const Vector3& second,
                     Vector3& normal) {
      normal = first.cross(second);
      if (!std::isfinite(normal.mod()) || normal.mod() <= EPS) return false;
      normal = normal.unit();
      return true;
    }

    bool signedPlaneAngle(const Vector3& previous, const Vector3& current,
                          const Vector3& axis, double& angle) {
      if (previous.mod() <= EPS || current.mod() <= EPS || axis.mod() <= EPS) {
        return false;
      }
      const Vector3 first = previous.unit();
      const Vector3 second = current.unit();
      const Vector3 direction = axis.unit();
      angle = std::atan2(direction.dot(first.cross(second)),
                         clampUnit(first.dot(second)));
      return std::isfinite(angle);
    }

    struct PlaneSplit {
      fastjet::PseudoJet hard;
      fastjet::PseudoJet soft;
      Vector3 normal;
      double z = 0.0;
      double kt = 0.0;
      double psi = 0.0;
    };

    std::vector<PlaneSplit> declusterHardBranch(
        const fastjet::PseudoJet& root, const Vector3& referenceNormal,
        double referencePsi) {
      std::vector<PlaneSplit> sequence;
      fastjet::PseudoJet current = root;
      Vector3 previousNormal = referenceNormal;
      double psi = referencePsi;
      fastjet::PseudoJet first, second;
      while (current.has_parents(first, second)) {
        if (first.perp() < second.perp()) std::swap(first, second);
        const double ptSum = first.perp() + second.perp();
        if (ptSum <= EPS) break;
        Vector3 normal;
        if (!planeNormal(pseudoVector(first), pseudoVector(second), normal)) {
          current = first;
          continue;
        }
        double increment = 0.0;
        if (signedPlaneAngle(previousNormal, normal, pseudoVector(first),
                             increment)) {
          psi += increment;
        }
        PlaneSplit split;
        split.hard = first;
        split.soft = second;
        split.normal = normal;
        split.z = second.perp()/ptSum;
        split.kt = second.perp()*first.delta_R(second);
        split.psi = psi;
        sequence.push_back(split);
        previousNormal = normal;
        current = first;
      }
      return sequence;
    }

    bool highestKtSplit(const std::vector<PlaneSplit>& sequence,
                        double zMin, double ktMin, PlaneSplit& selected) {
      bool found = false;
      for (const PlaneSplit& split : sequence) {
        if (split.z <= zMin || split.kt <= ktMin) continue;
        if (!found || split.kt > selected.kt) {
          selected = split;
          found = true;
        }
      }
      return found;
    }

    struct JetAngles {
      bool primary05 = false;
      bool primary10 = false;
      bool loose = false;
      bool symmetric = false;
      bool perturbative = false;
      double primaryPsi05 = 0.0;
      double primaryPsi10 = 0.0;
      double loosePsi = 0.0;
      double symmetricPsi = 0.0;
      double perturbativePsi = 0.0;
    };

    JetAngles declusteringAngles(const Jet& jet, double radius) {
      JetAngles result;
      std::vector<fastjet::PseudoJet> particles;
      particles.reserve(jet.particles().size());
      for (const Particle& particle : jet.particles()) {
        const FourMomentum& momentum = particle.momentum();
        if (particle.pT() <= 0.0) continue;
        particles.emplace_back(momentum.px()/GeV, momentum.py()/GeV,
                               momentum.pz()/GeV, momentum.E()/GeV);
      }
      if (particles.size() < 2) return result;

      const fastjet::JetDefinition definition(
          fastjet::cambridge_algorithm, radius, fastjet::E_scheme);
      const fastjet::ClusterSequence clustering(particles, definition);
      std::vector<fastjet::PseudoJet> roots =
          fastjet::sorted_by_pt(clustering.inclusive_jets(0.0));
      if (roots.empty()) return result;

      Vector3 referenceNormal;
      if (!planeNormal(Vector3::mkZ(), jet.p3(), referenceNormal)) {
        return result;
      }
      const std::vector<PlaneSplit> primary = declusterHardBranch(
          roots.front(), referenceNormal, 0.0);
      PlaneSplit selected05, selected10;
      result.primary05 = highestKtSplit(primary, 0.1, 0.5, selected05);
      result.primary10 = highestKtSplit(primary, 0.1, 1.0, selected10);
      if (result.primary05) {
        result.primaryPsi05 = wrapToPi(selected05.psi);
        const std::vector<PlaneSplit> secondary = declusterHardBranch(
            selected05.soft, selected05.normal, selected05.psi);
        PlaneSplit loose, symmetric;
        result.loose = highestKtSplit(secondary, 0.1, 0.5, loose);
        result.symmetric = highestKtSplit(secondary, 0.3, 0.5, symmetric);
        if (result.loose) {
          result.loosePsi = wrapToPi(loose.psi - selected05.psi);
        }
        if (result.symmetric) {
          result.symmetricPsi = wrapToPi(
              symmetric.psi - selected05.psi);
        }
      }
      if (result.primary10) {
        result.primaryPsi10 = wrapToPi(selected10.psi);
        const std::vector<PlaneSplit> secondary = declusterHardBranch(
            selected10.soft, selected10.normal, selected10.psi);
        PlaneSplit perturbative;
        result.perturbative = highestKtSplit(
            secondary, 0.1, 1.0, perturbative);
        if (result.perturbative) {
          result.perturbativePsi = wrapToPi(
              perturbative.psi - selected10.psi);
        }
      }
      return result;
    }

    size_t angleBin(double angle) {
      const double wrapped = wrapToPi(angle);
      const double scaled = (wrapped + M_PI)/(2.0*M_PI);
      return std::min(NANGLE - 1,
                      static_cast<size_t>(std::floor(NANGLE*scaled)));
    }

  }


  /// Particle-level jet shapes designed to expose shower-spin correlations.
  class MC_POLJETSHAPES : public Analysis {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(MC_POLJETSHAPES);

    void init() {
      _radius = getOption<double>("R", 0.5);
      _leadingPtMin = getOption<double>("PTJ1MIN", 5.0)*GeV;
      _subleadingPtMin = getOption<double>("PTJ2MIN", 4.0)*GeV;
      _etaMax = getOption<double>("ETAMAX", 1.5);

      const VisibleFinalState visible;
      declare(FastJets(visible, JetAlg::ANTIKT, _radius,
                       JetMuons::ALL, JetInvisibles::NONE), "Jets");

      _angularNames = {
        "dpsi12_j1_loose", "dpsi12_j2_loose",
        "dpsi12_j1_symmetric_secondary",
        "dpsi12_j2_symmetric_secondary",
        "dpsi12_j1_perturbative", "dpsi12_j2_perturbative",
        "hardplane_primary_j1_kt05", "hardplane_primary_j2_kt05",
        "hardplane_primary_j1_kt10", "hardplane_primary_j2_kt10",
        "interjet_dpsi11_kt05", "interjet_dpsi11_kt10",
        "eeec_squeezed_j1", "eeec_squeezed_j2", "bz_angle"
      };
      for (const std::string& name : _angularNames) {
        if (name == "bz_angle") {
          book(_physics[name], name + "_Yield", NANGLE, 0.0, M_PI);
        } else {
          book(_physics[name], name + "_Yield", NANGLE, -M_PI, M_PI);
        }
        _acceptedIndex[name] = _acceptedIndex.size();
      }

      for (const std::string jet : {"j1", "j2"}) {
        for (const std::string beta : {"beta1", "beta2"}) {
          book(_physics["q2_" + beta + "_" + jet],
               "q2_" + beta + "_" + jet + "_Yield", 24, -1.0, 1.0);
          book(_physics["s2_" + beta + "_" + jet],
               "s2_" + beta + "_" + jet + "_Yield", 24, -1.0, 1.0);
        }
      }

      const std::vector<double> rateEdges =
        {1.5, 2.5, 3.5, 4.5, 5.5, 7.0, 9.0, 11.0};
      const std::vector<double> tailEdges =
        {0.05, 0.15, 0.25, 0.35, 0.45, 0.55};
      book(_physics["dijet_threshold_denominator"],
           "dijet_threshold_denominator_Yield", rateEdges);
      book(_physics["ge3_threshold"], "ge3_threshold_Yield", rateEdges);
      book(_physics["ge4_threshold"], "ge4_threshold_Yield", rateEdges);
      book(_physics["pt31_cumulative_tail"],
           "pt31_cumulative_tail_Yield", tailEdges);
      book(_physics["pt41_cumulative_tail"],
           "pt41_cumulative_tail_Yield", tailEdges);

      book(_cutflow, "CutflowFraction", 8, 0.0, 8.0);
      book(_accepted, "AcceptedEntriesPerEvent",
           _angularNames.size(), 0.0, double(_angularNames.size()));
    }

    void analyze(const Event& event) {
      _cutflow->fill(0.5);
      const Jets jets = apply<FastJets>(event, "Jets").jetsByPt(
          Cuts::abseta < _etaMax);
      if (jets.size() < 2) return;
      _cutflow->fill(1.5);
      if (jets[0].pT() <= _leadingPtMin ||
          jets[1].pT() <= _subleadingPtMin) return;
      _cutflow->fill(2.5);

      const std::array<double, 7> thresholds = {{2, 3, 4, 5, 6, 8, 10}};
      for (double threshold : thresholds) {
        _physics["dijet_threshold_denominator"]->fill(threshold);
        if (jets.size() >= 3 && jets[2].pT()/GeV > threshold) {
          _physics["ge3_threshold"]->fill(threshold);
        }
        if (jets.size() >= 4 && jets[3].pT()/GeV > threshold) {
          _physics["ge4_threshold"]->fill(threshold);
        }
      }
      const std::array<double, 5> ratioCuts = {{0.1, 0.2, 0.3, 0.4, 0.5}};
      if (jets.size() >= 3) {
        const double ratio = jets[2].pT()/jets[0].pT();
        for (double cut : ratioCuts) {
          if (ratio > cut) _physics["pt31_cumulative_tail"]->fill(cut);
        }
      }
      if (jets.size() >= 4) {
        const double ratio = jets[3].pT()/jets[0].pT();
        for (double cut : ratioCuts) {
          if (ratio > cut) _physics["pt41_cumulative_tail"]->fill(cut);
        }
      }

      const JetAngles first = declusteringAngles(jets[0], _radius);
      const JetAngles second = declusteringAngles(jets[1], _radius);
      fillJetAngles(first, "j1");
      fillJetAngles(second, "j2");
      if (first.primary05) _cutflow->fill(3.5);
      if (first.loose) _cutflow->fill(4.5);
      if (first.primary05 && second.primary05) {
        _cutflow->fill(5.5);
        fillAngular("interjet_dpsi11_kt05", wrapToPi(
            first.primaryPsi05 - second.primaryPsi05));
        _cutflow->fill(6.5);
      }
      if (first.primary10 && second.primary10) {
        fillAngular("interjet_dpsi11_kt10", wrapToPi(
            first.primaryPsi10 - second.primaryPsi10));
      }

      fillQuadrupoles(jets[0], "j1");
      fillQuadrupoles(jets[1], "j2");
      fillEEEC(jets[0], "eeec_squeezed_j1");
      fillEEEC(jets[1], "eeec_squeezed_j2");

      if (jets.size() >= 4 && jets[2].pT() > 2.0*GeV &&
          jets[3].pT() > 2.0*GeV) {
        _cutflow->fill(7.5);
        double angle = 0.0;
        if (bengtssonZerwas(jets, angle)) fillAngular("bz_angle", angle);
      }
    }

    void finalize() {
      if (sumW() == 0.0) return;
      const double crossSectionFactor = crossSection()/picobarn/sumW();
      for (auto& entry : _physics) scale(entry.second, crossSectionFactor);
      scale(_cutflow, 1.0/sumW());
      scale(_accepted, 1.0/sumW());
    }

  private:
    void fillAngular(const std::string& name, double angle,
                     double weight = 1.0, double entries = 1.0) {
      if (!std::isfinite(angle) || !std::isfinite(weight)) return;
      _physics[name]->fill(angle, weight);
      const auto found = _acceptedIndex.find(name);
      if (found != _acceptedIndex.end() && entries > 0.0) {
        _accepted->fill(double(found->second) + 0.5, entries);
      }
    }

    void fillJetAngles(const JetAngles& angles, const std::string& jet) {
      if (angles.loose) {
        fillAngular("dpsi12_" + jet + "_loose", angles.loosePsi);
      }
      if (angles.symmetric) {
        fillAngular("dpsi12_" + jet + "_symmetric_secondary",
                    angles.symmetricPsi);
      }
      if (angles.perturbative) {
        fillAngular("dpsi12_" + jet + "_perturbative",
                    angles.perturbativePsi);
      }
      if (angles.primary05) {
        fillAngular("hardplane_primary_" + jet + "_kt05",
                    angles.primaryPsi05);
      }
      if (angles.primary10) {
        fillAngular("hardplane_primary_" + jet + "_kt10",
                    angles.primaryPsi10);
      }
    }

    void fillQuadrupoles(const Jet& jet, const std::string& label) {
      const Vector3 axis = jet.p3().unit();
      Vector3 e1 = Vector3::mkZ() - axis*axis.dot(Vector3::mkZ());
      if (e1.mod() <= EPS) return;
      e1 = e1.unit();
      const Vector3 e2 = axis.cross(e1).unit();
      double sumPt = 0.0;
      for (const Particle& particle : jet.particles()) {
        sumPt += particle.pT();
      }
      if (sumPt <= 0.0) return;
      double q1 = 0.0, s1 = 0.0, q2 = 0.0, s2 = 0.0;
      for (const Particle& particle : jet.particles()) {
        const double z = particle.pT()/sumPt;
        const double radius = deltaR(particle.momentum(), jet.momentum())/_radius;
        Vector3 transverse = particle.p3().unit();
        transverse -= axis*transverse.dot(axis);
        if (transverse.mod() <= EPS) continue;
        const double phi = std::atan2(transverse.dot(e2), transverse.dot(e1));
        const double cosine = std::cos(2.0*phi);
        const double sine = std::sin(2.0*phi);
        q1 += z*radius*cosine;
        s1 += z*radius*sine;
        q2 += z*radius*radius*cosine;
        s2 += z*radius*radius*sine;
      }
      _physics["q2_beta1_" + label]->fill(q1);
      _physics["s2_beta1_" + label]->fill(s1);
      _physics["q2_beta2_" + label]->fill(q2);
      _physics["s2_beta2_" + label]->fill(s2);
    }

    void fillEEEC(const Jet& jet, const std::string& name) {
      const Particles& particles = jet.particles();
      if (particles.size() < 3) return;
      double sumPt = 0.0;
      for (const Particle& particle : particles) sumPt += particle.pT();
      if (sumPt <= 0.0) return;

      std::array<double, NANGLE> eventWeights = {{0.0}};
      size_t acceptedTriplets = 0;
      for (size_t i = 0; i + 2 < particles.size(); ++i) {
        for (size_t j = i + 1; j + 1 < particles.size(); ++j) {
          for (size_t k = j + 1; k < particles.size(); ++k) {
            const std::array<const Particle*, 3> input =
              {{&particles[i], &particles[j], &particles[k]}};
            const std::array<std::array<size_t, 3>, 3> choices = {{
              {{0, 1, 2}}, {{0, 2, 1}}, {{1, 2, 0}}
            }};
            double smallest = 1.0e100;
            std::array<size_t, 3> selected = choices[0];
            for (const auto& choice : choices) {
              const double distance = deltaR(
                  input[choice[0]]->momentum(), input[choice[1]]->momentum());
              if (distance < smallest) {
                smallest = distance;
                selected = choice;
              }
            }
            const Particle& first = *input[selected[0]];
            const Particle& second = *input[selected[1]];
            const Particle& third = *input[selected[2]];
            const FourMomentum pair = first.momentum() + second.momentum();
            const double thetaS = smallest/_radius;
            const double thetaL = deltaR(pair, third.momentum())/_radius;
            if (!(thetaS > 0.02 && thetaS < 0.10 &&
                  thetaL > std::sqrt(0.1) && thetaL < 1.0)) continue;

            Vector3 smallNormal, largeNormal;
            if (!planeNormal(first.p3(), second.p3(), smallNormal) ||
                !planeNormal(pair.p3(), third.p3(), largeNormal)) continue;
            double angle = 0.0;
            if (!signedPlaneAngle(smallNormal, largeNormal, pair.p3(), angle)) {
              continue;
            }
            const double zFirst = first.pT()/sumPt;
            const double zSecond = second.pT()/sumPt;
            const double zThird = third.pT()/sumPt;
            eventWeights[angleBin(angle)] +=
                6.0*zFirst*zSecond*zThird;
            ++acceptedTriplets;
          }
        }
      }
      if (acceptedTriplets == 0) return;
      for (size_t bin = 0; bin < NANGLE; ++bin) {
        if (eventWeights[bin] == 0.0) continue;
        const double centre = -M_PI + (double(bin) + 0.5)*2.0*M_PI/NANGLE;
        _physics[name]->fill(centre, eventWeights[bin]);
      }
      const auto found = _acceptedIndex.find(name);
      if (found != _acceptedIndex.end()) {
        _accepted->fill(double(found->second) + 0.5,
                        double(acceptedTriplets));
      }
    }

    bool bengtssonZerwas(const Jets& jets, double& angle) const {
      FourMomentum total;
      for (size_t index = 0; index < 4; ++index) {
        total += jets[index].momentum();
      }
      if (total.mass2() <= 0.0 || total.E() <= 0.0) return false;
      const LorentzTransform toCentreOfMass =
          LorentzTransform::mkFrameTransform(total);
      std::array<FourMomentum, 4> momenta;
      for (size_t index = 0; index < 4; ++index) {
        momenta[index] = toCentreOfMass.transform(jets[index].momentum());
      }
      std::sort(momenta.begin(), momenta.end(),
                [](const FourMomentum& first, const FourMomentum& second) {
                  return first.E() > second.E();
                });
      Vector3 firstNormal, secondNormal;
      if (!planeNormal(momenta[0].p3(), momenta[1].p3(), firstNormal) ||
          !planeNormal(momenta[2].p3(), momenta[3].p3(), secondNormal)) {
        return false;
      }
      angle = std::acos(clampUnit(firstNormal.dot(secondNormal)));
      return std::isfinite(angle);
    }

    double _radius = 0.5;
    double _leadingPtMin = 5.0*GeV;
    double _subleadingPtMin = 4.0*GeV;
    double _etaMax = 1.5;
    std::map<std::string, Histo1DPtr> _physics;
    Histo1DPtr _cutflow, _accepted;
    std::vector<std::string> _angularNames;
    std::map<std::string, size_t> _acceptedIndex;
  };

  RIVET_DECLARE_PLUGIN(MC_POLJETSHAPES);
}
