// Particle-level plane geometry, matching the frozen MC_POLJETSHAPES definitions.
#ifndef HERWIGPHENO_POLJETSHAPESHARD_HH
#define HERWIGPHENO_POLJETSHAPESHARD_HH
#include "Rivet/Analysis.hh"
#include "fastjet/ClusterSequence.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>
namespace Rivet { namespace PolJetShapesHard {
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
                        double zMin, double ktMin, PlaneSplit& selected,
                        double zMax = 1.0) {
      bool found = false;
      for (const PlaneSplit& split : sequence) {
        if (split.z <= zMin || split.z >= zMax || split.kt <= ktMin) continue;
        if (!found || split.kt > selected.kt) {
          selected = split;
          found = true;
        }
      }
      return found;
    }


    // Half-open, disjoint particle-jet windows. This is not a generation cut.
    inline int ptWindow(double pt) {
      if (pt >= 20.0 && pt < 30.0) return 0;
      if (pt >= 30.0 && pt < 45.0) return 1;
      return -1;
    }
    struct Angles {
      std::array<bool, 2> primary{{false, false}}, pair{{false, false}},
                          hardshare{{false, false}};
      std::array<double, 2> primaryPsi{{0, 0}}, pairPsi{{0, 0}},
                            hardsharePsi{{0, 0}};
    };
    inline Angles angles(const Jet& jet, double radius) {
      Angles result;
      std::vector<fastjet::PseudoJet> particles;
      for (const Particle& particle : jet.particles()) {
        if (particle.pT() <= 0) continue;
        const FourMomentum& p = particle.momentum();
        particles.emplace_back(p.px()/GeV, p.py()/GeV, p.pz()/GeV, p.E()/GeV);
      }
      if (particles.size() < 2) return result;
      const fastjet::JetDefinition definition(fastjet::cambridge_algorithm,
                                              radius, fastjet::E_scheme);
      const fastjet::ClusterSequence cluster(particles, definition);
      const auto roots = fastjet::sorted_by_pt(cluster.inclusive_jets(0));
      Vector3 reference;
      if (roots.empty() || !planeNormal(Vector3(0,0,1), jet.p3(), reference))
        return result;
      const auto primary = declusterHardBranch(roots.front(), reference, 0);
      for (size_t k = 0; k < 2; ++k) {
        const double kt = k+1.0;
        PlaneSplit first, second;
        if (highestKtSplit(primary, 0.1, kt, first)) {
          result.primary[k] = true;
          result.primaryPsi[k] = wrapToPi(first.psi);
          const auto secondary = declusterHardBranch(first.soft, first.normal, first.psi);
          if (highestKtSplit(secondary, 0.1, kt, second)) {
            result.pair[k] = true;
            result.pairPsi[k] = wrapToPi(second.psi-first.psi);
          }
        }
        if (highestKtSplit(primary, 0.25, kt, first, 0.40)) {
          const auto secondary = declusterHardBranch(first.soft, first.normal, first.psi);
          if (highestKtSplit(secondary, 0.35, kt, second)) {
            result.hardshare[k] = true;
            result.hardsharePsi[k] = wrapToPi(second.psi-first.psi);
          }
        }
      }
      return result;
    }
} }
#endif
