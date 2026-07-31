// -*- C++ -*-
#pragma once

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Tools/RivetHepMC.hh"
#include "fastjet/ClusterSequence.hh"
#include <algorithm>
#include <cmath>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

namespace Rivet {
namespace STARPolarizedJets {

  struct ClusteredJet {
    FourMomentum momentum;
    double pT() const { return momentum.pT(); }
    double eta() const { return momentum.eta(); }
  };

  struct PartonSelection {
    bool provenanceEstablished = false;
    std::vector<ConstGenParticlePtr> partons;
  };

  template <typename...>
  using VoidT = void;

  template <typename EventT, typename = void>
  struct HasSignalProcessVertexAccessor : std::false_type {};

  template <typename EventT>
  struct HasSignalProcessVertexAccessor<
    EventT,
    VoidT<decltype(std::declval<const EventT&>().signal_process_vertex())>>
    : std::true_type {};

  inline bool isParton(int pid) {
    const int apid = std::abs(pid);
    return (apid >= 1 && apid <= 6) || pid == 21;
  }

  inline bool isRemnant(int pid) {
    return pid == 82;
  }

  inline ConstGenVertexPtr normalizeSignalVertex(
      const ConstGenVertexPtr& vertex) {
    return vertex;
  }

  inline ConstGenVertexPtr normalizeSignalVertex(
      const std::shared_ptr<RivetHepMC::GenVertex>& vertex) {
    return vertex;
  }

  inline ConstGenVertexPtr normalizeSignalVertex(
      const RivetHepMC::GenVertex* vertex) {
    return vertex ? vertex->shared_from_this() : ConstGenVertexPtr();
  }

  template <typename EventT>
  inline ConstGenVertexPtr directSignalVertex(
      const EventT& event, std::true_type) {
    return normalizeSignalVertex(event.signal_process_vertex());
  }

  template <typename EventT>
  inline ConstGenVertexPtr directSignalVertex(
      const EventT&, std::false_type) {
    return ConstGenVertexPtr();
  }

  inline ConstGenVertexPtr signalVertex(const GenEvent& event) {
    const ConstGenVertexPtr direct = directSignalVertex(
      event, HasSignalProcessVertexAccessor<GenEvent>{});
    if (direct) return direct;
    const auto attribute =
      event.attribute<RivetHepMC::IntAttribute>("signal_process_vertex");
    if (!attribute) return ConstGenVertexPtr();
    for (ConstGenVertexPtr vertex : HepMCUtils::vertices(&event)) {
      if (vertex && vertex->id() == attribute->value()) return vertex;
    }
    return ConstGenVertexPtr();
  }

  inline std::vector<ConstGenParticlePtr> children(
      ConstGenParticlePtr particle) {
    std::vector<ConstGenParticlePtr> output;
    if (!particle || !particle->end_vertex()) return output;
    const int parentId = HepMCUtils::uniqueId(particle);
    for (ConstGenParticlePtr child : particle->end_vertex()->particles_out()) {
      if (!child || HepMCUtils::uniqueId(child) == parentId) continue;
      output.push_back(child);
    }
    return output;
  }

  inline bool hasRemnantAncestor(ConstGenParticlePtr particle) {
    std::unordered_set<int> visited;
    std::vector<ConstGenParticlePtr> stack{particle};
    while (!stack.empty()) {
      ConstGenParticlePtr current = stack.back();
      stack.pop_back();
      if (!current || !current->production_vertex()) continue;
      for (ConstGenParticlePtr parent :
           current->production_vertex()->particles_in()) {
        if (!parent) continue;
        const int id = HepMCUtils::uniqueId(parent);
        if (!visited.insert(id).second) continue;
        if (isRemnant(parent->pdg_id())) return true;
        stack.push_back(parent);
      }
    }
    return false;
  }

  inline void terminalPartons(
      ConstGenParticlePtr particle,
      std::unordered_set<int>& accepted,
      std::vector<ConstGenParticlePtr>& output) {
    if (!particle || !isParton(particle->pdg_id()) ||
        hasRemnantAncestor(particle)) {
      return;
    }
    std::vector<ConstGenParticlePtr> partonChildren;
    for (ConstGenParticlePtr child : children(particle)) {
      if (child && isParton(child->pdg_id()) &&
          !hasRemnantAncestor(child)) {
        partonChildren.push_back(child);
      }
    }
    if (partonChildren.empty()) {
      const int id = HepMCUtils::uniqueId(particle);
      if (accepted.insert(id).second) output.push_back(particle);
      return;
    }
    for (ConstGenParticlePtr child : partonChildren) {
      terminalPartons(child, accepted, output);
    }
  }

  /// Select terminal shower partons attached to the hard signal process.
  ///
  /// No all-event fallback is permitted: without a resolvable signal-process
  /// vertex the requested exclusion of beam-remnant descendants is not
  /// demonstrable, and the analysis records a provenance failure.
  inline PartonSelection selectHardShowerPartons(const Event& event) {
    PartonSelection result;
    const GenEvent* generated = event.genEvent();
    if (!generated) return result;
    const ConstGenVertexPtr signal = signalVertex(*generated);
    if (!signal) return result;
    // Only outgoing hard-process lines may seed final-state jet inputs.
    // The explicitly identified signal vertex establishes their provenance.
    // Incoming beam lines (which Herwig's HepMC conversion does not always
    // attach to that vertex) must never be traversed into the clustered
    // final state.
    std::unordered_set<int> seedIds;
    std::vector<ConstGenParticlePtr> seeds;
    for (ConstGenParticlePtr outgoing : signal->particles_out()) {
      if (!outgoing || !isParton(outgoing->pdg_id())) continue;
      const int id = HepMCUtils::uniqueId(outgoing);
      if (seedIds.insert(id).second) seeds.push_back(outgoing);
    }
    if (seeds.empty()) return result;

    std::unordered_set<int> terminalIds;
    for (ConstGenParticlePtr seed : seeds) {
      terminalPartons(seed, terminalIds, result.partons);
    }
    result.provenanceEstablished = !result.partons.empty();
    return result;
  }

  inline std::vector<ClusteredJet> cluster(
      const std::vector<FourMomentum>& momenta, double radius) {
    std::vector<fastjet::PseudoJet> inputs;
    inputs.reserve(momenta.size());
    for (const FourMomentum& momentum : momenta) {
      inputs.emplace_back(
        momentum.px(), momentum.py(), momentum.pz(), momentum.E());
    }
    if (inputs.empty()) return {};
    const fastjet::JetDefinition definition(
      fastjet::antikt_algorithm, radius, fastjet::E_scheme, fastjet::Best);
    const fastjet::ClusterSequence sequence(inputs, definition);
    const std::vector<fastjet::PseudoJet> clustered =
      fastjet::sorted_by_pt(sequence.inclusive_jets());
    std::vector<ClusteredJet> jets;
    jets.reserve(clustered.size());
    for (const fastjet::PseudoJet& jet : clustered) {
      jets.push_back({FourMomentum(jet.E(), jet.px(), jet.py(), jet.pz())});
    }
    return jets;
  }

  inline std::vector<ClusteredJet> hardPartonJets(
      const Event& event, double radius, bool& valid) {
    const PartonSelection selection = selectHardShowerPartons(event);
    valid = selection.provenanceEstablished;
    std::vector<FourMomentum> inputs;
    inputs.reserve(selection.partons.size());
    for (ConstGenParticlePtr particle : selection.partons) {
      inputs.push_back(Particle(particle).momentum());
    }
    return valid ? cluster(inputs, radius) : std::vector<ClusteredJet>();
  }

  inline std::vector<ClusteredJet> stableParticleJets(
      const Particles& particles, double radius) {
    std::vector<FourMomentum> inputs;
    inputs.reserve(particles.size());
    for (const Particle& particle : particles) {
      inputs.push_back(particle.momentum());
    }
    return cluster(inputs, radius);
  }

  inline double deltaPhi(const ClusteredJet& first, const ClusteredJet& second) {
    double value = std::abs(
      first.momentum.phi()-second.momentum.phi());
    if (value > M_PI) value = 2.0*M_PI-value;
    return value;
  }

} // namespace STARPolarizedJets
} // namespace Rivet
