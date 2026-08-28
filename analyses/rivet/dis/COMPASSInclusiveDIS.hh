// -*- C++ -*-
#pragma once

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Tools/Beams.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace Rivet {
namespace COMPASSInclusiveDIS {

  struct Kinematics {
    bool valid = false;
    int targetPid = 0;
    int beamPid = 0;
    FourMomentum k, kp, target, q;
    double Q2 = -1.0;
    double x = -1.0;
    double y = -1.0;
    double W2 = -1.0;
    double theta = -1.0;
    double outgoingEnergy = -1.0;
  };

  inline Particle scatteredMuon(const Particles& finalState,
                                const Particle& incoming) {
    Particles candidates;
    for (const Particle& particle : finalState) {
      if (particle.pid() == incoming.pid()) candidates.push_back(particle);
    }
    if (candidates.empty()) return Particle();
    return *std::max_element(
      candidates.begin(), candidates.end(),
      [](const Particle& a, const Particle& b) { return a.E() < b.E(); });
  }

  inline Kinematics reconstructForLeptons(
      const Event& event, const Particles& finalState,
      const std::vector<int>& allowedBeamPids) {
    Kinematics out;
    const ParticlePair incoming = Rivet::beams(event);
    Particle lepton;
    Particle hadron;
    if (PID::isLepton(incoming.first.pid()) &&
        PID::isHadron(incoming.second.pid())) {
      lepton = incoming.first;
      hadron = incoming.second;
    } else if (PID::isLepton(incoming.second.pid()) &&
               PID::isHadron(incoming.first.pid())) {
      lepton = incoming.second;
      hadron = incoming.first;
    } else {
      return out;
    }
    if (std::find(allowedBeamPids.begin(), allowedBeamPids.end(),
                  lepton.pid()) == allowedBeamPids.end() ||
        (hadron.pid() != 2212 && hadron.pid() != 2112)) return out;

    const Particle outgoing = scatteredMuon(finalState, lepton);
    if (outgoing.pid() == PID::ANY) return out;
    const FourMomentum k = lepton.momentum();
    const FourMomentum kp = outgoing.momentum();
    const FourMomentum P = hadron.momentum();
    const FourMomentum q = k - kp;
    const double Pdotq = P * q;
    const double Pdotk = P * k;
    if (Pdotq <= 0.0 || Pdotk <= 0.0) return out;

    out.targetPid = hadron.pid();
    out.beamPid = lepton.pid();
    out.k = k;
    out.kp = kp;
    out.target = P;
    out.q = q;
    out.Q2 = -q.mass2() / GeV2;
    out.x = (-q.mass2()) / (2.0 * Pdotq);
    out.y = Pdotq / Pdotk;
    out.W2 = (P + q).mass2() / GeV2;
    out.theta = k.angle(kp);
    out.outgoingEnergy = kp.E() / GeV;
    out.valid = std::isfinite(out.Q2) && std::isfinite(out.x) &&
                std::isfinite(out.y) && std::isfinite(out.W2) &&
                std::isfinite(out.theta) &&
                std::isfinite(out.outgoingEnergy);
    return out;
  }

  inline Kinematics reconstruct(const Event& event,
                                const Particles& finalState) {
    // Backward-compatible COMPASS contract: lepton.pid() != -13 is rejected.
    return reconstructForLeptons(event, finalState, {-13});
  }

  /// E143 R1998: arithmetic mean of the published Ra, Rb, and Rc fits.
  inline double r1998(double x, double Q2) {
    if (!(x > 0.0) || !(Q2 > 0.04))
      return std::numeric_limits<double>::quiet_NaN();
    const std::array<double, 6> a =
      {{0.0485, 0.5470, 2.0621, -0.3804, 0.5090, -0.0285}};
    const std::array<double, 6> b =
      {{0.0481, 0.6114, -0.3509, -0.4611, 0.7172, -0.0317}};
    const std::array<double, 6> c =
      {{0.0577, 0.4644, 1.8288, 12.3708, -43.1043, 41.7415}};
    const double theta = 1.0 + 12.0 * Q2 / (Q2 + 1.0)
      * std::pow(0.125, 2) / (std::pow(0.125, 2) + x*x);
    const double logarithm = std::log(Q2 / 0.04);
    const double Ra = a[0] / logarithm * theta
      + a[1] / std::pow(std::pow(Q2, 4) + std::pow(a[2], 4), 0.25)
      * (1.0 + a[3]*x + a[4]*x*x) * std::pow(x, a[5]);
    const double Rb = b[0] / logarithm * theta
      + (b[1] / Q2 + b[2] / (Q2*Q2 + std::pow(0.3, 2)))
      * (1.0 + b[3]*x + b[4]*x*x) * std::pow(x, b[5]);
    const double Q2threshold = c[3]*x + c[4]*x*x + c[5]*x*x*x;
    const double Rc = c[0] / logarithm * theta
      + c[1] / std::sqrt(std::pow(Q2 - Q2threshold, 2) + c[2]*c[2]);
    return (Ra + Rb + Rc) / 3.0;
  }

  inline double eta(double x, double y, double Q2) {
    constexpr double nucleonMass = 0.938918754;
    constexpr double muonMass = 0.1056583755;
    const double gamma = 2.0 * nucleonMass * x / std::sqrt(Q2);
    const double massTerm = y*y*muonMass*muonMass/Q2;
    const double numerator = gamma *
      (1.0 - y - 0.25*gamma*gamma*y*y - massTerm);
    const double denominator =
      (1.0 + 0.5*gamma*gamma*y)*(1.0 - 0.5*y) - massTerm;
    return numerator / denominator;
  }

  inline double depolarization(double x, double y, double Q2) {
    constexpr double nucleonMass = 0.938918754;
    constexpr double muonMass = 0.1056583755;
    const double gamma2 = 4.0*nucleonMass*nucleonMass*x*x/Q2;
    const double muon2OverQ2 = muonMass*muonMass/Q2;
    const double numerator = y *
      ((1.0 + 0.5*gamma2*y)*(2.0-y) - 2.0*y*y*muon2OverQ2);
    const double denominator =
      y*y*(1.0 - 2.0*muon2OverQ2)*(1.0 + gamma2)
      + 2.0*(1.0 + r1998(x, Q2))*
        (1.0-y-0.25*gamma2*y*y);
    return numerator / denominator;
  }

  inline void fillMeasurement(const Kinematics& dis, double D,
                              const Histo1DPtr& ordinary,
                              const Histo1DPtr& weighted,
                              const Histo1DPtr& covariance) {
    ordinary->fill(dis.x);
    weighted->fill(dis.x, 1.0/D);
    // sumW2 is sum(eventWeight^2/D): the ordinary--1/D covariance.
    covariance->fill(dis.x, 1.0/std::sqrt(D));
  }

  inline void fillDiagnostics(const Kinematics& dis,
                              const Histo1DPtr& x,
                              const Histo1DPtr& q2,
                              const Histo1DPtr& y,
                              const Histo1DPtr& w2,
                              const Histo1DPtr& theta,
                              const Histo1DPtr& outgoingEnergy) {
    x->fill(dis.x);
    q2->fill(dis.Q2);
    y->fill(dis.y);
    w2->fill(dis.W2);
    theta->fill(dis.theta);
    outgoingEnergy->fill(dis.outgoingEnergy);
  }

} // namespace COMPASSInclusiveDIS
} // namespace Rivet
