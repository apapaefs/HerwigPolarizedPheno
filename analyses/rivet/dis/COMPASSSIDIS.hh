// -*- C++ -*-
#pragma once

#include "COMPASSInclusiveDIS.hh"
#include "COMPASSSIDISBinning.hh"
#include "Rivet/Projections/FinalState.hh"
#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>

namespace Rivet {
namespace COMPASSSIDIS {

  using DISKinematics = COMPASSInclusiveDIS::Kinematics;

  struct HadronKinematics {
    bool valid = false;
    double z = -1.0;
    double momentum = -1.0;
    double theta = -1.0;
    double transverseMomentum = -1.0;
    double transverseMomentum2 = -1.0;
  };

  inline HadronKinematics hadronKinematics(
      const DISKinematics& dis, const Particle& particle,
      bool pionMassAssumption = false) {
    HadronKinematics out;
    const FourMomentum observed = particle.momentum();
    FourMomentum momentum = observed;
    if (pionMassAssumption) {
      constexpr double pionMass = 0.13957039*GeV;
      const double energy = std::sqrt(observed.p3().mod2() + pionMass*pionMass);
      momentum = FourMomentum(
        energy, observed.px(), observed.py(), observed.pz());
    }
    const double denominator = dis.target * dis.q;
    if (!(denominator > 0.0) || dis.k.p3().mod() == 0.0) return out;
    out.z = (dis.target * momentum) / denominator;
    out.momentum = observed.p3().mod() / GeV;
    out.theta = dis.k.angle(observed);
    if (dis.q.p3().mod() > 0.0) {
      out.transverseMomentum =
        observed.p3().cross(dis.q.p3()).mod()/dis.q.p3().mod()/GeV;
      out.transverseMomentum2 = out.transverseMomentum*out.transverseMomentum;
    }
    out.valid = std::isfinite(out.z) && std::isfinite(out.momentum) &&
                std::isfinite(out.theta) &&
                std::isfinite(out.transverseMomentum) &&
                std::isfinite(out.transverseMomentum2);
    return out;
  }

  inline bool fixed160GeVBeam(const DISKinematics& dis) {
    return std::abs(dis.k.E()/GeV - 160.0) < 1.0e-6;
  }

  inline bool fixedBeamEnergy(const DISKinematics& dis, double energyGeV) {
    return std::abs(dis.k.E()/GeV - energyGeV) < 1.0e-6;
  }

  inline void fillEventCount(
      const Histo1DPtr& histogram, double coordinate, double count) {
    if (count > 0.0) histogram->fill(coordinate, count);
  }

  inline void fillMultiplicityCovariance(
      const Histo1DPtr& histogram, double coordinate, double count) {
    // sumW2 of sqrt(n) is the same-event numerator--DIS covariance n.
    if (count > 0.0) histogram->fill(coordinate, std::sqrt(count));
  }

  inline void fillCountCovariance(
      const Histo1DPtr& histogram, double coordinate,
      double firstCount, double secondCount) {
    // The squared fill weight stores the same-event product needed for the
    // covariance of two event-aggregated identified-hadron yields.
    const double product = firstCount*secondCount;
    if (product > 0.0) histogram->fill(coordinate, std::sqrt(product));
  }

} // namespace COMPASSSIDIS
} // namespace Rivet
