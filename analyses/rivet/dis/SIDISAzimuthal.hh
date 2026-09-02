// -*- C++ -*-
#pragma once

#include "COMPASSSIDIS.hh"
#include <algorithm>
#include <cmath>
#include <vector>

namespace Rivet {
namespace SIDISAzimuthal {

  struct Angle {
    bool valid = false;
    double cosine = 0.0;
    double sine = 0.0;
  };

  /// Trento-style hadron azimuth about q, measured from the lepton plane.
  ///
  /// The present analyses use only cosine harmonics, so the result is
  /// insensitive to the overall sign convention for phi.  Returning the
  /// signed sine as well makes the convention explicit and testable.
  inline Angle hadronAzimuth(
      const COMPASSSIDIS::DISKinematics& dis, const Particle& particle) {
    Angle out;
    if (dis.q.p3().mod() == 0.0) return out;
    const Vector3 qhat = dis.q.p3().unit();
    const Vector3 leptonTransverse =
      dis.k.p3() - qhat*dis.k.p3().dot(qhat);
    const Vector3 hadronTransverse =
      particle.momentum().p3()
      - qhat*particle.momentum().p3().dot(qhat);
    const double denominator =
      leptonTransverse.mod()*hadronTransverse.mod();
    if (!(denominator > 0.0)) return out;
    out.cosine = std::max(-1.0, std::min(
      1.0, leptonTransverse.dot(hadronTransverse)/denominator));
    out.sine = std::max(-1.0, std::min(
      1.0, qhat.dot(leptonTransverse.cross(hadronTransverse))/denominator));
    out.valid = std::isfinite(out.cosine) && std::isfinite(out.sine);
    return out;
  }

  inline double cosineHarmonic(const Angle& angle, int harmonic) {
    if (harmonic == 1) return angle.cosine;
    if (harmonic == 2) return 2.0*angle.cosine*angle.cosine - 1.0;
    return std::cos(double(harmonic)*std::atan2(angle.sine, angle.cosine));
  }

  /// COMPASS depolarization factors multiplying A_UU^{cos(n phi)}.
  inline double epsilon(double y, int harmonic) {
    const double denominator = 1.0 + (1.0-y)*(1.0-y);
    if (!(y >= 0.0 && y <= 1.0) || !(denominator > 0.0)) return 0.0;
    if (harmonic == 1)
      return 2.0*(2.0-y)*std::sqrt(std::max(0.0, 1.0-y))/denominator;
    if (harmonic == 2)
      return 2.0*(1.0-y)/denominator;
    return 0.0;
  }

  inline int binIndex(double value, const std::vector<double>& edges) {
    if (edges.size() < 2 || value < edges.front() || value >= edges.back())
      return -1;
    return int(std::upper_bound(edges.begin(), edges.end(), value)
               - edges.begin()) - 1;
  }

  inline double feynmanX(
      const COMPASSSIDIS::DISKinematics& dis, const Particle& particle) {
    const FourMomentum hadronicSystem = dis.target + dis.q;
    if (!(hadronicSystem.mass() > 0.0)) return -1.0;
    LorentzTransform boost;
    boost.setBetaVec(-hadronicSystem.betaVec());
    const FourMomentum qcm = boost.transform(dis.q);
    const FourMomentum hcm = boost.transform(particle.momentum());
    if (qcm.p3().mod() == 0.0) return -1.0;
    return 2.0*hcm.p3().dot(qcm.p3().unit())/hadronicSystem.mass();
  }

  inline void fillEventBins(
      const Histo1DPtr& histogram, const std::vector<double>& values) {
    for (size_t index = 0; index < values.size(); ++index)
      if (values[index] != 0.0)
        histogram->fill(double(index) + 0.5, values[index]);
  }

  inline void fillSignedCovariance(
      const Histo1DPtr& positive, const Histo1DPtr& negative,
      const std::vector<double>& numerator,
      const std::vector<double>& denominator) {
    for (size_t index = 0; index < numerator.size(); ++index) {
      const double covariance = numerator[index]*denominator[index];
      if (covariance > 0.0)
        positive->fill(double(index) + 0.5, std::sqrt(covariance));
      else if (covariance < 0.0)
        negative->fill(double(index) + 0.5, std::sqrt(-covariance));
    }
  }

} // namespace SIDISAzimuthal
} // namespace Rivet
