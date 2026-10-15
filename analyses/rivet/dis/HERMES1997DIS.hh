// -*- C++ -*-
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

// HERMES, Phys. Lett. B404 (1997) 383, Table I (hep-ex/9703005).
// These windows are a theory-comparison prescription, not detector acceptance.
namespace HERMES1997DIS {
  constexpr std::array<double, 10> xEdges =
    {{0.023, 0.04, 0.055, 0.075, 0.1, 0.14, 0.2, 0.3, 0.4, 0.6}};
  constexpr std::array<double, 9> q2Means =
    {{1.22, 1.47, 1.73, 1.99, 2.3, 2.65, 3.07, 3.86, 5.25}};
  constexpr std::array<double, 3> windowFractions = {{0.025, 0.05, 0.10}};
  constexpr std::array<const char*, 3> windowNames = {{"W025", "W050", "W100"}};

  inline bool passesFiducial(double q2, double w2, double y, double theta) {
    // Paper p. 3: Q2 > 1, W2 > 4, y < .85, 40 < theta < 220 mrad.
    // Q2 < 15 is the published outer range; all windows lie below 5.775.
    return std::isfinite(q2) && std::isfinite(w2) && std::isfinite(y) && std::isfinite(theta)
      && q2 > 1.0 && q2 < 15.0 && w2 > 4.0
      && y > 0.0 && y < 0.85 && theta > 0.04 && theta < 0.22;
  }

  inline int xBin(double x) {
    // Half-open bins agree with YODA, including exclusion of the final edge.
    if (!std::isfinite(x) || x < xEdges.front() || x >= xEdges.back()) return -1;
    return int(std::upper_bound(xEdges.begin(), xEdges.end(), x) - xEdges.begin()) - 1;
  }

  inline bool inWindow(double x, double q2, double fraction) {
    const int bin = xBin(x);
    if (bin < 0 || !std::isfinite(q2) || q2 <= 1.0 ||
        !std::isfinite(fraction) || fraction <= 0.0) return false;
    const double centre = q2Means[bin];
    return q2 >= centre*(1.0-fraction) && q2 <= centre*(1.0+fraction);
  }

  inline double r1990(double x, double q2) {
    // Whitlow R1990 fit B, as used by the inclusive HERMES template.
    // The HERMES 1997 publication cites Whitlow (1990), ref. 32.
    // This implementation uses fit B, consistently with the existing HERMES template.
    if (!(x > 0.0) || !(q2 > 0.04)) return std::numeric_limits<double>::quiet_NaN();
    const double s = 0.125*0.125;
    const double theta = 1.0 + 12.0*q2/(q2+1.0)*s/(s+x*x);
    return 0.0635/std::log(q2/0.04)*theta + 0.5747/q2 - 0.3534/(q2*q2+0.09);
  }

  inline double depolarization(double x, double y, double q2, double mass) {
    if (!(x > 0.0) || !(y > 0.0 && y < 1.0) || !(q2 > 0.04) || !(mass > 0.0))
      return std::numeric_limits<double>::quiet_NaN();
    const double gamma2 = 4.0*mass*mass*x*x/q2;
    const double numerator = y*(2.0-y)*(1.0+0.5*gamma2*y);
    const double denominator = y*y*(1.0+gamma2)
      + 2.0*(1.0+r1990(x,q2))*(1.0-y-0.25*gamma2*y*y);
    return numerator/denominator;
  }
}
