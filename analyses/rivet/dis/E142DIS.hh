// -*- C++ -*-
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

// E142, Phys. Rev. D54 (1996) 6620, Tables 7 and 9 and Eqs. 39--47.
// These windows are a theory-comparison prescription, not detector acceptance.
namespace E142DIS {
  constexpr std::array<double, 9> xEdges =
    {{0.03, 0.04, 0.06, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60}};
  constexpr std::array<double, 8> q2Means =
    {{1.1, 1.2, 1.8, 2.5, 3.1, 3.7, 4.4, 5.5}};
  constexpr std::array<double, 3> windowFractions = {{0.025, 0.05, 0.10}};
  constexpr std::array<const char*, 3> windowNames = {{"W025", "W050", "W100"}};

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
    // E142 uses the SLAC global R fit (paper ref. 71); this fit is an
    // explicit choice in the approximate prediction's extraction convention.
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
