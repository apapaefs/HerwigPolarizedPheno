// -*- C++ -*-
#pragma once
#include "COMPASSSIDIS.hh"

namespace Rivet {
namespace COMPASSSIDIS {

  constexpr double pionMassGeV = 0.13957039;
  constexpr double kaonMassGeV = 0.493677;

  /// Generated-level multiplicity domain, COMPASS 2025 Sec. 3.1 (also 2017).
  /// The entire z bin must lie in the 12--40 GeV hadron momentum interval.
  /// RICH entrance geometry is part of the acceptance correction, not a
  /// vertex-angle selection on generated hadrons.
  template <typename Cell>
  inline bool multiplicityDISCell(const Cell& cell, const DISKinematics& dis,
                                  double massGeV) {
    if (!(dis.x >= cell.xLow && dis.x < cell.xHigh &&
          dis.y >= cell.yLow && dis.y < cell.yHigh) ||
        !(cell.zLow > 0. && cell.zHigh > cell.zLow) ||
        !(dis.target.mass() > 0.)) return false;
    const double nu = (dis.target*dis.q)/dis.target.mass()/GeV;
    return std::isfinite(nu) &&
      nu > std::sqrt(12.*12. + massGeV*massGeV)/cell.zLow &&
      nu < std::sqrt(40.*40. + massGeV*massGeV)/cell.zHigh;
  }

} // namespace COMPASSSIDIS
} // namespace Rivet
