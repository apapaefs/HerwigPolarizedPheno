// -*- C++ -*-
#pragma once
#include "SIDISAzimuthal.hh"

namespace Rivet {
namespace SIDISAzimuthal {

  /// Store off-diagonal event second moments in triangular cell order.
  /// Inputs are nonnegative event counts or epsilon-weighted counts. The
  /// proxy's sumW2 stores sum(w_event^2 * input_i * input_j), including NEG.
  inline void fillCellCovariance(const Histo1DPtr& covariance,
                                 const std::vector<double>& inputs,
                                 size_t dimensions) {
    const size_t pairs = dimensions*(dimensions-1)/2;
    for (size_t offset = 0; offset < inputs.size(); offset += dimensions) {
      for (size_t i = 0; i < dimensions; ++i) {
        if (inputs[offset+i] == 0.) continue;
        for (size_t j = i+1; j < dimensions; ++j) {
          const double product = inputs[offset+i]*inputs[offset+j];
          if (product <= 0.) continue;
          const size_t pair = i*(2*dimensions-i-1)/2 + j-i-1;
          covariance->fill(double((offset/dimensions)*pairs + pair)+.5,
                           std::sqrt(product));
        }
      }
    }
  }

} // namespace SIDISAzimuthal
} // namespace Rivet
