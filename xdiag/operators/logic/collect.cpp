// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "collect.hpp"

#include <map>
#include <xdiag/math/scalar.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag {

OpSum collect(OpSum const &ops, double tol) try {
  // Group terms by monomial, summing scalar coefficients
  // Use a vector to preserve insertion order (for deterministic output)
  std::vector<Monomial> order;
  std::map<Monomial, Scalar> sums;

  for (auto const &[coeff, mono] : ops.plain()) {
    if (sums.find(mono) == sums.end()) {
      order.push_back(mono);
      sums[mono] = coeff.scalar();
    } else {
      sums[mono] += coeff.scalar();
    }
  }

  OpSum result;
  for (auto const &mono : order) {
    Scalar c = sums.at(mono);
    if (abs(c) > tol) {
      result += c * mono;
    }
  }
  return result;
}
XDIAG_CATCH

} // namespace xdiag
