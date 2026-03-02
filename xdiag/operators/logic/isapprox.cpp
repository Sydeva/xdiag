// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "isapprox.hpp"

#include <xdiag/math/matrix.hpp>
#include <xdiag/operators/logic/order.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag {

bool isapprox(Op const &op1, Op const &op2, double rtol, double atol) try {
  auto [a1, o1] = order(op1);
  auto [a2, o2] = order(op2);

  if (o1.type() == o2.type()) {
    if (o1.hassites() && o2.hassites()) {
      auto s1 = o1.sites();
      auto s2 = o2.sites();
      if (s1 == s2) {
        if (o1.hasmatrix() && o2.hasmatrix()) {
          return isapprox(o1.matrix() * a1, o2.matrix() * a2, rtol, atol);
        } else if (!o1.hasmatrix() && !o2.hasmatrix()) {
          return isapprox(a1, a2, rtol, atol);
        } else {
          return false;
        }
      } else {
        return false;
      }
    } else if (!o1.hassites() && !o2.hassites()) {
      return isapprox(a1, a2, rtol, atol);
    } else {
      return false;
    }
  } else {
    return false;
  }
}
XDIAG_CATCH

bool isapprox(OpSum const &ops1, OpSum const &ops2, double rtol,
              double atol) try {
  auto t1z = order(ops1).terms();
  auto t2z = order(ops2).terms();

  // Remove zero entries
  OpSum t1o;
  for (auto const &[coeff, mono] : t1z) {
    if (!isapprox(coeff.scalar(), Scalar(0.), rtol, atol)) {
      t1o += coeff * mono;
    }
  }
  OpSum t2o;
  for (auto const &[coeff, mono] : t2z) {
    if (!isapprox(coeff.scalar(), Scalar(0.), rtol, atol)) {
      t2o += coeff * mono;
    }
  }
  auto const &t1 = t1o.terms();
  auto const &t2 = t2o.terms();

  if (t1.size() != t2.size()) {
    return false;
  }
  for (int64_t i = 0; i < (int64_t)t1.size(); ++i) {
    Scalar a1 = t1[i].coeff.scalar();
    Scalar a2 = t2[i].coeff.scalar();
    // For length-1 monomials compare the Op; for general monomials compare Op
    // by Op
    if (t1[i].monomial.size() != t2[i].monomial.size()) {
      return false;
    }
    if (!isapprox(a1, a2, rtol, atol)) {
      return false;
    }
    for (int64_t j = 0; j < t1[i].monomial.size(); ++j) {
      if (!isapprox(t1[i].monomial[j], t2[i].monomial[j], rtol, atol)) {
        return false;
      }
    }
  }
  return true;
}
XDIAG_CATCH

std::optional<Scalar> isapprox_multiple(OpSum const &ops1, OpSum const &ops2,
                                        double rtol, double atol) try {
  auto const &t1 = order(ops1).terms();
  auto const &t2 = order(ops2).terms();

  if ((t1.size() != t2.size()) || (t1.size() == 0)) {
    return std::nullopt;
  }

  // Determine first non-zero coefficient of ops1 and ratio of the terms
  Scalar ratio = 0.;
  int64_t idx0 = 0;
  for (; idx0 < (int64_t)t1.size(); ++idx0) {
    Scalar a01 = t1[idx0].coeff.scalar();
    if (abs(a01) > 1e-12) {
      Scalar a02 = t2[idx0].coeff.scalar();
      ratio = a02 / a01;
      break;
    }
  }
  if (idx0 == (int64_t)t1.size()) {
    XDIAG_THROW("All coefficients in first operator are zero.");
  }

  for (int64_t i = 0; i < (int64_t)t1.size(); ++i) {
    Scalar a1 = t1[i].coeff.scalar();
    Scalar a2 = t2[i].coeff.scalar();
    if (!isapprox(a1 * ratio, a2, rtol, atol)) {
      return std::nullopt;
    }
    if (t1[i].monomial.size() != t2[i].monomial.size()) {
      return std::nullopt;
    }
    for (int64_t j = 0; j < t1[i].monomial.size(); ++j) {
      if (!isapprox(t1[i].monomial[j], t2[i].monomial[j], rtol, atol)) {
        return std::nullopt;
      }
    }
  }
  return ratio;
}
XDIAG_CATCH

} // namespace xdiag
