// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "hc.hpp"

#include <xdiag/operators/logic/isapprox.hpp>
#include <xdiag/operators/logic/types.hpp>
#include <xdiag/operators/logic/valid.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag {

Op hc(Op const &op) try {
  check_valid(op);

  std::string type = op.type();
  if (type == "S+") {
    return Op("S-", op.sites());
  } else if (type == "S-") {
    return Op("S+", op.sites());
  } else if (type == "Cdagup") {
    return Op("Cup", op.sites());
  } else if (type == "Cup") {
    return Op("Cdagup", op.sites());
  } else if (type == "Cdagdn") {
    return Op("Cdn", op.sites());
  } else if (type == "Cdn") {
    return Op("Cdagdn", op.sites());
  } else { // default: the type does not change
    if (op.hassites()) {
      if (op.hasmatrix()) {
        return Op(op.type(), op.sites(), op.matrix().hc());
      } else {
        return Op(op.type(), op.sites());
      }
    } else {
      return Op(op.type());
    }
  }
}
XDIAG_CATCH

OpSum hc(OpSum const &ops) try {
  OpSum ops_hc;
  for (auto const &[coeff, mono] : ops.plain()) {
    // For length-1 monomials of real-coupling types, coefficient stays real
    if (mono.size() == 1) {
      std::string type = mono[0].type();
      if ((type == "Exchange") || (type == "Hop") || (type == "Hopup") ||
          (type == "Hopdn")) {
        ops_hc += coeff.scalar() * mono.hc();
      } else {
        ops_hc += conj(coeff.scalar()) * mono.hc();
      }
    } else {
      // General monomial: conjugate the coefficient
      ops_hc += conj(coeff.scalar()) * mono.hc();
    }
  }
  return ops_hc;
}
XDIAG_CATCH

bool ishermitian(Op const &op) { return isapprox(op, hc(op)); }
bool ishermitian(OpSum const &ops) { return isapprox(ops, hc(ops)); }

} // namespace xdiag
