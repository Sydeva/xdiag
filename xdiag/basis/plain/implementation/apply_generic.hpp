// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <string>

#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/type_name.hpp>

#include <xdiag/basis/plain/implementation/apply_exchange.hpp>
#include <xdiag/basis/plain/implementation/apply_szsz.hpp>

namespace xdiag::basis::plain {

template <typename coeff_t, typename basis_t, typename fill_f>
void apply_generic(OpSum const &ops, basis_t const &basis_in,
                   basis_t const &basis_out, fill_f fill) try {
  for (auto const &[c, monomial] : ops) {
    assert(monomial.size() == 1);

    Op op = monomial[0];
    std::string type = op.type();
    if (type == "Exchange") {
      plain::apply_exchange<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "SzSz") {
      plain::apply_szsz<coeff_t>(c, op, basis_in, fill);
    }
  }
}
XDIAG_CATCH

} // namespace xdiag::basis::plain
