// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <string>

#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>
#include <xdiag/utils/type_name.hpp>

#include <xdiag/basis/apply_identity.hpp>
#include <xdiag/basis/plain/implementation/apply_exchange.hpp>
#include <xdiag/basis/plain/implementation/apply_matrix.hpp>
#include <xdiag/basis/plain/implementation/apply_scalar_chirality.hpp>
#include <xdiag/basis/plain/implementation/apply_spsm.hpp>
#include <xdiag/basis/plain/implementation/apply_sz.hpp>
#include <xdiag/basis/plain/implementation/apply_szsz.hpp>

namespace xdiag::basis::plain {

template <typename coeff_t, typename basis_t, typename fill_f>
void apply_generic(OpSum const &ops, basis_t const &basis_in,
                   basis_t const &basis_out, fill_f fill) try {

  for (auto const &[c, monomial] : ops) {
    assert(monomial.size() == 1);

    Op op = monomial[0];
    std::string type = op.type();

    if (type == "Id") {
      basis::apply_identity<coeff_t>(c, basis_in, fill);
    } else if (type == "Exchange") {
      plain::apply_exchange<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "SzSz") {
      plain::apply_szsz<coeff_t>(c, op, basis_in, fill);
    } else if (type == "Sz") {
      plain::apply_sz<coeff_t>(c, op, basis_in, fill);
    } else if (type == "S+") {
      plain::apply_spsm<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "S-") {
      plain::apply_spsm<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "ScalarChirality") {
      plain::apply_scalar_chirality<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "Matrix") {
      plain::apply_matrix<coeff_t>(c, op, basis_in, basis_out, fill);
    } else {
      XDIAG_THROW(fmt::format("Unknown Op type for plain basis \"{}\"", type));
    }
  }
}
XDIAG_CATCH

} // namespace xdiag::basis::plain
