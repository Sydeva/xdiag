// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/basis/non_branching_op.hpp>
#include <xdiag/basis/plain/implementation/apply_diag.hpp>
#include <xdiag/basis/plain/implementation/apply_offdiag.hpp>

namespace xdiag::basis::plain {

template <typename coeff_t, class basis_t, class fill_f>
void apply_matrix(Coeff const &c, Op const &op, basis_t const &basis_in,
                  basis_t const &basis_out, fill_f fill) {
  using bit_t = typename basis_t::bit_t;

  // Decompose into sum of non-branching operators
  auto ops_nb = basis::non_branching_ops<bit_t, coeff_t>(c, op);

  // Loop over sum of non-branching operators
  for (auto const &op_nb : ops_nb) {

    // Diagonal terms
    if (op_nb.isdiagonal()) {
      auto term_coeff = [&](bit_t spins) -> coeff_t {
        auto local_spins = op_nb.extract(spins);
        return op_nb.coeff(local_spins);
      };
      apply_diag(basis_in, term_coeff, fill);
    } else { // Offdiagonal terms
      auto non_zero_term = [&](bit_t spins) -> bool {
        auto local_spins = op_nb.extract(spins);
        return op_nb.non_zero_term(local_spins);
      };
      auto term_action = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
        auto local_spins = op_nb.extract(spins);
        auto [local_spins_new, coeff] = op_nb.state_coeff(local_spins);
        auto spins_new = op_nb.deposit(local_spins_new, spins);
        return {spins_new, coeff};
      };
      apply_offdiag(basis_in, basis_out, non_zero_term, term_action, fill);
    }
  } // for (auto const &op_nb : ops_nb)
}

} // namespace xdiag::basis::plain
