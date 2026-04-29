// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <functional>

#include <xdiag/basis/spinhalf/apply/apply_term_offdiag_no_sym.hpp>
#include <xdiag/basis/spinhalf/apply/apply_term_offdiag_sym.hpp>
#include <xdiag/bits/bitops.hpp>
#include <xdiag/common.hpp>

namespace xdiag::basis::spinhalf {

// S+S+ term: J S^+_i S^+_j   (both spins must be down; flips both up)
// S-S- term: J S^-_i S^-_j   (both spins must be up;   flips both down)

template <typename coeff_t, bool symmetric, class basis_t, class fill_f>
void apply_spsp_smsm(Coupling const &cpl, Op const &op, basis_t const &basis_in,
                     basis_t const &basis_out, fill_f fill) {
  using bit_t = typename basis_t::bit_t;

  coeff_t J = cpl.scalar().as<coeff_t>();
  int64_t s1 = op[0];
  int64_t s2 = op[1];
  bit_t mask1 = ((bit_t)1 << s1);
  bit_t mask2 = ((bit_t)1 << s2);
  bit_t flipmask = mask1 | mask2;

  std::function<bool(bit_t)> non_zero_term;
  std::function<std::pair<bit_t, coeff_t>(bit_t)> term_action;

  if (op.type() == "S+S+") {
    // non-zero only when both sites are spin-down (bits 0)
    non_zero_term = [&](bit_t spins) -> bool {
      return !(spins & mask1) && !(spins & mask2);
    };
    term_action = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
      return {spins | flipmask, J};
    };
  } else { // "S-S-"
    // non-zero only when both sites are spin-up (bits 1)
    non_zero_term = [&](bit_t spins) -> bool {
      return (spins & mask1) && (spins & mask2);
    };
    term_action = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
      return {spins & ~flipmask, J};
    };
  }

  if constexpr (symmetric) {
    spinhalf::apply_term_offdiag_sym<coeff_t>(basis_in, basis_out,
                                              non_zero_term, term_action, fill);
  } else {
    spinhalf::apply_term_offdiag_no_sym<coeff_t>(
        basis_in, basis_out, non_zero_term, term_action, fill);
  }
}

} // namespace xdiag::basis::spinhalf

