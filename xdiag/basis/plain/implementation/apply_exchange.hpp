// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <utility>

#include <xdiag/basis/plain/implementation/apply_offdiag.hpp>
#include <xdiag/bits/get_set_bit.hpp>
#include <xdiag/bits/popcount.hpp>

namespace xdiag::basis::plain {

template <typename coeff_t, class basis_t, class fill_f>
void apply_exchange(Coeff const &c, Op const &op, basis_t const &basis_in,
                    basis_t const &basis_out, fill_f fill) {
  using bit_t = typename basis_t::bit_t;

  coeff_t J = c.scalar().as<coeff_t>();
  bit_t mask = bit_t();
  int64_t s1 = op[0];
  int64_t s2 = op[1];
  bits::set_bit(mask, s1);
  bits::set_bit(mask, s2);

  auto non_zero_term = [&](bit_t spins) -> bool {
    return bits::popcount(spins & mask) & 1;
  };

  coeff_t Jhalf = J / 2.0;
  if constexpr (isreal<coeff_t>()) {
    apply_offdiag(
        basis_in, basis_out, non_zero_term,
        [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
          return {spins ^ mask, Jhalf};
        },
        fill);
  } else {
    coeff_t Jhalf_conj = conj(Jhalf);
    apply_offdiag(
        basis_in, basis_out, non_zero_term,
        [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
          return {spins ^ mask, bits::get_bit(spins, s1) ? Jhalf : Jhalf_conj};
        },
        fill);
  }
}

} // namespace xdiag::basis::plain
