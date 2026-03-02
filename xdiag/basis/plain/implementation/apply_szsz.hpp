// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/basis/plain/implementation/apply_diag.hpp>
#include <xdiag/bits/get_set_bit.hpp>
#include <xdiag/bits/popcount.hpp>

namespace xdiag::basis::plain {

template <typename coeff_t, class basis_t, class fill_f>
void apply_szsz(Coeff const &c, Op const &op, basis_t const &basis,
                fill_f fill) {
  assert(Op.size() == 2);
  using bit_t = typename basis_t::bit_t;

  coeff_t J = c.scalar().as<coeff_t>();
  coeff_t val_same = J / 4.0;
  coeff_t val_diff = -J / 4.0;

  bit_t mask = bit_t();
  bits::set_bit(mask, op[0]);
  bits::set_bit(mask, op[1]);

  apply_diag(
      basis,
      [&](bit_t spins) {
        return bits::popcount(spins & mask) & 1 ? val_diff : val_same;
      },
      fill);
}

} // namespace xdiag::basis::plain
