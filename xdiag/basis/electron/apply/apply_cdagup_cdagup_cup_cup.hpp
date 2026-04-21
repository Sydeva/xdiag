// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cassert>
#include <utility>

#include <xdiag/basis/electron/apply/generic_term_ups.hpp>
#include <xdiag/bits/bitops.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/op.hpp>

namespace xdiag::basis::electron {

namespace detail {
template <typename bit_t>
inline bool apply_cup(int64_t site, bit_t &ups, int &sign) {
  bit_t site_mask = ((bit_t)1 << site);
  if ((ups & site_mask) == 0) {
    return false;
  }
  if (bits::popcnt(ups & (site_mask - 1)) & 1) {
    sign = -sign;
  }
  ups ^= site_mask;
  return true;
}

template <typename bit_t>
inline bool apply_cdagup(int64_t site, bit_t &ups, int &sign) {
  bit_t site_mask = ((bit_t)1 << site);
  if (ups & site_mask) {
    return false;
  }
  if (bits::popcnt(ups & (site_mask - 1)) & 1) {
    sign = -sign;
  }
  ups ^= site_mask;
  return true;
}
} // namespace detail

template <bool symmetric, typename coeff_t, typename basis_t, typename fill_f>
void apply_cdagup_cdagup_cup_cup(Coupling const &cpl, Op const &op,
                                 basis_t const &basis_in,
                                 basis_t const &basis_out, fill_f fill) {
  using bit_t = typename basis_t::bit_t;

  coeff_t c = cpl.scalar().as<coeff_t>();
  int64_t i = op[0];
  int64_t j = op[1];
  int64_t k = op[2];
  int64_t l = op[3];

  auto non_zero_term = [&](bit_t ups_in) -> bool {
    bit_t ups = ups_in;
    int sign = 1;
    return detail::apply_cup(i, ups, sign) && detail::apply_cup(j, ups, sign) &&
           detail::apply_cdagup(k, ups, sign) &&
           detail::apply_cdagup(l, ups, sign);
  };

  auto term_action = [&](bit_t ups_in) -> std::pair<bit_t, coeff_t> {
    bit_t ups = ups_in;
    int sign = 1;
    [[maybe_unused]] bool ok =
        detail::apply_cup(i, ups, sign) && detail::apply_cup(j, ups, sign) &&
        detail::apply_cdagup(k, ups, sign) &&
        detail::apply_cdagup(l, ups, sign);
    assert(ok);
    return {ups, sign > 0 ? c : -c};
  };

  generic_term_ups<symmetric, coeff_t>(basis_in, basis_out, non_zero_term,
                                       term_action, fill);
}

template <bool symmetric, typename coeff_t, typename basis_t, typename fill_f>
void apply_cdagup_cdagup_cup_cup_hc(Coupling const &cpl, Op const &op,
                                    basis_t const &basis_in,
                                    basis_t const &basis_out, fill_f fill) {
  using bit_t = typename basis_t::bit_t;

  coeff_t c = cpl.scalar().as<coeff_t>();
  int64_t i = op[0];
  int64_t j = op[1];
  int64_t k = op[2];
  int64_t l = op[3];

  auto non_zero_term = [&](bit_t ups_in) -> bool {
    bit_t ups = ups_in;
    int sign = 1;
    return detail::apply_cup(l, ups, sign) && detail::apply_cup(k, ups, sign) &&
           detail::apply_cdagup(j, ups, sign) &&
           detail::apply_cdagup(i, ups, sign);
  };

  auto term_action = [&](bit_t ups_in) -> std::pair<bit_t, coeff_t> {
    bit_t ups = ups_in;
    int sign = 1;
    [[maybe_unused]] bool ok =
        detail::apply_cup(l, ups, sign) && detail::apply_cup(k, ups, sign) &&
        detail::apply_cdagup(j, ups, sign) &&
        detail::apply_cdagup(i, ups, sign);
    assert(ok);
    return {ups, sign > 0 ? c : -c};
  };

  generic_term_ups<symmetric, coeff_t>(basis_in, basis_out, non_zero_term,
                                       term_action, fill);
}

} // namespace xdiag::basis::electron



