// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "norm.hpp"

#include <cmath>
#include <cstdint>

#include <xdiag/bits/bitarray.hpp>
#include <xdiag/bits/bitset.hpp>

namespace xdiag::symmetries {

// determines whether a state is a representative
template <typename bit_t, typename coeff_t>
double norm(bit_t state, SitePermutation const &action,
            arma::Col<coeff_t> const &characters) {
  coeff_t amplitude = 0.0;
  for (int64_t sym = 0; sym < action.size(); ++sym) {
    bit_t tstate = action.apply(sym, state);
    if (tstate == state) {
      amplitude += characters(sym);
    }
  }
  return std::sqrt(std::abs(amplitude));
}

#define INSTANTIATE_NORM(BIT_TYPE)                                             \
  template double norm(BIT_TYPE, SitePermutation const &, arma::vec const &);  \
  template double norm(BIT_TYPE, SitePermutation const &,                      \
                       arma::cx_vec const &);                                  \
  using namespace bits;

INSTANTIATE_NORM(uint16_t);
INSTANTIATE_NORM(uint32_t);
INSTANTIATE_NORM(uint64_t);
INSTANTIATE_NORM(BitsetDynamic);
INSTANTIATE_NORM(BitsetStatic1);
INSTANTIATE_NORM(BitsetStatic2);
INSTANTIATE_NORM(BitsetStatic4);
INSTANTIATE_NORM(BitsetStatic8);

INSTANTIATE_NORM(BitArray1);
INSTANTIATE_NORM(BitArray2);
INSTANTIATE_NORM(BitArray3);
INSTANTIATE_NORM(BitArray4);
INSTANTIATE_NORM(BitArray5);
INSTANTIATE_NORM(BitArray6);
INSTANTIATE_NORM(BitArray7);
INSTANTIATE_NORM(BitArray8);

INSTANTIATE_NORM(BitArrayLong1);
INSTANTIATE_NORM(BitArrayLong2);
INSTANTIATE_NORM(BitArrayLong3);
INSTANTIATE_NORM(BitArrayLong4);
INSTANTIATE_NORM(BitArrayLong5);
INSTANTIATE_NORM(BitArrayLong6);
INSTANTIATE_NORM(BitArrayLong7);
INSTANTIATE_NORM(BitArrayLong8);

#undef INSTANTIATE_NORM

} // namespace xdiag::symmetries
