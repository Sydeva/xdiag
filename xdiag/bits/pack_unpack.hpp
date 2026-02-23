// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <xdiag/bits/bitarray.hpp>

namespace xdiag::bits {

// returns the q-ary represenation of number as a bitarray
template <typename bit_t, int nbits>
BitArray<bit_t, nbits> unpack(int64_t number, int64_t q);

// inverse of unpack: returns the integer whose q-ary representation is the
// first n slots of array
template <typename bit_t, int nbits>
int64_t pack(BitArray<bit_t, nbits> array, int64_t q, int64_t n);

} // namespace xdiag::bits
