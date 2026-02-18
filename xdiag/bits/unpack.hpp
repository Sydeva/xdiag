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

} // namespace xdiag::bits
