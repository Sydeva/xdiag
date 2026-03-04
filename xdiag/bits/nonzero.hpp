// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <xdiag/bits/bitset.hpp>

namespace xdiag::bits {
constexpr bool nonzero(uint8_t s) { return (bool)s; }
constexpr bool nonzero(uint16_t s) { return (bool)s; }
constexpr bool nonzero(uint32_t s) { return (bool)s; }
constexpr bool nonzero(uint64_t s) { return (bool)s; }
constexpr bool nonzero(int8_t s) { return (bool)s; }
constexpr bool nonzero(int16_t s) { return (bool)s; }
constexpr bool nonzero(int32_t s) { return (bool)s; }
constexpr bool nonzero(int64_t s) { return (bool)s; }

template <typename chunk_t, int64_t n_chunks>
constexpr bool nonzero(Bitset<chunk_t, n_chunks> const &bits) {
  return bits.any();
}

} // namespace xdiag::bits
