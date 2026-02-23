// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

namespace xdiag::combinatorics {

template <typename bit_t> bit_t next_combination(bit_t v) noexcept;

// next_combination(v, n): like next_combination(v) for integral bit_t, but for
// Bitset types uses only test/set/reset (O(1) each) instead of multi-word
// arithmetic, which is O(nchunks) per operation. Prefer this overload from
// iterators that already store n.
template <typename bit_t> bit_t next_combination(bit_t v, int64_t n) noexcept;

template <class bit_t> bit_t nth_combination(int64_t n, int64_t k, int64_t idx);
template <class bit_t> int64_t rank_combination(bit_t bits, int64_t n);

// template <typename bit_t>
// bit_t get_nth_pattern(int64_t n, int64_t nsites, int64_t nupspins);

// template <typename bit_t>
// int64_t get_n_for_pattern(bit_t pattern, int64_t nsites, int64_t nupspins);

} // namespace xdiag::combinatorics
