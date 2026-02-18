// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "bit_patterns.hpp"

#include <type_traits>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/combinatorics/binomial.hpp>

namespace xdiag::combinatorics {

// Helper functions for generic bit operations (zero overhead for native types)
template <typename bit_t>
static constexpr bit_t bit_zero(int64_t nbits = 64) noexcept {
  if constexpr (std::is_integral<bit_t>::value) {
    return 0;
  } else {
    return bit_t(nbits);
  }
}

template <typename bit_t>
static constexpr bit_t bit_one(int64_t nbits = 64) noexcept {
  if constexpr (std::is_integral<bit_t>::value) {
    return 1;
  } else {
    bit_t result(nbits);
    result.set(0);
    return result;
  }
}

template <typename bit_t>
bit_t get_next_pattern(bit_t v) noexcept {
  // Bit twiddling Hack from
  // http://graphics.stanford.edu/~seander/bithacks.html
  // #NextBitPermutation

  // DONT USE FAST VERSION

  // // Fast version (needs __builtin_ctz(v)) (some problem with 0)
  // bit_t t = v | (v - 1); // t gets v's least significant 0
  // return ((t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(v) + 1)));

  // // Fast version (needs __builtin_ctz(v)) (some problem with 0)
  // int_t t = v | (v - 1); // t gets v's least significant 0
  // return v == 0 ? ~v :((t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(v) +
  // 1)));

  // Slow version (should work everywhere)
  const bit_t zero = bit_zero<bit_t>();
  const bit_t one = bit_one<bit_t>();
  bit_t t = (v | (v - one)) + one;
  return v == zero ? ~v : t | ((((t & -t) / (v & -v)) >> 1) - one);
}

template <typename bit_t>
bit_t get_nth_pattern(int64_t n, int64_t nsites, int64_t nupspins) {
  bit_t state = bit_zero<bit_t>(nsites + 1);
  int64_t counter = n;
  for (int64_t n_varying_bits = nupspins - 1; n_varying_bits >= 0;
       --n_varying_bits) {
    int64_t n_combinations = 0;
    for (int64_t n_allowed_pos = n_varying_bits; n_allowed_pos <= nsites;
         ++n_allowed_pos) {
      n_combinations += binomial(n_allowed_pos, n_varying_bits);

      if (n_combinations > counter) {
        counter -= n_combinations - binomial(n_allowed_pos, n_varying_bits);
        state |= (bit_one<bit_t>(nsites + 1) << n_allowed_pos);
        break;
      }
    }
  }
  return state;
}

template <typename bit_t>
int64_t get_n_for_pattern(bit_t pattern, int64_t nsites, int64_t nupspins) {
  int64_t n = 0;
  bit_t workpattern = pattern;
  for (int64_t n_varying_bits = nupspins - 1; n_varying_bits >= 0;
       --n_varying_bits) {
    for (int64_t i = 0; i <= nsites; ++i) {
      // MSB is at 2^i
      if ((bit_one<bit_t>() << (i + 1)) > workpattern) {
        n += binomial(i, n_varying_bits + 1);
        workpattern ^= (bit_one<bit_t>() << i);
        break;
      }
    }
  }
  return n;
}

template uint16_t get_next_pattern<uint16_t>(uint16_t v) noexcept;
template uint32_t get_next_pattern<uint32_t>(uint32_t v) noexcept;
template uint64_t get_next_pattern<uint64_t>(uint64_t v) noexcept;

template uint16_t get_nth_pattern<uint16_t>(int64_t n, int64_t nsites,
                                            int64_t nupspins);
template uint32_t get_nth_pattern<uint32_t>(int64_t n, int64_t nsites,
                                            int64_t nupspins);
template uint64_t get_nth_pattern<uint64_t>(int64_t n, int64_t nsites,
                                            int64_t nupspins);

template int64_t get_n_for_pattern<uint16_t>(uint16_t pattern, int64_t nsites,
                                             int64_t nupspins);
template int64_t get_n_for_pattern<uint32_t>(uint32_t pattern, int64_t nsites,
                                             int64_t nupspins);
template int64_t get_n_for_pattern<uint64_t>(uint64_t pattern, int64_t nsites,
                                             int64_t nupspins);

// Bitset instantiations
#define INSTANTIATE_BIT_PATTERNS(CHUNK_T, NCHUNKS)                            \
  template bits::Bitset<CHUNK_T, NCHUNKS>                                     \
  get_next_pattern<bits::Bitset<CHUNK_T, NCHUNKS>>(                           \
      bits::Bitset<CHUNK_T, NCHUNKS> v) noexcept;                             \
  template bits::Bitset<CHUNK_T, NCHUNKS>                                     \
  get_nth_pattern<bits::Bitset<CHUNK_T, NCHUNKS>>(int64_t n, int64_t nsites, \
                                                   int64_t nupspins);          \
  template int64_t get_n_for_pattern<bits::Bitset<CHUNK_T, NCHUNKS>>(        \
      bits::Bitset<CHUNK_T, NCHUNKS> pattern, int64_t nsites,                 \
      int64_t nupspins);

#define INSTANTIATE_BIT_PATTERNS_FOR_NCHUNKS(CHUNK_T)                         \
  INSTANTIATE_BIT_PATTERNS(CHUNK_T, 0)                                         \
  INSTANTIATE_BIT_PATTERNS(CHUNK_T, 1)                                         \
  INSTANTIATE_BIT_PATTERNS(CHUNK_T, 2)                                         \
  INSTANTIATE_BIT_PATTERNS(CHUNK_T, 4)                                         \
  INSTANTIATE_BIT_PATTERNS(CHUNK_T, 8)

INSTANTIATE_BIT_PATTERNS_FOR_NCHUNKS(uint8_t)
INSTANTIATE_BIT_PATTERNS_FOR_NCHUNKS(uint16_t)
INSTANTIATE_BIT_PATTERNS_FOR_NCHUNKS(uint32_t)
INSTANTIATE_BIT_PATTERNS_FOR_NCHUNKS(uint64_t)

#undef INSTANTIATE_BIT_PATTERNS_FOR_NCHUNKS
#undef INSTANTIATE_BIT_PATTERNS

} // namespace xdiag::combinatorics
