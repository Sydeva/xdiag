// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../catch.hpp"

#include <xdiag/bits/bitset.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/utils/logger.hpp>

template <typename bit_t> void test_combinations() {
  using namespace xdiag;
  using namespace xdiag::combinatorics;
  using namespace xdiag::bits;

  for (int n = 0; n < 7; ++n) {
    for (int k = 0; k <= n; ++k) {
      Combinations<bit_t> combs(n, k);
      REQUIRE(n == combs.n());
      REQUIRE(k == combs.k());
      int64_t ctr = 0;
      bit_t current = 0;
      for (auto comb : combs) {

        if (ctr != 0)
          REQUIRE(comb > current);
        current = comb;
        ++ctr;
        REQUIRE(popcnt(comb) == k);
        REQUIRE(comb < ((bit_t)1 << n));
      }
      REQUIRE(ctr == combs.size());
    }
  }
}

// Test Bitset combinations against native uint64_t reference
template <typename chunk_t, int64_t nchunks> void test_combinations_bitset() {
  using namespace xdiag;
  using namespace xdiag::combinatorics;
  using namespace xdiag::bits;

  // Calculate maximum n based on Bitset capacity
  constexpr int64_t chunk_bits = std::numeric_limits<chunk_t>::digits;
  constexpr int64_t max_bits = (nchunks == 0) ? 64 : (nchunks * chunk_bits);
  constexpr int max_n = std::min(int64_t(16), max_bits - 1);

  // Test various (n, k) pairs
  for (int n = 0; n <= max_n; ++n) {
    for (int k = 0; k <= n; ++k) {
      Combinations<Bitset<chunk_t, nchunks>> combs_bitset(n, k);
      Combinations<uint64_t> combs_uint64(n, k);

      REQUIRE(combs_bitset.n() == combs_uint64.n());
      REQUIRE(combs_bitset.k() == combs_uint64.k());
      REQUIRE(combs_bitset.size() == combs_uint64.size());

      // Compare each generated pattern
      auto it_bitset = combs_bitset.begin();
      auto it_uint64 = combs_uint64.begin();

      int64_t count = 0;
      while (it_bitset != combs_bitset.end() &&
             it_uint64 != combs_uint64.end()) {
        auto pattern_bitset = *it_bitset;
        auto pattern_uint64 = *it_uint64;

        // Convert Bitset to uint64_t for comparison
        uint64_t bitset_as_uint64 = to_uint64(pattern_bitset);
        REQUIRE(bitset_as_uint64 == pattern_uint64);

        // Verify bit count
        REQUIRE(pattern_bitset.count() == k);
        REQUIRE(popcnt(pattern_uint64) == k);

        ++it_bitset;
        ++it_uint64;
        ++count;
      }

      REQUIRE(count == combs_bitset.size());
      REQUIRE(it_bitset == combs_bitset.end());
      REQUIRE(it_uint64 == combs_uint64.end());
    }
  }
}

// Test Bitset with larger n values
template <typename chunk_t, int64_t nchunks>
void test_combinations_bitset_large() {
  using namespace xdiag;
  using namespace xdiag::combinatorics;
  using namespace xdiag::bits;

  // Test cases that exceed single uint64_t capacity (n > 64)
  if constexpr (nchunks >= 2) {
    for (int n : {65, 70, 80, 100}) {
      for (int k : {0, 1, 2, n - 1, n}) {
        if (k > n)
          continue;

        Combinations<Bitset<chunk_t, nchunks>> combs(n, k);

        REQUIRE(combs.n() == n);
        REQUIRE(combs.k() == k);

        int64_t count = 0;
        Bitset<chunk_t, nchunks> prev(n);

        for (auto pattern : combs) {
          // Verify bit count
          REQUIRE(pattern.count() == k);

          // Verify patterns are increasing
          if (count > 0) {
            REQUIRE(pattern > prev);
          }

          prev = pattern;
          ++count;
        }

        REQUIRE(count == combs.size());
      }
    }
  }
}

TEST_CASE("Combinations", "[combinatorics]") {
  using namespace xdiag::bits;

  xdiag::Log("Testing Combinations");

  SECTION("native integers") {
    test_combinations<uint16_t>();
    test_combinations<uint32_t>();
    test_combinations<uint64_t>();
  }

  SECTION("bitset vs uint64_t") {
    xdiag::Log("Testing Bitset<uint8_t, 0> combinations");
    test_combinations_bitset<uint8_t, 0>();

    xdiag::Log("Testing Bitset<uint8_t, 1> combinations");
    test_combinations_bitset<uint8_t, 1>();

    xdiag::Log("Testing Bitset<uint16_t, 1> combinations");
    test_combinations_bitset<uint16_t, 1>();

    xdiag::Log("Testing Bitset<uint32_t, 1> combinations");
    test_combinations_bitset<uint32_t, 1>();

    xdiag::Log("Testing Bitset<uint64_t, 1> combinations");
    test_combinations_bitset<uint64_t, 1>();

    xdiag::Log("Testing Bitset<uint64_t, 2> combinations");
    test_combinations_bitset<uint64_t, 2>();
  }

  SECTION("bitset large n") {
    xdiag::Log("Testing Bitset<uint64_t, 2> large n");
    test_combinations_bitset_large<uint64_t, 2>();

    xdiag::Log("Testing Bitset<uint64_t, 4> large n");
    test_combinations_bitset_large<uint64_t, 4>();

    xdiag::Log("Testing Bitset<uint64_t, 8> large n");
    test_combinations_bitset_large<uint64_t, 8>();
  }
}
