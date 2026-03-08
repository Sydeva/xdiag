// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../catch.hpp"

#include <iostream>
#include <random>

#include <xdiag/symmetries/permutation.hpp>
#include <xdiag/utils/logger.hpp>

TEST_CASE("permutation", "[symmetries]") {
  using namespace xdiag;
  Log("Test Permutation");

  std::random_device rd;
  std::mt19937 g(rd());

  // Test if identity is correct
  for (int64_t nsites = 1; nsites < 8; ++nsites) {
    std::vector<int64_t> pv(nsites);
    for (int64_t i = 0; i < nsites; ++i) {
      pv[i] = i;
    }
    auto p1 = Permutation(pv);
    auto p2 = Permutation(nsites);
    REQUIRE(p1 == p2);
  }

  for (int64_t nsites = 1; nsites < 8; ++nsites) {

    // Test identity multiplies
    for (int64_t i = 0; i < 5; ++i) {
      auto id = Permutation(nsites);
      auto a = id.array();
      std::shuffle(a.begin(), a.end(), g);
      auto p1 = Permutation(a);
      auto pi = Permutation(nsites);
      REQUIRE(p1 * pi == p1);
      REQUIRE(pi * p1 == p1);
    }

    // Test inverse
    for (int64_t i = 0; i < 20; ++i) {
      auto id = Permutation(nsites);
      auto a = id.array();
      std::shuffle(a.begin(), a.end(), g);
      auto p = Permutation(a);
      auto pinv = inv(p);
      REQUIRE(p * pinv == id);
    }
  }

  Log("done");
}
