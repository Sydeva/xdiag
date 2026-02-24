// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include <cstdint>

#include <xdiag/basis/apply.hpp>
#include <xdiag/basis/plain/basis_onthefly.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

TEST_CASE("basis_onthefly", "[basis]") try {
  using namespace xdiag;
  using namespace xdiag::basis;
  using namespace xdiag::combinatorics;

  Basis *b1 = new BasisOnTheFly(Subsets<uint32_t>(3));
  Basis *b2 = new BasisOnTheFly(Subsets<uint64_t>(3));
  Basis *b3 = new BasisOnTheFly(Combinations<uint32_t>(3, 2));
  Basis *b4 = new BasisOnTheFly(Combinations<uint64_t>(3, 2));

  Log("type: {} name: {}", b1->type(), b1->name());
  Log("type: {} name: {}", b2->type(), b2->name());
  Log("type: {} name: {}", b3->type(), b3->name());
  Log("type: {} name: {}", b4->type(), b4->name());

  auto ops = OpSum();
  apply(ops, b1, b1);
  apply(ops, b2, b2);
  apply(ops, b3, b3);
  apply(ops, b4, b4);
  apply(ops, b1, b4);
} catch (xdiag::Error e) {
  error_trace(e);
}
