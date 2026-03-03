// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include "testcases_spinhalf.hpp"

#include <xdiag/algorithms/sparse_diag.hpp>
#include <xdiag/blocks/spinhalf/spinhalf.hpp>
#include <xdiag/operators/logic/algebra.hpp>
#include <xdiag/operators/logic/isapprox.hpp>
#include <xdiag/operators/logic/normal_order.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/xdiag_show.hpp>

TEST_CASE("spinhalf", "[blocks]") try {
  using namespace xdiag;
  auto b = Spinhalf(3);
  for (auto s : b) {
    Log("{}", to_string(s));
  }

  auto b2 = Spinhalf(4, 2);
  for (auto s : b2) {
    Log("{}", to_string(s));
  }

  // OpSum ops = complex(0.0, 1.0) * (Op("SdotS", {0, 1}) * Op("SdotS", {1, 2})
  // -
  //                                  Op("SdotS", {1, 2}) * Op("SdotS", {0,
  //                                  1}));
  // OpSum o1 = normal_order(ops, operators::spin_algebra(), 3);
  // OpSum o2 = normal_order(Op("ScalarChirality", {0, 1, 2}),
  //                         operators::spin_algebra(), 3);
  // XDIAG_SHOW(o1);
  // XDIAG_SHOW(o2);
  // XDIAG_SHOW(o2 +
  //            3.21 * Op("SdotS", {0, 1}) *
  //                Op("Matrix", 0, arma::mat{{0, 1}, {1, 0}}) +
  //            "U" * Op("HubbardU"));

  // XDIAG_SHOW(isapprox(o1, o2));

  auto ops = testcases::HBchain(4, 1.0);
  Log("hello");
  Log.set_verbosity(2);

  double e0 = eigval0(ops, b2);
  Log("e0: {}", e0);

} catch (xdiag::Error e) {
  error_trace(e);
}
