// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include <xdiag/blocks/spinhalf/spinhalf.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

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

} catch (xdiag::Error e) {
  error_trace(e);
}
