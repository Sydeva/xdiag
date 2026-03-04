// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "compile.hpp"

#include <xdiag/operators/logic/algebra.hpp>
#include <xdiag/operators/logic/normal_order.hpp>
#include <xdiag/operators/logic/valid.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::blocks {

OpSum compile(OpSum const &ops, Spinhalf const &block) try {
  operators::check_valid(ops);
  auto impl_alg = operators::spinhalf_implementation_algebra();
  return normal_order(ops.plain(), impl_alg, block.nsites());
}
XDIAG_CATCH

} // namespace xdiag::blocks
