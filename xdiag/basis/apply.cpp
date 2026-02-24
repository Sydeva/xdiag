// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "apply.hpp"

#include <xdiag/basis/dispatcher.hpp>
#include <xdiag/basis/plain/apply.hpp>
#include <xdiag/basis/plain/basis_onthefly.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::basis {

void apply(OpSum const &ops, std::shared_ptr<Basis> const &in,
           std::shared_ptr<Basis> const &out) try {
  using namespace combinatorics;

  Dispatcher d;
#define ADD_DISPATCH(BASIS, FUNCTION)                                          \
  d.add<BASIS>(                                                                \
      [&](BASIS const &in, BASIS const &out) { FUNCTION(ops, in, out); });

  ADD_DISPATCH(BasisOnTheFly<Subsets<uint32_t>>, plain::apply);
  ADD_DISPATCH(BasisOnTheFly<Subsets<uint64_t>>, plain::apply);
  ADD_DISPATCH(BasisOnTheFly<Combinations<uint32_t>>, plain::apply);
  ADD_DISPATCH(BasisOnTheFly<Combinations<uint64_t>>, plain::apply);
#undef ADD_DISPATCH

  d.dispatch(in, out);
}
XDIAG_CATCH

} // namespace xdiag::basis
