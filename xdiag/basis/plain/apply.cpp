// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "apply.hpp"

#include <xdiag/basis/plain/basis_onthefly.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

namespace xdiag::basis::plain {

template <typename basis_t>
void apply(OpSum const &ops, basis_t const &in, basis_t const &out) try {
  Log("{}\n{}\n", in.name(), out.name());
}
XDIAG_CATCH

using namespace combinatorics;
template void apply(OpSum const &, BasisOnTheFly<Subsets<uint32_t>> const &,
                    BasisOnTheFly<Subsets<uint32_t>> const &);
template void apply(OpSum const &, BasisOnTheFly<Subsets<uint64_t>> const &,
                    BasisOnTheFly<Subsets<uint64_t>> const &);
template void apply(OpSum const &,
                    BasisOnTheFly<Combinations<uint32_t>> const &,
                    BasisOnTheFly<Combinations<uint32_t>> const &);
template void apply(OpSum const &,
                    BasisOnTheFly<Combinations<uint64_t>> const &,
                    BasisOnTheFly<Combinations<uint64_t>> const &);
} // namespace xdiag::basis::plain
