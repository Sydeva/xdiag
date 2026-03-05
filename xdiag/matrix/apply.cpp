// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "apply.hpp"

#include <xdiag/armadillo.hpp>
#include <xdiag/basis/basis_onthefly.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/matrix/spinhalf/apply.hpp>
#include <xdiag/matrix/utils/dispatcher.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/variants.hpp>

namespace xdiag {

template <typename op_t, typename mat_t>
void apply(op_t const &ops, Block const &block_in, mat_t const &vec_in,
           Block const &block_out, mat_t &vec_out) try {
  utils::visit_same_type(
      block_in, block_out,
      [&](auto const &bin, auto const &bout) {
        apply(OpSum(ops), bin, vec_in, bout, vec_out);
      },
      "Type mismatch of Block types");
}
XDIAG_CATCH

template <typename vec_t>
XDIAG_API void apply(OpSum const &ops, Spinhalf const &block_in,
                     vec_t const &vec_in, Spinhalf const &block_out,
                     vec_t &vec_out) {
  using namespace basis;
  using namespace combinatorics;
  using namespace bits;

  matrix::Dispatcher d;
#define ADD_DISPATCH(BASIS)                                                    \
  d.add<BASIS>([&](BASIS const &basis_in, BASIS const &basis_out) {            \
    matrix::spinhalf::apply(ops, basis_in, vec_in, basis_out, vec_out);        \
  });
  ADD_DISPATCH(BasisOnTheFly<Subsets<uint32_t>>);
  ADD_DISPATCH(BasisOnTheFly<Subsets<uint64_t>>);
  ADD_DISPATCH(BasisOnTheFly<Combinations<uint32_t>>);
  ADD_DISPATCH(BasisOnTheFly<Combinations<uint64_t>>);
  ADD_DISPATCH(BasisOnTheFly<LinTable<uint32_t>>);
  ADD_DISPATCH(BasisOnTheFly<LinTable<uint64_t>>);
  ADD_DISPATCH(BasisOnTheFly<Combinations<BitsetDynamic>>);
  ADD_DISPATCH(BasisOnTheFly<Combinations<BitsetStatic2>>);
  ADD_DISPATCH(BasisOnTheFly<Combinations<BitsetStatic4>>);
  ADD_DISPATCH(BasisOnTheFly<Combinations<BitsetStatic8>>);
#undef ADD_DISPATCH
  d.dispatch(block_in.basis(), block_out.basis());
}

} // namespace xdiag

using namespace xdiag;
using namespace arma;

#define INSTANTIATE_XDIAG_APPLY(OP_TYPE, MAT_TYPE)                             \
  template void xdiag::apply(OP_TYPE const &, Block const &, MAT_TYPE const &, \
                             Block const &, MAT_TYPE &);

INSTANTIATE_XDIAG_APPLY(Op, vec)
INSTANTIATE_XDIAG_APPLY(Op, cx_vec)
INSTANTIATE_XDIAG_APPLY(Op, mat)
INSTANTIATE_XDIAG_APPLY(Op, cx_mat)

INSTANTIATE_XDIAG_APPLY(Monomial, vec)
INSTANTIATE_XDIAG_APPLY(Monomial, cx_vec)
INSTANTIATE_XDIAG_APPLY(Monomial, mat)
INSTANTIATE_XDIAG_APPLY(Monomial, cx_mat)

INSTANTIATE_XDIAG_APPLY(OpSum, vec)
INSTANTIATE_XDIAG_APPLY(OpSum, cx_vec)
INSTANTIATE_XDIAG_APPLY(OpSum, mat)
INSTANTIATE_XDIAG_APPLY(OpSum, cx_mat)
