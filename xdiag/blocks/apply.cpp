// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "apply.hpp"

#include <xdiag/armadillo.hpp>
#include <xdiag/blocks/apply_block.hpp>
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
        apply_block(ops, bin, vec_in, bout, vec_out);
      },
      "Type mismatch of Block types");
}
XDIAG_CATCH

} // namespace xdiag

using namespace xdiag;
using namespace arma;

#define INSTANTIATE_XDIAG_APPLY(OP_TYPE, MAT_TYPE)                             \
  template void xdiag::apply(OP_TYPE const &, Block const &, MAT_TYPE const &, \
                             Block const &, MAT_TYPE &);

// BEGIN_INSTANTIATION_GROUP(op)
INSTANTIATE_XDIAG_APPLY(Op, vec)
INSTANTIATE_XDIAG_APPLY(Op, cx_vec)
INSTANTIATE_XDIAG_APPLY(Op, mat)
INSTANTIATE_XDIAG_APPLY(Op, cx_mat)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(opsum)
INSTANTIATE_XDIAG_APPLY(OpSum, vec)
INSTANTIATE_XDIAG_APPLY(OpSum, cx_vec)
INSTANTIATE_XDIAG_APPLY(OpSum, mat)
INSTANTIATE_XDIAG_APPLY(OpSum, cx_mat)
// END_INSTANTIATION_GROUP
