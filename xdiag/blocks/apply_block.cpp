// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "apply_block.hpp"

#include <xdiag/armadillo.hpp>
#include <xdiag/basis/apply.hpp>
#include <xdiag/blocks/spinhalf/compile.hpp>
#include <xdiag/blocks/spinhalf/spinhalf.hpp>
#include <xdiag/blocks/spinhalf/valid.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

namespace xdiag {

template <typename op_t, typename mat_t, typename block_t>
void apply_block(op_t const &ops, block_t const &block_in, mat_t const &mat_in,
                 block_t const &block_out, mat_t &mat_out) try {
  using coeff_t = typename mat_t::elem_type;
  mat_out.zeros();
  blocks::check_valid(block_in, ops);
  // OpSum opsc = blocks::compile(block_in, ops, "implementation");
  basis::apply(ops, block_in.basis(), mat_in, block_out.basis(), mat_out);
}
XDIAG_CATCH

} // namespace xdiag

using namespace xdiag;
using namespace arma;

#define INSTANTIATE_XDIAG_APPLY_BLOCK(OP_TYPE, BLOCK_TYPE, MAT_TYPE)           \
  template void xdiag::apply_block(OP_TYPE const &, BLOCK_TYPE const &,        \
                                   MAT_TYPE const &, BLOCK_TYPE const &,       \
                                   MAT_TYPE &);

// BEGIN_INSTANTIATION_GROUP(spinhalf)
INSTANTIATE_XDIAG_APPLY_BLOCK(Op, Spinhalf, vec)
INSTANTIATE_XDIAG_APPLY_BLOCK(Op, Spinhalf, cx_vec)
INSTANTIATE_XDIAG_APPLY_BLOCK(Op, Spinhalf, mat)
INSTANTIATE_XDIAG_APPLY_BLOCK(Op, Spinhalf, cx_mat)
INSTANTIATE_XDIAG_APPLY_BLOCK(OpSum, Spinhalf, vec)
INSTANTIATE_XDIAG_APPLY_BLOCK(OpSum, Spinhalf, cx_vec)
INSTANTIATE_XDIAG_APPLY_BLOCK(OpSum, Spinhalf, mat)
INSTANTIATE_XDIAG_APPLY_BLOCK(OpSum, Spinhalf, cx_mat)
// END_INSTANTIATION_GROUP

#undef INSTANTIATE_XDIAG_APPLY_BLOCK
