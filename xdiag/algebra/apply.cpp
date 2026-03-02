// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "apply.hpp"

#include <xdiag/basis/apply.hpp>

namespace xdiag {

template <typename op_t, typename mat_t, typename block_t>
void apply(op_t const &ops, block_t const &block_in, mat_t const &mat_in,
           block_t const &block_out, mat_t &mat_out) try {
  using coeff_t = typename mat_t::elem_type;

  check_valid(ops, block_in.nsites());
  mat_out.zeros();
  // OpSum opsc = operators::compile<block_t>(ops);
  OpSum opsc = ops;
  basis::apply(opsc, block_in.basis(), mat_in, block_out.basis(), mat_out);
}
XDIAG_CATCH

} // namespace xdiag

using namespace xdiag;
using namespace arma;

#define INSTANTIATE_XDIAG_APPLY(OP_TYPE, BLOCK_TYPE, MAT_TYPE)                 \
  template void xdiag::apply(OP_TYPE const &, BLOCK_TYPE const &,              \
                             MAT_TYPE const &, BLOCK_TYPE const &,             \
                             MAT_TYPE &);

// BEGIN_INSTANTIATION_GROUP(spinhalf)
INSTANTIATE_XDIAG_APPLY(Op, Spinhalf, vec)
INSTANTIATE_XDIAG_APPLY(Op, Spinhalf, cx_vec)
INSTANTIATE_XDIAG_APPLY(Op, Spinhalf, mat)
INSTANTIATE_XDIAG_APPLY(Op, Spinhalf, cx_mat)
INSTANTIATE_XDIAG_APPLY(OpSum, Spinhalf, vec)
INSTANTIATE_XDIAG_APPLY(OpSum, Spinhalf, cx_vec)
INSTANTIATE_XDIAG_APPLY(OpSum, Spinhalf, mat)
INSTANTIATE_XDIAG_APPLY(OpSum, Spinhalf, cx_mat)
// END_INSTANTIATION_GROUP

#undef INSTANTIATE_XDIAG_APPLY
