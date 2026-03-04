// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "apply_block.hpp"

#include <xdiag/armadillo.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/matrix/apply/apply_basis.hpp>
#include <xdiag/matrix/implementation/compile.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

namespace xdiag::matrix {

template <typename mat_t, typename block_t>
void apply_block(OpSum const &ops, block_t const &block_in, mat_t const &mat_in,
                 block_t const &block_out, mat_t &mat_out) try {
  using coeff_t = typename mat_t::elem_type;
  mat_out.zeros();
  OpSum ops_compiled = compile(ops, block_in);
  apply_basis(ops_compiled, block_in.basis(), mat_in, block_out.basis(),
              mat_out);
}
XDIAG_CATCH

} // namespace xdiag::matrix

using namespace xdiag;
using namespace arma;

#define INSTANTIATE_XDIAG_MATRIX_APPLY_BLOCK(BLOCK_TYPE, MAT_TYPE)             \
  template void xdiag::matrix::apply_block(OpSum const &, BLOCK_TYPE const &,  \
                                           MAT_TYPE const &,                   \
                                           BLOCK_TYPE const &, MAT_TYPE &);

// BEGIN_INSTANTIATION_GROUP(spinhalf)
INSTANTIATE_XDIAG_MATRIX_APPLY_BLOCK(Spinhalf, vec)
INSTANTIATE_XDIAG_MATRIX_APPLY_BLOCK(Spinhalf, cx_vec)
INSTANTIATE_XDIAG_MATRIX_APPLY_BLOCK(Spinhalf, mat)
INSTANTIATE_XDIAG_MATRIX_APPLY_BLOCK(Spinhalf, cx_mat)
// END_INSTANTIATION_GROUP

#undef INSTANTIATE_XDIAG_APPLY_BLOCK
