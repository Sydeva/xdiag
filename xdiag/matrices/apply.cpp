// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "apply.hpp"

#include <xdiag/armadillo.hpp>
#include <xdiag/matrices/spinhalf/dispatch_basis.hpp>
#include <xdiag/matrices/spinhalf/matrix_generic.hpp>
#include <xdiag/matrices/utils/fill_functions.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/variants.hpp>

// Dispatch overview
// -----------------
// The apply function routes through two layers of type erasure before reaching
// the actual numerical kernel:
//
// Layer 1 — Block type (variant dispatch via visit_same_type):
//   apply(op_t, Block, mat_t, Block, mat_t)
//     Block is a std::variant<Spinhalf, ...>. visit_same_type unwraps both
//     block_in and block_out to their concrete type (enforcing they match),
//     then calls apply(OpSum, ConcreteBlock, ...). Op/Monomial are first
//     promoted to OpSum here.
//
// Layer 2 — Basis type (runtime dispatch via dispatch_basis):
//   apply(OpSum, Spinhalf, mat_t, Spinhalf, mat_t)
//     Spinhalf stores its basis as a shared_ptr<Basis>. dispatch_basis builds
//     a lookup table keyed on the basis type-id and invokes the matching entry,
//     giving the lambda the concrete BasisOnTheFly<...> objects.
//
// Kernel — matrix_generic<coeff_t>(ops, basis_in, basis_out, fill_f):
//     The innermost lambda closes over vec_in/vec_out and calls fill_apply,
//     which performs vec_out[idx_out] += val * vec_in[idx_in].

namespace xdiag {

template <typename op_t, typename mat_t>
void apply(op_t const &ops, Block const &block_in, mat_t const &vec_in,
           Block const &block_out, mat_t &vec_out) try {
  // Layer 1: unwrap the Block variant; op_t is promoted to OpSum inside.
  utils::visit_same_type(
      block_in, block_out,
      [&](auto const &bin, auto const &bout) {
        apply(OpSum(ops), bin, vec_in, bout, vec_out);
      },
      "Type mismatch of Block types");
}
XDIAG_CATCH

template <typename vec_t>
void apply(OpSum const &ops, Spinhalf const &block_in, vec_t const &vec_in,
           Spinhalf const &block_out, vec_t &vec_out) try {
  using coeff_t = typename vec_t::elem_type;
  vec_out.zeros();

  // Layer 2: unwrap the basis pointer; the lambda receives the concrete
  // BasisOnTheFly<...> type, allowing matrix_generic to be instantiated
  // at compile time for each basis specialization.
  matrices::spinhalf::dispatch_basis(
      block_in, block_out, [&](auto const &basis_in, auto const &basis_out) {
        // Kernel: fill_apply accumulates vec_out[idx_out] += val*vec_in[idx_in]
        matrices::spinhalf::matrix_generic<coeff_t>(
            ops, basis_in, basis_out,
            [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
              matrices::fill_apply(vec_in, vec_out, idx_in, idx_out, val);
            });
      });
}
XDIAG_CATCH

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
