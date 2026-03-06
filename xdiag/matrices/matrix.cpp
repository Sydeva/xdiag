// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "matrix.hpp"

#include <xdiag/armadillo.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/matrices/spinhalf/dispatch_basis.hpp>
#include <xdiag/matrices/spinhalf/matrix_generic.hpp>
#include <xdiag/matrices/utils/fill_functions.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/operators/qns/block.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/variants.hpp>

// Dispatch overview
// -----------------
// The matrix function routes through two layers of type erasure before
// reaching the actual numerical kernel, mirroring apply.cpp exactly:
//
// Layer 1 — Block type (variant dispatch via visit_same_type):
//   matrix(op_t, Block)  →  matrix(op_t, Block, Block)
//     The output Block is first determined from the operator quantum numbers,
//     then visit_same_type unwraps both Block variants to their concrete type
//     and calls matrix(OpSum, ConcreteBlock, ConcreteBlock, coeff_t*).
//     Op/Monomial are promoted to OpSum here. The output arma matrix is
//     allocated and zeroed at this layer; a raw pointer is passed down so the
//     kernel can fill it without owning the storage.
//
// Layer 2 — Basis type (runtime dispatch via dispatch_basis):
//   matrix(OpSum, Spinhalf, Spinhalf, coeff_t*)
//     dispatch_basis resolves the shared_ptr<Basis> stored in each Spinhalf
//     to a concrete BasisOnTheFly<...> type via a type-id lookup table, then
//     calls the lambda with the concrete basis objects.
//
// Kernel — matrix_generic<coeff_t>(ops, basis_in, basis_out, fill_f):
//     The innermost lambda closes over the raw matrix pointer and calls
//     fill_matrix, which performs mat[idx_out + idx_in*m] += val, building
//     the dense matrix in column-major order.

namespace xdiag {

template <typename op_t>
arma::mat matrix(op_t const &op, Block const &blocki) try {
  auto blockr = block(OpSum(op), blocki);
  return matrix(op, blocki, blockr);
}
XDIAG_CATCH

template <typename op_t>
arma::cx_mat matrixC(op_t const &op, Block const &blocki) try {
  auto blockr = block(OpSum(op), blocki);
  return matrixC(op, blocki, blockr);
}
XDIAG_CATCH

template <typename op_t>
arma::mat matrix(op_t const &op, Block const &block_in,
                 Block const &block_out) try {
  int64_t m = size(block_out);
  int64_t n = dim(block_in);
  arma::mat mat(m, n, arma::fill::zeros); // zeroed here; kernel uses +=
  // Layer 1: unwrap the Block variant; op_t is promoted to OpSum inside.
  utils::visit_same_type(
      block_in, block_out,
      [&](auto const &bin, auto const &bout) {
        matrix(OpSum(op), bin, bout, mat.memptr());
      },
      "Type mismatch of Block types");
  return mat;
}
XDIAG_CATCH

template <typename op_t>
arma::cx_mat matrixC(op_t const &op, Block const &block_in,
                     Block const &block_out) try {
  int64_t m = size(block_out);
  int64_t n = dim(block_in);
  arma::cx_mat mat(m, n, arma::fill::zeros); // zeroed here; kernel uses +=
  // Layer 1: unwrap the Block variant; op_t is promoted to OpSum inside.
  utils::visit_same_type(
      block_in, block_out,
      [&](auto const &bin, auto const &bout) {
        matrix(OpSum(op), bin, bout, mat.memptr());
      },
      "Type mismatch of Block types");
  return mat;
}
XDIAG_CATCH

template <typename coeff_t>
void matrix(OpSum const &ops, Spinhalf const &block_in,
            Spinhalf const &block_out, coeff_t *mat) try {
  // Layer 2: unwrap the basis pointer; the lambda receives the concrete
  // BasisOnTheFly<...> type, allowing matrix_generic to be instantiated
  // at compile time for each basis specialization.
  matrices::spinhalf::dispatch_basis(
      block_in, block_out, [&](auto const &basis_in, auto const &basis_out) {
        int64_t m = basis_out.size();
        // Kernel: fill_matrix accumulates mat[idx_out + idx_in*m] += val,
        // building the dense matrix in column-major order.
        matrices::spinhalf::matrix_generic<coeff_t>(
            ops, basis_in, basis_out,
            [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
              matrices::fill_matrix(mat, m, idx_in, idx_out, val);
            });
      });
}
XDIAG_CATCH

template arma::mat matrix(Op const &, Block const &);
template arma::mat matrix(Monomial const &, Block const &);
template arma::mat matrix(OpSum const &, Block const &);

template arma::cx_mat matrixC(Op const &, Block const &);
template arma::cx_mat matrixC(Monomial const &, Block const &);
template arma::cx_mat matrixC(OpSum const &, Block const &);

template arma::mat matrix(Op const &op, Block const &, Block const &);
template arma::mat matrix(Monomial const &op, Block const &, Block const &);
template arma::mat matrix(OpSum const &op, Block const &, Block const &);

template arma::cx_mat matrixC(Op const &op, Block const &, Block const &);
template arma::cx_mat matrixC(Monomial const &op, Block const &, Block const &);
template arma::cx_mat matrixC(OpSum const &op, Block const &, Block const &);

template void matrix(OpSum const &, Spinhalf const &, Spinhalf const &,
                     double *);
template void matrix(OpSum const &, Spinhalf const &, Spinhalf const &,
                     complex *);

} // namespace xdiag
