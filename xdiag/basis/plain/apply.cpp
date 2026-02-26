// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "apply.hpp"

#include <cstdint>

#include <xdiag/armadillo.hpp>
#include <xdiag/basis/fill_functions.hpp>
#include <xdiag/basis/plain/basis_onthefly.hpp>
#include <xdiag/basis/plain/implementation/apply_generic.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::basis::plain {

template <typename mat_t, typename basis_t>
void apply(OpSum const &ops, basis_t const &basis_in, mat_t const &mat_in,
           basis_t const &basis_out, mat_t &mat_out) try {
  using coeff_t = typename mat_t::elem_type;
  apply_generic<coeff_t>(ops, basis_in, basis_out,
                         [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
                           fill_apply(mat_in, mat_out, idx_in, idx_out, val);
                         });
}
XDIAG_CATCH

using namespace combinatorics;
void apply(OpSum const &ops, BasisOnTheFly<Subsets<uint32_t>> const &basis_in,
           arma::vec const &mat_in,
           BasisOnTheFly<Subsets<uint32_t>> const &basis_out,
           arma::vec &mat_out);
void apply(OpSum const &ops, BasisOnTheFly<Subsets<uint32_t>> const &basis_in,
           arma::cx_vec const &mat_in,
           BasisOnTheFly<Subsets<uint32_t>> const &basis_out,
           arma::cx_vec &mat_out);
void apply(OpSum const &ops, BasisOnTheFly<Subsets<uint64_t>> const &basis_in,
           arma::vec const &mat_in,
           BasisOnTheFly<Subsets<uint64_t>> const &basis_out,
           arma::vec &mat_out);
void apply(OpSum const &ops, BasisOnTheFly<Subsets<uint64_t>> const &basis_in,
           arma::cx_vec const &mat_in,
           BasisOnTheFly<Subsets<uint64_t>> const &basis_out,
           arma::cx_vec &mat_out);

void apply(OpSum const &ops,
           BasisOnTheFly<Combinations<uint32_t>> const &basis_in,
           arma::mat const &mat_in,
           BasisOnTheFly<Combinations<uint32_t>> const &basis_out,
           arma::mat &mat_out);
void apply(OpSum const &ops,
           BasisOnTheFly<Combinations<uint32_t>> const &basis_in,
           arma::cx_mat const &mat_in,
           BasisOnTheFly<Combinations<uint32_t>> const &basis_out,
           arma::cx_mat &mat_out);
void apply(OpSum const &ops,
           BasisOnTheFly<Combinations<uint64_t>> const &basis_in,
           arma::mat const &mat_in,
           BasisOnTheFly<Combinations<uint64_t>> const &basis_out,
           arma::mat &mat_out);
void apply(OpSum const &ops,
           BasisOnTheFly<Combinations<uint64_t>> const &basis_in,
           arma::cx_mat const &mat_in,
           BasisOnTheFly<Combinations<uint64_t>> const &basis_out,
           arma::cx_mat &mat_out);

void apply(OpSum const &ops, BasisOnTheFly<LinTable<uint32_t>> const &basis_in,
           arma::mat const &mat_in,
           BasisOnTheFly<LinTable<uint32_t>> const &basis_out,
           arma::mat &mat_out);
void apply(OpSum const &ops, BasisOnTheFly<LinTable<uint32_t>> const &basis_in,
           arma::cx_mat const &mat_in,
           BasisOnTheFly<LinTable<uint32_t>> const &basis_out,
           arma::cx_mat &mat_out);
void apply(OpSum const &ops, BasisOnTheFly<LinTable<uint64_t>> const &basis_in,
           arma::mat const &mat_in,
           BasisOnTheFly<LinTable<uint64_t>> const &basis_out,
           arma::mat &mat_out);
void apply(OpSum const &ops, BasisOnTheFly<LinTable<uint64_t>> const &basis_in,
           arma::cx_mat const &mat_in,
           BasisOnTheFly<LinTable<uint64_t>> const &basis_out,
           arma::cx_mat &mat_out);

void apply(OpSum const &ops, BasisOnTheFly<Subsets<uint32_t>> const &basis_in,
           arma::mat const &mat_in,
           BasisOnTheFly<Subsets<uint32_t>> const &basis_out,
           arma::mat &mat_out);
void apply(OpSum const &ops, BasisOnTheFly<Subsets<uint32_t>> const &basis_in,
           arma::cx_mat const &mat_in,
           BasisOnTheFly<Subsets<uint32_t>> const &basis_out,
           arma::cx_mat &mat_out);
void apply(OpSum const &ops, BasisOnTheFly<Subsets<uint64_t>> const &basis_in,
           arma::mat const &mat_in,
           BasisOnTheFly<Subsets<uint64_t>> const &basis_out,
           arma::mat &mat_out);
void apply(OpSum const &ops, BasisOnTheFly<Subsets<uint64_t>> const &basis_in,
           arma::cx_mat const &mat_in,
           BasisOnTheFly<Subsets<uint64_t>> const &basis_out,
           arma::cx_mat &mat_out);

void apply(OpSum const &ops,
           BasisOnTheFly<Combinations<uint32_t>> const &basis_in,
           arma::mat const &mat_in,
           BasisOnTheFly<Combinations<uint32_t>> const &basis_out,
           arma::mat &mat_out);
void apply(OpSum const &ops,
           BasisOnTheFly<Combinations<uint32_t>> const &basis_in,
           arma::cx_mat const &mat_in,
           BasisOnTheFly<Combinations<uint32_t>> const &basis_out,
           arma::cx_mat &mat_out);
void apply(OpSum const &ops,
           BasisOnTheFly<Combinations<uint64_t>> const &basis_in,
           arma::mat const &mat_in,
           BasisOnTheFly<Combinations<uint64_t>> const &basis_out,
           arma::mat &mat_out);
void apply(OpSum const &ops,
           BasisOnTheFly<Combinations<uint64_t>> const &basis_in,
           arma::cx_mat const &mat_in,
           BasisOnTheFly<Combinations<uint64_t>> const &basis_out,
           arma::cx_mat &mat_out);

void apply(OpSum const &ops, BasisOnTheFly<LinTable<uint32_t>> const &basis_in,
           arma::mat const &mat_in,
           BasisOnTheFly<LinTable<uint32_t>> const &basis_out,
           arma::mat &mat_out);
void apply(OpSum const &ops, BasisOnTheFly<LinTable<uint32_t>> const &basis_in,
           arma::cx_mat const &mat_in,
           BasisOnTheFly<LinTable<uint32_t>> const &basis_out,
           arma::cx_mat &mat_out);
void apply(OpSum const &ops, BasisOnTheFly<LinTable<uint64_t>> const &basis_in,
           arma::mat const &mat_in,
           BasisOnTheFly<LinTable<uint64_t>> const &basis_out,
           arma::mat &mat_out);
void apply(OpSum const &ops, BasisOnTheFly<LinTable<uint64_t>> const &basis_in,
           arma::cx_mat const &mat_in,
           BasisOnTheFly<LinTable<uint64_t>> const &basis_out,
           arma::cx_mat &mat_out);

} // namespace xdiag::basis::plain
