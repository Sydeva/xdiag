// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "non_branching_op.hpp"

#include <limits>
#include <tuple>
#include <type_traits>
#include <vector>

#include <xdiag/bits/bitset.hpp>
#include <xdiag/bits/get_set_bit.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

namespace xdiag::basis {

template <typename coeff_t>
static bool is_non_branching_matrix(arma::Mat<coeff_t> const &mat,
                                    double precision) {
  for (arma::uword i = 0; i < mat.n_rows; ++i) {
    int64_t non_zero_in_row = 0;
    for (arma::uword j = 0; j < mat.n_cols; ++j) {
      if (std::abs(mat(i, j)) > precision) {
        ++non_zero_in_row;
      }
    }
    if (non_zero_in_row > 1) {
      return false;
    }
  }
  return true;
}

template <typename bit_t, typename coeff_t>
NonBranchingOp<bit_t, coeff_t>::NonBranchingOp(
    std::vector<int64_t> const &sites, arma::Mat<coeff_t> const &matrix,
    double precision) try
    : sites_(sites), mask_(0), diagonal_(true) {

  if (sites.size() > std::numeric_limits<cbit_t>::digits) {
    XDIAG_THROW(
        "Number of sites in NonBranchingOp too large for given bit type");
  }

  if (!is_non_branching_matrix(matrix, precision)) {
    XDIAG_THROW("Trying to create a NonBranchingOp from a matrix which is "
                "branching (i.e. more than one entry per row/column)");
  }

  // Set mask
  for (int64_t s : sites) {
    bits::set_bit(mask_, s);
  }
  mask_ = ~mask_;

  // Matrix dimension is 2**(no. sites of op)
  int64_t dim = (int64_t)1 << sites.size();
  if ((matrix.n_cols != dim) || (matrix.n_rows != dim)) {
    XDIAG_THROW(fmt::format(
        "Invalid matrix dimension for non-branching Op matrix. Expected "
        "dim={}, but received n_rows={} and n_cols={}.",
        dim, matrix.n_rows, matrix.n_cols));
  }

  // Set arrays where state is mapped to
  non_zero_term_ = std::vector<bool>(dim, false);
  state_applied_ = std::vector<cbit_t>(dim);
  coeff_ = std::vector<coeff_t>(dim, 0.);
  for (cbit_t in = 0; in < dim; ++in) {
    for (cbit_t out = 0; out < dim; ++out) {
      if (std::abs(matrix(out, in)) > precision) {
        non_zero_term_[in] = true;
        state_applied_[in] = out;
        coeff_[in] = matrix(out, in);
        break;
      }
    }
  }

  // Determine if diagonal
  for (cbit_t i = 0; i < dim; ++i) {
    if ((non_zero_term_[i]) && (state_applied_[i] != i)) {
      diagonal_ = false;
      break;
    }
  }
}
XDIAG_CATCH

template <typename bit_t, typename coeff_t>
bool NonBranchingOp<bit_t, coeff_t>::isdiagonal() const {
  return diagonal_;
}

template <typename bit_t, typename coeff_t>
bool NonBranchingOp<bit_t, coeff_t>::non_zero_term(cbit_t local_state) const {
  return non_zero_term_[local_state];
}
template <typename bit_t, typename coeff_t>
coeff_t NonBranchingOp<bit_t, coeff_t>::coeff(cbit_t local_state) const {
  return coeff_[local_state];
}

template <typename bit_t, typename coeff_t>
std::pair<typename NonBranchingOp<bit_t, coeff_t>::cbit_t, coeff_t>
NonBranchingOp<bit_t, coeff_t>::state_coeff(cbit_t local_state) const {
  return {state_applied_[local_state], coeff_[local_state]};
}

template <typename bit_t, typename coeff_t>
typename NonBranchingOp<bit_t, coeff_t>::cbit_t
NonBranchingOp<bit_t, coeff_t>::extract(bit_t state) const {
  cbit_t local_state = 0;
  for (int64_t i = 0; i < (int64_t)sites_.size(); ++i) {
    bits::set_bit(local_state, i, bits::get_bit(state, sites_[i]));
  }
  return local_state;
}

template <typename bit_t, typename coeff_t>
bit_t NonBranchingOp<bit_t, coeff_t>::deposit(cbit_t local_state,
                                              bit_t state) const {
  state &= mask_; // clear bits on site
  for (int64_t i = 0; i < (int64_t)sites_.size(); ++i) {
    bits::set_bit(state, sites_[i], bits::get_bit(local_state, i));
  }
  return state;
}

template <typename coeff_t>
static std::vector<arma::Mat<coeff_t>>
decompose_matrix_to_nonbranching(arma::Mat<coeff_t> const &mat,
                                 double precision) try {
  int64_t m = (int64_t)mat.n_rows;
  int64_t n = (int64_t)mat.n_cols;
  if (m != n) {
    XDIAG_THROW("Error: Op matrix is not square");
  }

  std::vector<std::tuple<int64_t, int64_t, coeff_t>> all_entries;

  // Get diagonal elements
  for (int64_t i = 0; i < n; ++i) {
    if (std::abs(mat(i, i)) > precision) {
      all_entries.push_back({i, i, mat(i, i)});
    }
  }

  // Get offidagonal elements
  for (int64_t n_diag = 1; n_diag < n; ++n_diag) {
    for (int64_t i = 0; i < n - n_diag; ++i) {
      if (std::abs(mat(i, i + n_diag)) > precision) {
        all_entries.push_back({i, i + n_diag, mat(i, i + n_diag)});
      }
      if (std::abs(mat(i + n_diag, i)) > precision) {
        all_entries.push_back({i + n_diag, i, mat(i + n_diag, i)});
      }
    }
  }

  // Reduce to minimal number of non-branching terms
  std::vector<arma::Mat<coeff_t>> mats_nb;

  while (all_entries.size() != 0) {
    std::vector<int64_t> forbidden_columns;
    std::vector<int64_t> forbidden_rows;
    std::vector<std::tuple<int64_t, int64_t, coeff_t>> current_entries;
    std::vector<int64_t> delete_entries;
    int64_t i = 0;
    for (auto [row, column, coeff] : all_entries) {

      if ((std::find(forbidden_rows.begin(), forbidden_rows.end(), row) ==
           forbidden_rows.end()) &&
          (std::find(forbidden_columns.begin(), forbidden_columns.end(),
                     column) == forbidden_columns.end())) {
        current_entries.push_back({row, column, coeff});
        forbidden_rows.push_back(row);
        forbidden_columns.push_back(column);
        delete_entries.push_back(i);
      }
      ++i;
    }

    for (int64_t i = delete_entries.size() - 1; i >= 0; --i)
      all_entries.erase(all_entries.begin() + delete_entries[i]);

    // Create non-branching matrix
    arma::Mat<coeff_t> mat_nb(m, n, arma::fill::zeros);
    for (auto [i, j, coeff] : current_entries) {
      mat_nb(i, j) = coeff;
    }
    mats_nb.push_back(mat_nb);
  }
  return mats_nb;
}
XDIAG_CATCH

static std::vector<Matrix>
decompose_matrix_to_nonbranching(Matrix const &mat, double precision) try {
  std::vector<Matrix> mats;
  if (mat.isreal()) {
    std::vector<arma::mat> mats2 =
        decompose_matrix_to_nonbranching(mat.as<arma::mat>(), precision);
    for (auto const &mat : mats2) {
      mats.push_back(Matrix(mat));
    }
  } else {
    std::vector<arma::cx_mat> mats2 =
        decompose_matrix_to_nonbranching(mat.as<arma::cx_mat>(), precision);
    for (auto const &mat : mats2) {
      mats.push_back(Matrix(mat));
    }
  }
  return mats;
}
XDIAG_CATCH

template <typename bit_t, typename coeff_t>
std::vector<NonBranchingOp<bit_t, coeff_t>>
non_branching_ops(Coeff const &cpl, Op const &op, double precision) try {
  if (cpl.isstring()) {
    XDIAG_THROW("Cannot convert Op to NonBranchingOps. Coeff found to be a "
                "string. To convert it to a NonBranchingOp the coupling must "
                "be either a real or complex number.");
  }

  std::vector<NonBranchingOp<bit_t, coeff_t>> ops;
  if (op.hasmatrix()) {
    auto mat = op.matrix() * cpl.scalar();
    auto sites = op.sites();
    auto mats_nb = decompose_matrix_to_nonbranching(mat, precision);
    for (auto m : mats_nb) {
      if constexpr (isreal<coeff_t>()) {
        if (m.isreal()) {
          ops.push_back(NonBranchingOp<bit_t, coeff_t>(sites, m.as<arma::mat>(),
                                                       precision));
        } else {
          XDIAG_THROW(
              "Cannot create a real NonBranchingOp from a complex matrix.")
        }
      } else {
        ops.push_back(NonBranchingOp<bit_t, coeff_t>(
            sites, m.as<arma::cx_mat>(), precision));
      }
    }
  } else {
    XDIAG_THROW(
        "Cannot convert Op to NonBranchingOps. Op has no matrix defined.");
  }
  return ops;
}
XDIAG_CATCH

#define INSTANTIATE_XDIAG_BASIS_NON_BRANCHING_OP(BIT_TYPE, NUMBER_TYPE)        \
  template class NonBranchingOp<BIT_TYPE, NUMBER_TYPE>;                        \
  template std::vector<NonBranchingOp<BIT_TYPE, NUMBER_TYPE>>                  \
  non_branching_ops<BIT_TYPE, NUMBER_TYPE>(Coeff const &, Op const &, double);

using namespace bits;

INSTANTIATE_XDIAG_BASIS_NON_BRANCHING_OP(uint32_t, double);
INSTANTIATE_XDIAG_BASIS_NON_BRANCHING_OP(uint64_t, double);
INSTANTIATE_XDIAG_BASIS_NON_BRANCHING_OP(BitsetDynamic, double);
INSTANTIATE_XDIAG_BASIS_NON_BRANCHING_OP(BitsetStatic2, double);
INSTANTIATE_XDIAG_BASIS_NON_BRANCHING_OP(BitsetStatic4, double);
INSTANTIATE_XDIAG_BASIS_NON_BRANCHING_OP(BitsetStatic8, double);

INSTANTIATE_XDIAG_BASIS_NON_BRANCHING_OP(uint32_t, complex);
INSTANTIATE_XDIAG_BASIS_NON_BRANCHING_OP(uint64_t, complex);
INSTANTIATE_XDIAG_BASIS_NON_BRANCHING_OP(BitsetDynamic, complex);
INSTANTIATE_XDIAG_BASIS_NON_BRANCHING_OP(BitsetStatic2, complex);
INSTANTIATE_XDIAG_BASIS_NON_BRANCHING_OP(BitsetStatic4, complex);
INSTANTIATE_XDIAG_BASIS_NON_BRANCHING_OP(BitsetStatic8, complex);

} // namespace xdiag::basis
