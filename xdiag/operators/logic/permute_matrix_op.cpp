// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "permute_matrix_op.hpp"

#include <algorithm>
#include <numeric>
#include <vector>

#include <xdiag/math/ipow.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

namespace xdiag::operators {

// Maps a global Hilbert-space index to its image under the inverse of perm,
// where d is the local dimension per site.
static int64_t permuted_index(int64_t idx, int64_t d,
                              std::vector<int64_t> const &perm) {
  std::vector<int64_t> perm_inv(perm.size());
  for (int64_t i = 0; i < (int64_t)perm.size(); ++i) {
    perm_inv[perm[i]] = i;
  }
  int64_t exp = 0;
  int64_t idx_permuted = 0;
  while (idx) {
    int64_t local = idx % d;
    idx_permuted += local * math::ipow(d, perm_inv[exp++]);
    idx /= d;
  }
  return idx_permuted;
}

// Permutes rows and columns of a square matrix according to a site permutation.
template <typename T>
static arma::Mat<T> permute_matrix(arma::Mat<T> const &mat,
                                   std::vector<int64_t> const &perm) try {
  int64_t m = mat.n_rows;
  int64_t n = mat.n_cols;
  if (m != n) {
    XDIAG_THROW(fmt::format("Matrix is not square, size: ({}, {})", m, n));
  }

  // Determine local dimension d such that d^nsites == m
  int64_t d = 0;
  int64_t nsites = (int64_t)perm.size();
  while (math::ipow(d, nsites) < m) {
    ++d;
  }
  if (math::ipow(d, nsites) != m) {
    XDIAG_THROW(fmt::format(
        "Matrix dimensions are not of the form d^N, where N is the number "
        "of sites (here {})",
        nsites));
  }

  arma::Mat<T> mat_permuted(m, n);
  for (int64_t i = 0; i < m; ++i) {
    int64_t ip = permuted_index(i, d, perm);
    for (int64_t j = 0; j < m; ++j) {
      int64_t jp = permuted_index(j, d, perm);
      mat_permuted(ip, jp) = mat(i, j);
    }
  }
  return mat_permuted;
}
XDIAG_CATCH

Op permute_matrix_op(Op const &op) try {
  auto const &sites = op.sites();
  int64_t n = (int64_t)sites.size();

  // Compute the permutation that sorts sites in ascending order
  std::vector<int64_t> perm(n);
  std::iota(perm.begin(), perm.end(), 0);
  std::sort(perm.begin(), perm.end(),
            [&](int64_t a, int64_t b) { return sites[a] < sites[b]; });

  // Apply permutation to site list
  std::vector<int64_t> sites_sorted(n);
  for (int64_t i = 0; i < n; ++i) {
    sites_sorted[i] = sites[perm[i]];
  }

  Matrix const &mat = op.matrix();
  if (mat.isreal()) {
    return Op("Matrix", sites_sorted,
              permute_matrix(mat.as<arma::mat>(), perm));
  } else {
    return Op("Matrix", sites_sorted,
              permute_matrix(mat.as<arma::cx_mat>(), perm));
  }
}
XDIAG_CATCH

} // namespace xdiag::operators
