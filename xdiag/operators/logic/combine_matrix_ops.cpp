// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "combine_matrix_ops.hpp"

#include <map>
#include <vector>

#include <xdiag/armadillo.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/operators/logic/valid.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

namespace xdiag::operators {

// Embeds a single "Matrix" Op into the full dim_total × dim_total space.
//
// Indexing convention: full state bit i encodes the spin at all_sites[i].
// Op local state bit j encodes the spin at op.sites()[j].
//
// For a given full-space row r and column c, the entry is:
//   M[local_r, local_c]   if r and c agree on every bit NOT in the op's sites
//   0                     otherwise
// where local_r (local_c) is the local index extracted from r (c) by reading
// only those bits that correspond to the op's sites.
template <typename coeff_t>
static arma::Mat<coeff_t>
embed_op(Op const &op, std::vector<int64_t> const &all_sites,
         std::map<int64_t, int64_t> const &site_to_pos) {

  int64_t N = (int64_t)all_sites.size();
  int64_t dim_total = (int64_t)1 << N;

  auto const &op_sites = op.sites();
  int64_t k = (int64_t)op_sites.size();

  // Position of each of the op's sites within all_sites
  std::vector<int64_t> op_pos(k);
  for (int64_t j = 0; j < k; ++j) {
    op_pos[j] = site_to_pos.at(op_sites[j]);
  }

  // Bitmask of bits that belong to the op; complement = bits the op leaves
  // unchanged (must match between r and c for a non-zero entry)
  int64_t op_mask = 0;
  for (int64_t pos : op_pos) {
    op_mask |= (int64_t)1 << pos;
  }
  int64_t comp_mask = (~op_mask) & (dim_total - 1);

  auto const M = op.matrix().as<arma::Mat<coeff_t>>();
  arma::Mat<coeff_t> full(dim_total, dim_total, arma::fill::zeros);

  for (int64_t r = 0; r < dim_total; ++r) {
    for (int64_t c = 0; c < dim_total; ++c) {
      // Row and column must agree on all sites outside the op
      if ((r & comp_mask) != (c & comp_mask)) {
        continue;
      }
      // Extract local indices for the op's sites
      int64_t local_r = 0, local_c = 0;
      for (int64_t j = 0; j < k; ++j) {
        local_r |= ((r >> op_pos[j]) & 1) << j;
        local_c |= ((c >> op_pos[j]) & 1) << j;
      }
      full(r, c) = M(local_r, local_c);
    }
  }
  return full;
}

// Combines multiple "Matrix" Ops into a single "Matrix" Op.
//
// The combined op acts on the union of all sites (in first-appearance order).
// Its matrix is the product of the individual ops' matrices, each first
// embedded into the full tensor product space of those combined sites.
//
// Overlapping sites are handled correctly: the embedding of each op is an
// identity on all sites it does not touch, so the full matrix is well-defined
// even when two ops share a site.
Op combine_matrix_ops(std::vector<Op> const &ops) try {
  if (ops.empty()) {
    XDIAG_THROW("combine_matrix_ops: ops vector must not be empty");
  }

  for (auto const &op : ops) {
    if (op.type() != "Matrix") {
      XDIAG_THROW(fmt::format(
          "combine_matrix_ops: all ops must be of type \"Matrix\", "
          "got \"{}\"",
          op.type()));
    }
    must_have_sites(op);
    must_have_matrix(op);
  }

  // Collect unique sites in first-appearance order
  std::vector<int64_t> all_sites;
  std::map<int64_t, int64_t> site_to_pos;
  for (auto const &op : ops) {
    for (int64_t s : op.sites()) {
      if (site_to_pos.find(s) == site_to_pos.end()) {
        site_to_pos[s] = (int64_t)all_sites.size();
        all_sites.push_back(s);
      }
    }
  }

  // Result is complex if any matrix is complex
  bool all_real = true;
  for (auto const &op : ops) {
    if (!op.matrix().isreal()) {
      all_real = false;
      break;
    }
  }

  // Build the product of embedded matrices
  if (all_real) {
    arma::mat combined = embed_op<double>(ops[0], all_sites, site_to_pos);
    for (std::size_t i = 1; i < ops.size(); ++i) {
      combined = combined * embed_op<double>(ops[i], all_sites, site_to_pos);
    }
    return Op("Matrix", all_sites, combined);
  } else {
    arma::cx_mat combined =
        embed_op<complex>(ops[0], all_sites, site_to_pos);
    for (std::size_t i = 1; i < ops.size(); ++i) {
      combined =
          combined * embed_op<complex>(ops[i], all_sites, site_to_pos);
    }
    return Op("Matrix", all_sites, combined);
  }
}
XDIAG_CATCH

} // namespace xdiag::operators
