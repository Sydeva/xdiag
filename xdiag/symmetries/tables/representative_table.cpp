// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "representative_table.hpp"

#include <algorithm>
#include <cstdint>
#include <type_traits>

#include <xdiag/bits/log2.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>

namespace xdiag::symmetries {

template <typename state_indexing_t, typename coeff_t>
static void representative_table_initialize(
    state_indexing_t const &state_indexing, GroupAction const &group_action,
    std::vector<coeff_t> const &characters,
    bits::BitVector<typename state_indexing_t::bit_t> &representative,
    bits::BitVector<typename state_indexing_t::bit_t> &representative_index,
    bits::BitVector<typename state_indexing_t::bit_t> &representative_symmetry,
    bits::BitVector<typename state_indexing_t::bit_t>
        &representative_norm_index,
    std::vector<double> norm) try {
  using bit_t = typename state_indexing_t::bit_t;
  using bits::BitVector;

  // Create vector holding the indices for each state yielding the
  // representative
  auto states = state_indexing.states();
  int64_t nbits = states.n();
  int64_t size = state_indexing.size();
  try {
    representative_index = BitVector<bit_t>(nbits, size);
  } catch (...) {
    XDIAG_THROW("Unable to allocate representative index array");
  }

  // First pass, simply count number of representatives and number of different
  // norms
  int64_t nrepresentatives = 0;
  for (auto state : states) {
    if (is_representative(state, group_action)) {
      double nrm = norm(state, group_action, characters);

      if (std::fabs(nrm) > 1e-6) { // representative found
        ++nrepresentatives;

        // Check if norm already registered ...
        auto it = std::find_if(norm.begin(), norm.end(), [&](double n) {
          return std::fabs(n - nrm) < 1e-6;
        });

        // ... if not, register it
        if (it == norm.end()) {
          norm.push_back(nrm);
        }
      }
    }
  }
  std::sort(norm.begin(), norm.end());

  // Create vector holding the representatives and norm_index
  try {
    representative = BitVector<bit_t>(nbits, nrepresentatives);
  } catch (...) {
    XDIAG_THROW("Unable to allocate representative array");
  }

  int64_t nbits_for_norm = bits::ceillog2(norm.size());
  try {
    representative_norm_index =
        BitVector<bit_t>(nbits_for_norm, nrepresentatives);
  } catch (...) {
    XDIAG_THROW("Unable to allocate representative norm index array");
  }

  // Second pass, fill all vectors
  int64_t idx = 0;
  nrepresentatives = 0;
  for (auto state : states) {
    if (is_representative(state, group_action)) {
      double nrm = norm(state, group_action, characters);

      if (std::fabs(nrm) > 1e-6) { // representative found
        representative_index[idx] = nrepresentatives;
        representative[nrepresentatives] = state;

        // Get index of norm
        auto it = std::find_if(norm.begin(), norm.end(), [&](double n) {
          return std::fabs(n - nrm) < 1e-6;
        });
        representative_norm_index[nrepresentatives] =
            std::distance(it, norm.begin());
        ++nrepresentatives;
      }
    }
    ++idx;
  }

  // Compute the symmetries that yield the representative and fill
  // non-representative representative_index
  int64_t nbits_for_symmetry = bits::ceillog2(group_action.n_symmetries());
  try {
    representative_symmetry =
        BitVector<bit_t>(nbits_for_symmetry, nrepresentatives);
  } catch (...) {
    XDIAG_THROW("Unable to allocate representative symmetry array");
  }

  auto const &group = group_action.permutation_group();
  for (int64_t rep_idx = 0; rep_idx < representative.size(); ++rep_idx) {
    bit_t rep = representative[rep_idx];
    for (int64_t sym = 0; sym < group_action.n_symmetries(); ++sym) {
      bit_t state = group_action.apply(sym, rep);
      int64_t idx = state_indexing.index(state);
      representative_index[idx] = rep_idx;
      representative_symmetry[idx] = group.inv(sym);
    }
  }
}
XDIAG_CATCH

template <typename state_indexing_t>
RepresentativeTable<state_indexing_t>::RepresentativeTable(
    state_indexing_t const &state_indexing, GroupAction const &group_action,
    std::vector<double> const &characters) try {
  representative_table_initialize(state_indexing, group_action, characters,
                                  representative_, representative_index_,
                                  representative_symmetry_,
                                  representative_norm_index_, norm_);
}
XDIAG_CATCH

template <typename state_indexing_t>
RepresentativeTable<state_indexing_t>::RepresentativeTable(
    state_indexing_t const &state_indexing, GroupAction const &group_action,
    std::vector<complex> const &characters) try {
  representative_table_initialize(state_indexing, group_action, characters,
                                  representative_, representative_index_,
                                  representative_symmetry_,
                                  representative_norm_index_, norm_);
}
XDIAG_CATCH

int64_t RepresentativeTable<state_indexing_t>::size() const {
  return representative_.size();
}

template <typename state_indexing_t>
bool RepresentativeTable<state_indexing_t>::operator==(
    RepresentativeTable<state_indexing_t> const &rhs) const {
  return (representative_ == rhs.representative_) &&
         (representative_index_ == rhs.representative_index_) &&
         (representative_symmetry_ == rhs.representative_symmetry_) &&
         (representative_norm_index_ == rhs.representative_norm_index_) &&
         (norm_ == rhs.norm_);
}

template <typename state_indexing_t>
bool RepresentativeTable<state_indexing_t>::operator!=(
    RepresentativeTable<state_indexing_t> const &rhs) const {
  return ~operator==(rhs);
}

using namespace combinatorics;
template class RepresentativeTable<Subsets<uint32_t>>;
template class RepresentativeTable<Subsets<uint64_t>>;
template class RepresentativeTable<Combinations<uint32_t>>;
template class RepresentativeTable<Combinations<uint64_t>>;

} // namespace xdiag::symmetries
