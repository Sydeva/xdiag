// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "representative_table.hpp"

#include <algorithm>
#include <cstdint>
#include <type_traits>

#include <xdiag/bits/log2.hpp>
#include <xdiag/bits/nbits.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/symmetries/action/isrepresentative.hpp>
#include <xdiag/symmetries/action/norm.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::symmetries {

template <typename enumeration_t, typename coeff_t>
static void representative_table_initialize(
    enumeration_t const &enumeration, SitePermutation const &action,
    arma::Col<coeff_t> const &characters,
    bits::BitVector<typename enumeration_t::bit_t> &representative,
    bits::BitVector<typename enumeration_t::bit_t> &representative_index,
    bits::BitVector<typename enumeration_t::bit_t> &representative_symmetry,
    bits::BitVector<typename enumeration_t::bit_t> &representative_norm_index,
    std::vector<double> norms) try {
  using bit_t = typename enumeration_t::bit_t;
  using bits::BitVector;

  // First pass, simply count number of representatives and number of different
  // norms
  int64_t nrepresentatives = 0;
  for (auto state : enumeration) {
    if (isrepresentative(state, action)) {
      double nrm = norm(state, action, characters);

      if (std::fabs(nrm) > 1e-6) { // representative found
        ++nrepresentatives;

        // Check if norm already registered ...
        auto it = std::find_if(norms.begin(), norms.end(), [&](double n) {
          return std::fabs(n - nrm) < 1e-6;
        });

        // ... if not, register it
        if (it == norms.end()) {
          norms.push_back(nrm);
        }
      }
    }
  }
  std::sort(norms.begin(), norms.end());

  // Create vector holding the indices for each state yielding the
  // representative
  try {
    int64_t size = enumeration.size();
    int64_t nbits = bits::ceillog2(nrepresentatives);
    representative_index = BitVector<bit_t>(size, nbits);
  } catch (...) {
    XDIAG_THROW("Unable to allocate representative index array");
  }

  // Create vector holding the representatives and norm_index
  try {
    int64_t size = nrepresentatives;
    int64_t nbits = enumeration.bitwidth();
    representative = BitVector<bit_t>(size, nbits);
  } catch (...) {
    XDIAG_THROW("Unable to allocate representative array");
  }

  try {
    int64_t size = nrepresentatives;
    int64_t nbits = bits::ceillog2(norm.size());
    representative_norm_index = BitVector<bit_t>(size, nbits);
  } catch (...) {
    XDIAG_THROW("Unable to allocate representative norm index array");
  }

  // Second pass, fill all vectors
  int64_t idx = 0;
  nrepresentatives = 0;
  for (auto state : enumeration) {
    if (isrepresentative(state, action)) {
      double nrm = norm(state, action, characters);

      if (std::fabs(nrm) > 1e-6) { // representative found
        representative_index[idx] = nrepresentatives;
        representative[nrepresentatives] = state;

        // Get index of norm
        auto it = std::find_if(norms.begin(), norms.end(), [&](double n) {
          return std::fabs(n - nrm) < 1e-6;
        });
        representative_norm_index[nrepresentatives] =
            std::distance(it, norms.begin());
        ++nrepresentatives;
      }
    }
    ++idx;
  }

  // Compute the symmetries that yield the representative and fill
  // non-representative representative_index
  int64_t nbits_for_symmetry = bits::ceillog2(action.size());
  try {
    representative_symmetry =
        BitVector<bit_t>(nbits_for_symmetry, nrepresentatives);
  } catch (...) {
    XDIAG_THROW("Unable to allocate representative symmetry array");
  }

  auto const &group = action.group();
  for (int64_t rep_idx = 0; rep_idx < representative.size(); ++rep_idx) {
    bit_t rep = representative[rep_idx];
    for (int64_t sym = 0; sym < action.size(); ++sym) {
      bit_t state = action.apply(sym, rep);
      int64_t idx = enumeration.index(state);
      representative_index[idx] = rep_idx;
      representative_symmetry[idx] = group.inv(sym);
    }
  }
}
XDIAG_CATCH

template <typename enumeration_t>
RepresentativeTable<enumeration_t>::RepresentativeTable(
    enumeration_t const &enumeration, SitePermutation const &action,
    Representation const &irrep) try {
  if (action.group() != irrep.group()) {
    XDIAG_THROW("PermutationGroup of SitePermutation does not agree with group "
                "or Representation");
  }
  if (isreal(irrep)) {
    representative_table_initialize(
        enumeration, action, irrep.characters().as<arma::vec>(),
        representative_, representative_index_, representative_symmetry_,
        representative_norm_index_, norm_);
  } else {
    representative_table_initialize(
        enumeration, action, irrep.characters().as<arma::cx_vec>(),
        representative_, representative_index_, representative_symmetry_,
        representative_norm_index_, norm_);
  }
}
XDIAG_CATCH

template <typename enumeration_t>
int64_t RepresentativeTable<enumeration_t>::size() const {
  return representative_.size();
}

template <typename enumeration_t>
bool RepresentativeTable<enumeration_t>::operator==(
    RepresentativeTable<enumeration_t> const &rhs) const {
  return (representative_ == rhs.representative_) &&
         (representative_index_ == rhs.representative_index_) &&
         (representative_symmetry_ == rhs.representative_symmetry_) &&
         (representative_norm_index_ == rhs.representative_norm_index_) &&
         (norm_ == rhs.norm_);
}

template <typename enumeration_t>
bool RepresentativeTable<enumeration_t>::operator!=(
    RepresentativeTable<enumeration_t> const &rhs) const {
  return ~operator==(rhs);
}

using namespace combinatorics;
template class RepresentativeTable<Subsets<uint32_t>>;
template class RepresentativeTable<Subsets<uint64_t>>;
template class RepresentativeTable<Combinations<uint32_t>>;
template class RepresentativeTable<Combinations<uint64_t>>;

} // namespace xdiag::symmetries
