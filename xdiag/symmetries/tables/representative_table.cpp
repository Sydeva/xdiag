// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "representative_table.hpp"

#include <algorithm>
#include <cstdint>
#include <type_traits>

#include <xdiag/bits/bitarray.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/bits/log2.hpp>
#include <xdiag/bits/nbits.hpp>
#include <xdiag/combinatorics/bounded_multisets/bounded_multisets.hpp>
#include <xdiag/combinatorics/bounded_partitions/bounded_partitions.hpp>
#include <xdiag/combinatorics/bounded_partitions/schaefer_table.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
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
    bits::BitVector<uint64_t> &representative_index,
    bits::BitVector<uint64_t> &representative_symmetry,
    bits::BitVector<uint64_t> &representative_norm_index,
    std::vector<double> &norms) try {
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

  // --------------------------------------------------------------
  // Now come the allocations, since we know al the relevant sizes
  // --------------------------------------------------------------

  // Create vector holding the representatives
  try {
    int64_t size = nrepresentatives;
    int64_t nbits = enumeration.bitwidth();
    representative = BitVector<bit_t>(size, nbits);
  } catch (...) {
    XDIAG_THROW("Unable to allocate representative array");
  }

  // Create vector holding the indices for each state yielding the
  // representative
  try {
    int64_t size = enumeration.size();
    int64_t nbits = std::max(1u, bits::ceillog2(nrepresentatives));
    representative_index = BitVector<uint64_t>(size, nbits);
  } catch (...) {
    XDIAG_THROW("Unable to allocate representative index array");
  }

  // Create vector holding the symmetry which yields the representative
  try {
    int64_t size = enumeration.size();
    int64_t nbits = std::max(1u, bits::ceillog2(action.size()));
    representative_symmetry = BitVector<uint64_t>(size, nbits);
  } catch (...) {
    XDIAG_THROW("Unable to allocate representative symmetry array");
  }

  // Create vector holding the norm index of the states
  try {
    int64_t size = nrepresentatives;
    int64_t nbits = std::max(1u, bits::ceillog2(norms.size()));
    representative_norm_index = BitVector<uint64_t>(size, nbits);
  } catch (...) {
    XDIAG_THROW("Unable to allocate representative norm index array");
  }

  // --------------------------------------------------------------
  // Second pass, fill all vectors
  // --------------------------------------------------------------
  int64_t idx = 0;
  nrepresentatives = 0;
  for (auto state : enumeration) {
    if (isrepresentative(state, action)) {
      double nrm = norm(state, action, characters);

      if (std::fabs(nrm) > 1e-6) { // representative found
        representative_index[idx] = (uint64_t)nrepresentatives;
        representative[nrepresentatives] = state;

        // Get index of norm
        auto it = std::find_if(norms.begin(), norms.end(), [&](double n) {
          return std::fabs(n - nrm) < 1e-6;
        });
        representative_norm_index[nrepresentatives] =
            (uint64_t)std::distance(norms.begin(), it);
        ++nrepresentatives;
      }
    }
    ++idx;
  }

  // Compute the symmetries that yield the representative and fill
  // non-representative representative_index
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
  return !operator==(rhs);
}

} // namespace xdiag::symmetries

using namespace xdiag::combinatorics;
using namespace xdiag::bits;

#define INSTANTIATE_REPRESENTATIVE_TABLE(ENUMERATION_TYPE)                     \
  template class xdiag::symmetries::RepresentativeTable<ENUMERATION_TYPE>;


// BEGIN_INSTANTIATION_GROUP(subsets)
INSTANTIATE_REPRESENTATIVE_TABLE(Subsets<uint32_t>);
INSTANTIATE_REPRESENTATIVE_TABLE(Subsets<uint64_t>);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(combinations)
INSTANTIATE_REPRESENTATIVE_TABLE(Combinations<uint32_t>);
INSTANTIATE_REPRESENTATIVE_TABLE(Combinations<uint64_t>);
INSTANTIATE_REPRESENTATIVE_TABLE(Combinations<BitsetStatic2>);
INSTANTIATE_REPRESENTATIVE_TABLE(Combinations<BitsetStatic4>);
INSTANTIATE_REPRESENTATIVE_TABLE(Combinations<BitsetStatic8>);
INSTANTIATE_REPRESENTATIVE_TABLE(Combinations<BitsetDynamic>);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(lintable)
INSTANTIATE_REPRESENTATIVE_TABLE(LinTable<uint32_t>);
INSTANTIATE_REPRESENTATIVE_TABLE(LinTable<uint64_t>);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(bounded_multisets)
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArray1>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArray2>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArray3>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArray4>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArray5>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArray6>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArray7>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArray8>);

INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArrayLong1>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArrayLong2>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArrayLong3>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArrayLong4>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArrayLong5>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArrayLong6>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArrayLong7>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedMultisets<BitArrayLong8>);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(bounded_partitions)
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArray1>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArray2>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArray3>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArray4>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArray5>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArray6>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArray7>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArray8>);

INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArrayLong1>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArrayLong2>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArrayLong3>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArrayLong4>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArrayLong5>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArrayLong6>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArrayLong7>);
INSTANTIATE_REPRESENTATIVE_TABLE(BoundedPartitions<BitArrayLong8>);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(schaefer_table)
INSTANTIATE_REPRESENTATIVE_TABLE(SchaeferTable<BitArray1>);
INSTANTIATE_REPRESENTATIVE_TABLE(SchaeferTable<BitArray2>);
INSTANTIATE_REPRESENTATIVE_TABLE(SchaeferTable<BitArray3>);
INSTANTIATE_REPRESENTATIVE_TABLE(SchaeferTable<BitArray4>);
INSTANTIATE_REPRESENTATIVE_TABLE(SchaeferTable<BitArray5>);
INSTANTIATE_REPRESENTATIVE_TABLE(SchaeferTable<BitArray6>);
INSTANTIATE_REPRESENTATIVE_TABLE(SchaeferTable<BitArray7>);
INSTANTIATE_REPRESENTATIVE_TABLE(SchaeferTable<BitArray8>);
// END_INSTANTIATION_GROUP

#undef INSTANTIATE_REPRESENTATIVE_TABLE
