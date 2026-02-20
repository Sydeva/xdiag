// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "bounded_partitions.hpp"

#include <algorithm>
#include <xdiag/bits/bitarray.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/bits/log2.hpp>
#include <xdiag/combinatorics/binomial.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::combinatorics {

// ---------------------------------------------------------------------------
// Count weak compositions of total into n parts each in {0, ..., bound-1}.
// Inclusion-exclusion: sum_{k=0}^{K} (-1)^k C(n,k) C(total-k*bound+n-1, n-1)
// ---------------------------------------------------------------------------
static int64_t count_bounded_partitions(int64_t n, int64_t total,
                                         int64_t bound) {
  if (n == 0)
    return (total == 0) ? 1 : 0;
  if (total < 0 || total > n * (bound - 1))
    return 0;
  int64_t result = 0;
  for (int64_t k = 0; k * bound <= total; ++k) {
    int64_t b1 = binom(n, k);
    int64_t b2 = binom(total - k * bound + n - 1, n - 1);
    if (k % 2 == 0)
      result += b1 * b2;
    else
      result -= b1 * b2;
  }
  return result;
}

// ---------------------------------------------------------------------------
// Fill slots lo..hi with the rlex-first valid configuration summing to S.
// Greedy right-to-left: a[hi] = max(0, S - (hi-lo)*(bound-1)), then a[hi-1],
// etc. This minimises a[hi] first, then a[hi-1], ..., leaving as much as
// possible to the left — i.e. the smallest config when slots are read
// right-to-left (consistent with the little-endian BoundedMultisets order).
// ---------------------------------------------------------------------------
template <typename bitarray_t>
static void fill_rlex_first(bitarray_t &arr, int64_t lo, int64_t hi,
                             int64_t S, int64_t bound) {
  int64_t remaining = S;
  for (int64_t j = hi; j > lo; --j) {
    int64_t val = std::max((int64_t)0, remaining - (j - lo) * (bound - 1));
    arr.set(j, val);
    remaining -= val;
  }
  arr.set(lo, remaining);
}

// ---------------------------------------------------------------------------
// BoundedPartitions
// ---------------------------------------------------------------------------

template <typename bitarray_t>
BoundedPartitions<bitarray_t>::BoundedPartitions(int64_t n, int64_t total,
                                                  int64_t bound) try
    : n_(n), total_(total), bound_(bound),
      size_(count_bounded_partitions(n, total, bound)) {
  if (n < 0) {
    XDIAG_THROW("Error constructing BoundedPartitions: n < 0");
  }
  if (total < 0) {
    XDIAG_THROW("Error constructing BoundedPartitions: total < 0");
  }
  if (bound < 2) {
    XDIAG_THROW("Error constructing BoundedPartitions: bound < 2");
  }
  int64_t nbits_for_bound = bits::ceillog2(bound);
  if (nbits_for_bound > nbits) {
    XDIAG_THROW(fmt::format("Error constructing BoundedPartitions: bound ({}) "
                            "is too large for BitArray type with nbits ({})",
                            bound, nbits));
  }
  if (n > bitarray_t::maximum_size) {
    XDIAG_THROW(fmt::format(
        "Error constructing BoundedPartitions: n ({}) "
        "is too large for BitArray type with maximum_size ({})",
        n, bitarray_t::maximum_size));
  }
}
XDIAG_CATCH

template <typename bitarray_t>
int64_t BoundedPartitions<bitarray_t>::n() const { return n_; }

template <typename bitarray_t>
int64_t BoundedPartitions<bitarray_t>::total() const { return total_; }

template <typename bitarray_t>
int64_t BoundedPartitions<bitarray_t>::bound() const { return bound_; }

template <typename bitarray_t>
int64_t BoundedPartitions<bitarray_t>::size() const { return size_; }

template <typename bitarray_t>
BoundedPartitionsIterator<bitarray_t>
BoundedPartitions<bitarray_t>::begin() const {
  if (size_ == 0)
    return end();
  bitarray_t first{};
  if (n_ > 0)
    fill_rlex_first(first, 0, n_ - 1, total_, bound_);
  return BoundedPartitionsIterator<bitarray_t>(n_, bound_, 0, first);
}

template <typename bitarray_t>
BoundedPartitionsIterator<bitarray_t>
BoundedPartitions<bitarray_t>::end() const {
  return BoundedPartitionsIterator<bitarray_t>(n_, bound_, size_, bitarray_t{});
}

template <typename bitarray_t>
bool BoundedPartitions<bitarray_t>::operator==(
    BoundedPartitions<bitarray_t> const &rhs) const {
  return (n_ == rhs.n_) && (total_ == rhs.total_) && (bound_ == rhs.bound_);
}

template <typename bitarray_t>
bool BoundedPartitions<bitarray_t>::operator!=(
    BoundedPartitions<bitarray_t> const &rhs) const {
  return !operator==(rhs);
}

// ---------------------------------------------------------------------------
// BoundedPartitionsIterator
// ---------------------------------------------------------------------------

template <typename bitarray_t>
BoundedPartitionsIterator<bitarray_t>::BoundedPartitionsIterator(
    int64_t n, int64_t bound, int64_t idx, bitarray_t current)
    : n_(n), bound_(bound), idx_(idx), current_(current) {}

template <typename bitarray_t>
bool BoundedPartitionsIterator<bitarray_t>::operator==(
    BoundedPartitionsIterator<bitarray_t> const &rhs) const {
  return idx_ == rhs.idx_;
}

template <typename bitarray_t>
bool BoundedPartitionsIterator<bitarray_t>::operator!=(
    BoundedPartitionsIterator<bitarray_t> const &rhs) const {
  return idx_ != rhs.idx_;
}

// Advance to the next sequence in reverse-lex order (slot n-1 most significant,
// slot 0 least significant — consistent with BoundedMultisets).
// Scan from i = 1 up to n-1: at the leftmost position where a[i] can be
// incremented and the left-side sum a[0..i-1] reduced by 1 is achievable,
// do so and fill a[0..i-1] with the rlex-first config for the new left-sum.
template <typename bitarray_t>
BoundedPartitionsIterator<bitarray_t> &
BoundedPartitionsIterator<bitarray_t>::operator++() {
  for (int64_t i = 1; i < n_; ++i) {
    int64_t ai = current_.get(i);
    if (ai >= bound_ - 1)
      continue;
    // Sum of a[0..i-1]
    int64_t si_left = 0;
    for (int64_t j = 0; j < i; ++j)
      si_left += current_.get(j);
    // After incrementing a[i], left slots a[0..i-1] must sum to si_left - 1
    int64_t new_left = si_left - 1;
    if (new_left < 0 || new_left > i * (bound_ - 1))
      continue;
    current_.set(i, ai + 1);
    fill_rlex_first(current_, 0, i - 1, new_left, bound_);
    ++idx_;
    return *this;
  }
  ++idx_; // move to end
  return *this;
}

template <typename bitarray_t>
auto BoundedPartitionsIterator<bitarray_t>::operator*() const -> bitarray_t {
  return current_;
}

// ---------------------------------------------------------------------------
// Template instantiations
// ---------------------------------------------------------------------------
// clang-format off
#define INSTANTIATE_BP(BIT_T, NBITS)                                           \
  template class BoundedPartitions<bits::BitArray<BIT_T, NBITS>>;             \
  template class BoundedPartitionsIterator<bits::BitArray<BIT_T, NBITS>>;

#define INSTANTIATE_BP_ALL(BIT_T)                                              \
  INSTANTIATE_BP(BIT_T, 1)                                                     \
  INSTANTIATE_BP(BIT_T, 2)                                                     \
  INSTANTIATE_BP(BIT_T, 3)                                                     \
  INSTANTIATE_BP(BIT_T, 4)                                                     \
  INSTANTIATE_BP(BIT_T, 5)                                                     \
  INSTANTIATE_BP(BIT_T, 6)                                                     \
  INSTANTIATE_BP(BIT_T, 7)                                                     \
  INSTANTIATE_BP(BIT_T, 8)

INSTANTIATE_BP_ALL(uint16_t)
INSTANTIATE_BP_ALL(uint32_t)
INSTANTIATE_BP_ALL(uint64_t)
INSTANTIATE_BP_ALL(bits::BitsetDynamic)
INSTANTIATE_BP_ALL(bits::BitsetStatic1)
INSTANTIATE_BP_ALL(bits::BitsetStatic2)
INSTANTIATE_BP_ALL(bits::BitsetStatic4)
INSTANTIATE_BP_ALL(bits::BitsetStatic8)

#undef INSTANTIATE_BP_ALL
#undef INSTANTIATE_BP
// clang-format on

} // namespace xdiag::combinatorics
