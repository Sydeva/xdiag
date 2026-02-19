// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

namespace xdiag::combinatorics {

template <typename bitarray_t> class BoundedMultisetsIterator;

// BoundedMultisets<bitarray_t> enumerates all ordered sequences of length n
// with elements drawn from {0, ..., bound-1}. Each sequence is packed into a
// bitarray_t = BitArray<bit_t, nbits> using nbits bits per slot (runtime bound,
// compile-time packing width). The constructor checks that
// ceillog2(bound) <= nbits and n <= bitarray_t::maximum_size.
// Sequences are produced in little-endian base-bound order (slot 0 is the
// least significant digit); total count is bound^n.
// Requires: n >= 0, bound >= 2.
//
// Example:
//   using A = BitArray<uint64_t, 2>;     // 2-bit slots -> bound up to 4
//   BoundedMultisets<A> ms(3, 3);        // 27 triples from {0, 1, 2}
//   for (auto seq : ms)
//     use seq.get(0), seq.get(1), seq.get(2);
template <typename bitarray_t> class BoundedMultisets {
public:
  using bit_t = typename bitarray_t::bit_t;
  static constexpr int nbits = bitarray_t::nbits;
  using iterator_t = BoundedMultisetsIterator<bitarray_t>;

  BoundedMultisets() = default;
  BoundedMultisets(int64_t n, int64_t bound);

  int64_t n() const;
  int64_t bound() const;
  int64_t size() const;
  iterator_t begin() const;
  iterator_t end() const;

  bool operator==(BoundedMultisets<bitarray_t> const &rhs) const;
  bool operator!=(BoundedMultisets<bitarray_t> const &rhs) const;

private:
  int64_t n_ = 0;
  int64_t bound_ = 0;
  int64_t size_ = 0;
};

template <typename bitarray_t> class BoundedMultisetsIterator {
public:
  using bit_t = typename bitarray_t::bit_t;
  static constexpr int nbits = bitarray_t::nbits;

  BoundedMultisetsIterator() = default;
  BoundedMultisetsIterator(int64_t n, int64_t idx, int64_t bound);

  bool operator==(BoundedMultisetsIterator<bitarray_t> const &rhs) const;
  bool operator!=(BoundedMultisetsIterator<bitarray_t> const &rhs) const;
  BoundedMultisetsIterator &operator++();
  bitarray_t operator*() const;

private:
  int64_t idx_;
  int64_t bound_;
};

} // namespace xdiag::combinatorics
