// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <cstdint>
#include <limits>
#include <string>

namespace xdiag::bits {

// Compact storage for multiple small integers using bit packing.
//
// BitArray stores multiple small integers (1-8 bits each) packed into
// a single native integer or Bitset. It provides fast get/set operations
// with specialized implementations for both native integers and Bitsets.
//
// Template parameters:
//   bit_t: Storage type - uint16_t, uint32_t, uint64_t, or Bitset<uint64_t, X>
//   nbits: Bits per element (1-8). Determines how many elements fit.
//
// Example:
//   BitArray<uint64_t, 3> arr;  // Stores 21 3-bit values (0-7 each)
//   arr.set(0, 5);              // Set first element to 5
//   int64_t val = arr.get(0);   // Get first element (returns 5)
//
// Storage capacity:
//   maximum_size = total_bits / nbits
//   For uint64_t with nbits=3: 64/3 = 21 elements
//   For Bitset<uint64_t, 2> with nbits=3: 128/3 = 42 elements
template <typename bit_tt, int nbitss> class BitArray {
public:
  using bit_t = bit_tt;
  static constexpr int nbits = nbitss;

  // Maximum number of elements that fit in the storage
  static constexpr int64_t maximum_size =
      std::numeric_limits<bit_t>::digits / nbits;

  // Default constructor, initializes all bits to 0
  BitArray() = default;

  // Get element at index (returns nbits-bit value)
  int64_t get(int64_t idx) const noexcept;

  // Set element at index to value (stores lower nbits of value)
  void set(int64_t idx, int64_t value) noexcept;

  bool operator==(BitArray<bit_t, nbits> const &rhs) const noexcept;
  bool operator!=(BitArray<bit_t, nbits> const &rhs) const noexcept;

private:
  bit_t bits_{};
};

template <typename bit_t, int nbits>
std::string to_string(BitArray<bit_t, nbits> const &bits,
                      int64_t size = BitArray<bit_t, nbits>::maximum_size,
                      bool reverse = true);

template <typename bit_t, int nbits>
std::ostream &operator<<(std::ostream &out, BitArray<bit_t, nbits> const &bits);

} // namespace xdiag::bits
