// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "bitarray.hpp"

#include <type_traits>
#include <xdiag/bits/bitmask.hpp>
#include <xdiag/bits/bitset.hpp>

namespace xdiag::bits {

template <typename bit_t, int nbits>
int64_t BitArray<bit_t, nbits>::get(int64_t idx) const noexcept {
  if constexpr (std::is_integral<bit_t>::value) {
    return (bits_ >> (idx * nbits)) & bitmask<bit_t>(nbits);
  } else {
    return bits_.get_range(idx * nbits, nbits);
  }
}

template <typename bit_t, int nbits>
void BitArray<bit_t, nbits>::set(int64_t idx, int64_t value) noexcept {
  if constexpr (std::is_integral<bit_t>::value) {
    int64_t shift = idx * nbits;
    bit_t mask = bitmask<bit_t>(nbits) << shift;
    bits_ &= ~mask;                     // clear bits
    bits_ |= ((value << shift) & mask); // set bits
  } else {
    bits_.set_range(idx * nbits, nbits, (typename bit_t::chunk_t)value);
  }
}

template <typename bit_t, int nbits>
bool BitArray<bit_t, nbits>::operator==(
    BitArray<bit_t, nbits> const &rhs) const noexcept {
  return bits_ == rhs.bits_;
}

template <typename bit_t, int nbits>
bool BitArray<bit_t, nbits>::operator!=(
    BitArray<bit_t, nbits> const &rhs) const noexcept {
  return !operator==(rhs);
}

template <typename bit_t, int nbits>
std::string to_string(BitArray<bit_t, nbits> const &bits, int64_t size) {
  std::string str;
  for (int64_t i = 0; i < size; ++i) {
    str += std::to_string(bits.get(i)) + std::string(" ");
  }
  return str;
}

template <typename bit_t, int nbits>
std::ostream &operator<<(std::ostream &out,
                         BitArray<bit_t, nbits> const &bits) {
  out << to_string(bits);
  return out;
}

// Template instantiations
#define INSTANTIATE_BITARRAY(BIT_T, NBITS)                                     \
  template class BitArray<BIT_T, NBITS>;                                       \
  template std::string to_string(BitArray<BIT_T, NBITS> const &, int64_t);     \
  template std::ostream &operator<<(std::ostream &,                            \
                                    BitArray<BIT_T, NBITS> const &);

#define INSTANTIATE_BITARRAY_FOR_NBITS(BIT_T)                                  \
  INSTANTIATE_BITARRAY(BIT_T, 1)                                               \
  INSTANTIATE_BITARRAY(BIT_T, 2)                                               \
  INSTANTIATE_BITARRAY(BIT_T, 3)                                               \
  INSTANTIATE_BITARRAY(BIT_T, 4)                                               \
  INSTANTIATE_BITARRAY(BIT_T, 5)                                               \
  INSTANTIATE_BITARRAY(BIT_T, 6)                                               \
  INSTANTIATE_BITARRAY(BIT_T, 7)                                               \
  INSTANTIATE_BITARRAY(BIT_T, 8)

// Native integer types
INSTANTIATE_BITARRAY_FOR_NBITS(uint16_t)
INSTANTIATE_BITARRAY_FOR_NBITS(uint32_t)
INSTANTIATE_BITARRAY_FOR_NBITS(uint64_t)

// Bitset types
INSTANTIATE_BITARRAY_FOR_NBITS(BitsetDynamic)
INSTANTIATE_BITARRAY_FOR_NBITS(BitsetStatic1)
INSTANTIATE_BITARRAY_FOR_NBITS(BitsetStatic2)
INSTANTIATE_BITARRAY_FOR_NBITS(BitsetStatic4)
INSTANTIATE_BITARRAY_FOR_NBITS(BitsetStatic8)

#undef INSTANTIATE_BITARRAY_FOR_NBITS
#undef INSTANTIATE_BITARRAY

} // namespace xdiag::bits
