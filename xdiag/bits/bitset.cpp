// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "bitset.hpp"
#include <bitset>
#include <cassert>

#include <xdiag/bits/bitmask.hpp>
#include <xdiag/bits/gbit.hpp>
#include <xdiag/bits/popcnt.hpp>
#include <xdiag/utils/logger.hpp>

#include <xdiag/bits/log2.hpp>

namespace xdiag::bits {

template <typename chunk_t>
static constexpr int64_t n_chunks_for_bits(int64_t nbits) {
  constexpr size_t nchunkbits = std::numeric_limits<chunk_t>::digits;
  constexpr size_t chunkshift = floorlog2(nchunkbits);
  return nbits > 0 ? ((nbits - 1) >> chunkshift) + 1 : 0;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks>::Bitset(int64_t nbits)
    : chunks_([&]() {
        if constexpr (nchunks == 0) {
          // Dynamic: construct vector with calculated size (value-initialized
          // to 0)
          return storage_t(n_chunks_for_bits<chunk_t>(nbits));
        } else {
          // Static: value-initialize array (all elements to 0)
          return storage_t{};
        }
      }()) {
  // For static storage, verify nbits fits in nchunks
  if constexpr (nchunks > 0) {
    assert(n_chunks_for_bits<chunk_t>(nbits) <= nchunks);
  }
}

template <typename chunk_t, int64_t nchunks>
typename Bitset<chunk_t, nchunks>::storage_t const &
Bitset<chunk_t, nchunks>::chunks() const noexcept {
  return chunks_;
}

// Bit-level access
template <typename chunk_t, int64_t nchunks>
bool Bitset<chunk_t, nchunks>::test(int64_t pos) const noexcept {
  int64_t chunk_idx = pos >> chunkshift_;
  int64_t bit_idx = pos & chunkmask_;
  return gbit(chunks_[chunk_idx], bit_idx);
}

template <typename chunk_t, int64_t nchunks>
void Bitset<chunk_t, nchunks>::set(int64_t pos) noexcept {
  int64_t chunk_idx = pos >> chunkshift_;
  int64_t bit_idx = pos & chunkmask_;
  chunks_[chunk_idx] |= ((chunk_t)1 << bit_idx);
}

template <typename chunk_t, int64_t nchunks>
void Bitset<chunk_t, nchunks>::set(int64_t pos, bool value) noexcept {
  int64_t chunk_idx = pos >> chunkshift_;
  int64_t bit_idx = pos & chunkmask_;
  if (value) {
    chunks_[chunk_idx] |= ((chunk_t)1 << bit_idx);
  } else {
    chunks_[chunk_idx] &= ~((chunk_t)1 << bit_idx);
  }
}

template <typename chunk_t, int64_t nchunks>
void Bitset<chunk_t, nchunks>::reset(int64_t pos) noexcept {
  int64_t chunk_idx = pos >> chunkshift_;
  int64_t bit_idx = pos & chunkmask_;
  chunks_[chunk_idx] &= ~((chunk_t)1 << bit_idx);
}

template <typename chunk_t, int64_t nchunks>
void Bitset<chunk_t, nchunks>::flip(int64_t pos) noexcept {
  int64_t chunk_idx = pos >> chunkshift_;
  int64_t bit_idx = pos & chunkmask_;
  chunks_[chunk_idx] ^= ((chunk_t)1 << bit_idx);
}

// Bit-level access (ranged)
template <typename chunk_t, int64_t nchunks>
void Bitset<chunk_t, nchunks>::set_range(int64_t start, int64_t length,
                                         chunk_t bits) noexcept {
  assert(length <= nchunkbits_);
  if (!length) {
    return;
  }
  int64_t end = start + length;
  int64_t startchunk = start >> chunkshift_; // divide by nchunkbits
  int64_t startbit = start & chunkmask_;     // modulo    nchunkbits
  int64_t endchunk = end >> chunkshift_;
  int64_t endbit = end & chunkmask_;
  if ((endchunk == startchunk) || (endbit == 0)) {
    chunk_t mask = bitmask<chunk_t>(length) << startbit;
    chunks_[startchunk] &= ~mask;
    chunks_[startchunk] |= bits << startbit;
  } else {
    chunk_t negmask1 = bitmask<chunk_t>(startbit);
    chunks_[startchunk] &= negmask1;
    chunks_[startchunk] |= bits << startbit;
    chunk_t mask2 = bitmask<chunk_t>(endbit);
    chunks_[endchunk] &= ~mask2;
    chunks_[endchunk] |= bits >> (nchunkbits_ - startbit);
  }
}

template <typename chunk_t, int64_t nchunks>
chunk_t Bitset<chunk_t, nchunks>::get_range(int64_t start,
                                            int64_t length) const noexcept {
  assert(length <= nchunkbits_);
  if (!length) {
    return (chunk_t)0;
  }
  int64_t end = start + length;
  int64_t startchunk = start >> chunkshift_; // divide by nchunkbits
  int64_t startbit = start & chunkmask_;     // modulo    nchunkbits
  int64_t endchunk = end >> chunkshift_;
  int64_t endbit = end & chunkmask_;
  if ((endchunk == startchunk) || (endbit == 0)) {
    return (chunks_[startchunk] >> startbit) & bitmask<chunk_t>(length);
  } else {
    return ((chunks_[endchunk] & bitmask<chunk_t>(endbit))
            << (nchunkbits_ - startbit)) |
           (chunks_[startchunk] >> startbit);
  }
}

// Bitwise operations
template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks>
Bitset<chunk_t, nchunks>::operator&(Bitset const &rhs) const {
  Bitset result = *this;
  result &= rhs;
  return result;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks>
Bitset<chunk_t, nchunks>::operator|(Bitset const &rhs) const {
  Bitset result = *this;
  result |= rhs;
  return result;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks>
Bitset<chunk_t, nchunks>::operator^(Bitset const &rhs) const {
  Bitset result = *this;
  result ^= rhs;
  return result;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks> Bitset<chunk_t, nchunks>::operator~() const {
  Bitset result = *this;
  for (int64_t i = 0; i < std::size(result.chunks_); ++i) {
    result.chunks_[i] = ~result.chunks_[i];
  }
  return result;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks> &
Bitset<chunk_t, nchunks>::operator&=(Bitset const &rhs) noexcept {
  assert(std::size(chunks_) == std::size(rhs.chunks_));
  for (int64_t i = 0; i < std::size(chunks_); ++i) {
    chunks_[i] &= rhs.chunks_[i];
  }
  return *this;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks> &
Bitset<chunk_t, nchunks>::operator|=(Bitset const &rhs) noexcept {
  assert(std::size(chunks_) == std::size(rhs.chunks_));
  for (int64_t i = 0; i < std::size(chunks_); ++i) {
    chunks_[i] |= rhs.chunks_[i];
  }
  return *this;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks> &
Bitset<chunk_t, nchunks>::operator^=(Bitset const &rhs) noexcept {
  assert(std::size(chunks_) == std::size(rhs.chunks_));
  for (int64_t i = 0; i < std::size(chunks_); ++i) {
    chunks_[i] ^= rhs.chunks_[i];
  }
  return *this;
}

// Optimized shift-by-1 for division algorithm (private helper)
template <typename chunk_t, int64_t nchunks>
void Bitset<chunk_t, nchunks>::shift_left_by_1() noexcept {
  int64_t size = std::size(chunks_);
  chunk_t carry = 0;
  for (int64_t i = 0; i < size; ++i) {
    chunk_t next_carry = chunks_[i] >> (nchunkbits_ - 1);
    chunks_[i] = (chunks_[i] << 1) | carry;
    carry = next_carry;
  }
}

// Shift operations
template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks>
Bitset<chunk_t, nchunks>::operator<<(int64_t shift) const {
  Bitset result = *this;
  result <<= shift;
  return result;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks>
Bitset<chunk_t, nchunks>::operator>>(int64_t shift) const {
  Bitset result = *this;
  result >>= shift;
  return result;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks> &
Bitset<chunk_t, nchunks>::operator<<=(int64_t shift) noexcept {
  if (shift == 0)
    return *this;

  int64_t size = std::size(chunks_);
  int64_t chunk_shift = shift >> chunkshift_;
  int64_t bit_shift = shift & chunkmask_;

  if (chunk_shift >= size) {
    for (auto &chunk : chunks_) {
      chunk = 0;
    }
    return *this;
  }

  if (bit_shift == 0) {
    // Simple chunk-aligned shift
    for (int64_t i = size - 1; i >= chunk_shift; --i) {
      chunks_[i] = chunks_[i - chunk_shift];
    }
    for (int64_t i = 0; i < chunk_shift; ++i) {
      chunks_[i] = 0;
    }
  } else {
    // Shift with bit offset
    int64_t complement_shift = nchunkbits_ - bit_shift;
    for (int64_t i = size - 1; i > chunk_shift; --i) {
      chunks_[i] = (chunks_[i - chunk_shift] << bit_shift) |
                   (chunks_[i - chunk_shift - 1] >> complement_shift);
    }
    chunks_[chunk_shift] = chunks_[0] << bit_shift;
    for (int64_t i = 0; i < chunk_shift; ++i) {
      chunks_[i] = 0;
    }
  }
  return *this;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks> &
Bitset<chunk_t, nchunks>::operator>>=(int64_t shift) noexcept {
  if (shift == 0)
    return *this;

  int64_t size = std::size(chunks_);
  int64_t chunk_shift = shift >> chunkshift_;
  int64_t bit_shift = shift & chunkmask_;

  if (chunk_shift >= size) {
    for (auto &chunk : chunks_) {
      chunk = 0;
    }
    return *this;
  }

  if (bit_shift == 0) {
    // Simple chunk-aligned shift
    for (int64_t i = 0; i < size - chunk_shift; ++i) {
      chunks_[i] = chunks_[i + chunk_shift];
    }
    for (int64_t i = size - chunk_shift; i < size; ++i) {
      chunks_[i] = 0;
    }
  } else {
    // Shift with bit offset
    int64_t complement_shift = nchunkbits_ - bit_shift;
    for (int64_t i = 0; i < size - chunk_shift - 1; ++i) {
      chunks_[i] = (chunks_[i + chunk_shift] >> bit_shift) |
                   (chunks_[i + chunk_shift + 1] << complement_shift);
    }
    chunks_[size - chunk_shift - 1] = chunks_[size - 1] >> bit_shift;
    for (int64_t i = size - chunk_shift; i < size; ++i) {
      chunks_[i] = 0;
    }
  }
  return *this;
}

// Arithmetic operations
template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks>
Bitset<chunk_t, nchunks>::operator+(Bitset const &rhs) const {
  Bitset result = *this;
  result += rhs;
  return result;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks>
Bitset<chunk_t, nchunks>::operator-(Bitset const &rhs) const {
  Bitset result = *this;
  result -= rhs;
  return result;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks> Bitset<chunk_t, nchunks>::operator-() const {
  Bitset result = *this;
  // Two's complement: flip all bits and add 1 (correct for unsigned modular
  // arithmetic)
  for (int64_t i = 0; i < std::size(result.chunks_); ++i) {
    result.chunks_[i] = ~result.chunks_[i];
  }
  // Add 1
  chunk_t carry = 1;
  for (int64_t i = 0; i < std::size(result.chunks_) && carry; ++i) {
    chunk_t old_val = result.chunks_[i];
    result.chunks_[i] = old_val + carry;
    carry = (result.chunks_[i] < old_val) ? 1 : 0;
  }
  return result;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks>
Bitset<chunk_t, nchunks>::operator*(Bitset const &rhs) const {
  Bitset result = *this;
  result *= rhs;
  return result;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks>
Bitset<chunk_t, nchunks>::operator/(Bitset const &rhs) const {
  Bitset result = *this;
  result /= rhs;
  return result;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks> &
Bitset<chunk_t, nchunks>::operator+=(Bitset const &rhs) noexcept {
  assert(std::size(chunks_) == std::size(rhs.chunks_));
  chunk_t carry = 0;
  for (int64_t i = 0; i < std::size(chunks_); ++i) {
    chunk_t old_chunk = chunks_[i];
    chunks_[i] += rhs.chunks_[i] + carry;
    // Detect overflow: if result is less than either operand, we had overflow
    carry =
        (chunks_[i] < old_chunk || (carry && chunks_[i] == old_chunk)) ? 1 : 0;
  }
  return *this;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks> &
Bitset<chunk_t, nchunks>::operator-=(Bitset const &rhs) noexcept {
  assert(std::size(chunks_) == std::size(rhs.chunks_));
  chunk_t borrow = 0;
  for (int64_t i = 0; i < std::size(chunks_); ++i) {
    chunk_t old_chunk = chunks_[i];
    chunks_[i] -= rhs.chunks_[i] + borrow;
    // Detect underflow: if result is greater than original, we had underflow
    borrow =
        (chunks_[i] > old_chunk || (borrow && chunks_[i] == old_chunk)) ? 1 : 0;
  }
  return *this;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks> &
Bitset<chunk_t, nchunks>::operator*=(Bitset const &rhs) noexcept {
  assert(std::size(chunks_) == std::size(rhs.chunks_));
  Bitset<chunk_t, nchunks> result(std::size(chunks_) * nchunkbits_);

  // Grade-school multiplication without extended types
  // Split each chunk multiplication into high and low parts
  constexpr int half_bits = nchunkbits_ / 2;
  constexpr chunk_t low_mask = (static_cast<chunk_t>(1) << half_bits) - 1;

  for (int64_t i = 0; i < std::size(chunks_); ++i) {
    if (rhs.chunks_[i] == 0)
      continue;

    chunk_t a_lo = rhs.chunks_[i] & low_mask;
    chunk_t a_hi = rhs.chunks_[i] >> half_bits;

    chunk_t carry = 0;
    for (int64_t j = 0; j < std::size(chunks_) - i; ++j) {
      chunk_t b_lo = chunks_[j] & low_mask;
      chunk_t b_hi = chunks_[j] >> half_bits;

      // Compute a * b as (a_hi * 2^k + a_lo) * (b_hi * 2^k + b_lo)
      chunk_t p_ll = a_lo * b_lo;
      chunk_t p_lh = a_lo * b_hi;
      chunk_t p_hl = a_hi * b_lo;
      chunk_t p_hh = a_hi * b_hi;

      // Combine: p_ll + (p_lh + p_hl) * 2^k + p_hh * 2^(2k)
      chunk_t mid = p_lh + p_hl;
      chunk_t mid_carry = (mid < p_lh) ? 1 : 0;

      chunk_t low = p_ll + ((mid & low_mask) << half_bits);
      chunk_t low_carry = (low < p_ll) ? 1 : 0;

      chunk_t high =
          p_hh + (mid >> half_bits) + (mid_carry << half_bits) + low_carry;

      // Add to result with carry
      chunk_t old_val = result.chunks_[i + j];
      chunk_t temp = old_val + low;
      chunk_t add_carry1 = (temp < old_val) ? 1 : 0;
      result.chunks_[i + j] = temp + carry;
      chunk_t add_carry2 = (result.chunks_[i + j] < temp) ? 1 : 0;
      carry = high + add_carry1 + add_carry2;
    }
  }

  *this = result;
  return *this;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks> &
Bitset<chunk_t, nchunks>::operator/=(Bitset const &rhs) noexcept {
  assert(std::size(chunks_) == std::size(rhs.chunks_));

  // Division by zero check
  if (rhs.none()) {
    // Division by zero - return 0
    for (auto &chunk : chunks_) {
      chunk = 0;
    }
    return *this;
  }

  // Long division algorithm for multi-precision integers
  Bitset<chunk_t, nchunks> quotient(std::size(chunks_) * nchunkbits_);
  Bitset<chunk_t, nchunks> remainder(std::size(chunks_) * nchunkbits_);

  // Find the most significant bit position
  int64_t nbits = std::size(chunks_) * nchunkbits_;

  for (int64_t i = nbits - 1; i >= 0; --i) {
    // Shift remainder left by 1 (optimized)
    remainder.shift_left_by_1();
    // Set the least significant bit of remainder to bit i of dividend
    if (test(i)) {
      remainder.set(0);
    }

    // If remainder >= divisor, subtract divisor and set quotient bit
    if (remainder >= rhs) {
      remainder -= rhs;
      quotient.set(i);
    }
  }

  *this = quotient;
  return *this;
}

// Predicates
template <typename chunk_t, int64_t nchunks>
bool Bitset<chunk_t, nchunks>::all() const noexcept {
  for (auto chunk : chunks_) {
    if (chunk != std::numeric_limits<chunk_t>::max()) {
      return false;
    }
  }
  return true;
}

template <typename chunk_t, int64_t nchunks>
bool Bitset<chunk_t, nchunks>::any() const noexcept {
  for (auto chunk : chunks_) {
    if (chunk != 0) {
      return true;
    }
  }
  return false;
}

template <typename chunk_t, int64_t nchunks>
bool Bitset<chunk_t, nchunks>::none() const noexcept {
  return !any();
}

template <typename chunk_t, int64_t nchunks>
int64_t Bitset<chunk_t, nchunks>::count() const noexcept {
  int64_t total = 0;
  for (auto chunk : chunks_) {
    total += popcnt(chunk);
  }
  return total;
}

template <typename chunk_t, int64_t nchunks>
bool Bitset<chunk_t, nchunks>::operator==(
    Bitset<chunk_t, nchunks> const &rhs) const noexcept {
  return chunks_ == rhs.chunks_;
}

template <typename chunk_t, int64_t nchunks>
bool Bitset<chunk_t, nchunks>::operator!=(
    Bitset<chunk_t, nchunks> const &rhs) const noexcept {
  return !operator==(rhs);
}

template <typename chunk_t, int64_t nchunks>
bool Bitset<chunk_t, nchunks>::operator<(
    Bitset<chunk_t, nchunks> const &rhs) const noexcept {
  assert(std::size(chunks_) == std::size(rhs.chunks_));
  // Compare from most significant chunk to least significant
  for (int64_t i = std::size(chunks_) - 1; i >= 0; --i) {
    if (chunks_[i] < rhs.chunks_[i]) {
      return true;
    }
    if (chunks_[i] > rhs.chunks_[i]) {
      return false;
    }
  }
  return false; // equal
}

template <typename chunk_t, int64_t nchunks>
bool Bitset<chunk_t, nchunks>::operator<=(
    Bitset<chunk_t, nchunks> const &rhs) const noexcept {
  return !operator>(rhs);
}

template <typename chunk_t, int64_t nchunks>
bool Bitset<chunk_t, nchunks>::operator>(
    Bitset<chunk_t, nchunks> const &rhs) const noexcept {
  return rhs.operator<(*this);
}

template <typename chunk_t, int64_t nchunks>
bool Bitset<chunk_t, nchunks>::operator>=(
    Bitset<chunk_t, nchunks> const &rhs) const noexcept {
  return !operator<(rhs);
}

template <typename chunk_t, int64_t nchunks>
std::string to_string(Bitset<chunk_t, nchunks> const &bits) {
  std::string str;
  for (auto const &chunk : bits.chunks()) {
    str = std::bitset<bits.nchunkbits()>(chunk).to_string() + str;
  }
  return str;
}

template <typename chunk_t, int64_t nchunks>
std::ostream &operator<<(std::ostream &out,
                         Bitset<chunk_t, nchunks> const &bits) {
  out << to_string(bits);
  return out;
}

template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks> make_bitset(uint64_t value) {
  Bitset<chunk_t, nchunks> bits(64);
  for (int i = 0; i < 64; ++i) {
    if (value & (1ULL << i)) {
      bits.set(i);
    }
  }
  return bits;
}

template <typename chunk_t, int64_t nchunks>
uint64_t to_uint64(Bitset<chunk_t, nchunks> const &bits) {
  uint64_t result = 0;
  constexpr int64_t chunk_bits = std::numeric_limits<chunk_t>::digits;

  int64_t max_bits;
  if constexpr (nchunks == 0) {
    // Dynamic: runtime calculation
    max_bits =
        std::min(int64_t(64), int64_t(std::size(bits.chunks()) * chunk_bits));
  } else {
    // Static: compile-time calculation
    constexpr int64_t max_bits_static =
        std::min(int64_t(64), nchunks * chunk_bits);
    max_bits = max_bits_static;
  }

  for (int i = 0; i < max_bits; ++i) {
    if (bits.test(i)) {
      result |= (1ULL << i);
    }
  }
  return result;
}

// Explicit template instantiations
#define INSTANTIATE_XDIAG_BITS_BITSET(CHUNK_T, NCHUNKS)                        \
  template class Bitset<CHUNK_T, NCHUNKS>;                                     \
  template std::string to_string(Bitset<CHUNK_T, NCHUNKS> const &);            \
  template std::ostream &operator<<(std::ostream &,                            \
                                    Bitset<CHUNK_T, NCHUNKS> const &);         \
  template Bitset<CHUNK_T, NCHUNKS> make_bitset(uint64_t);                     \
  template uint64_t to_uint64(Bitset<CHUNK_T, NCHUNKS> const &);

#define INSTANTIATE_XDIAG_BITS_BITSET_FOR_NCHUNKS(CHUNK_T)                     \
  INSTANTIATE_XDIAG_BITS_BITSET(CHUNK_T, 0)                                    \
  INSTANTIATE_XDIAG_BITS_BITSET(CHUNK_T, 1)                                    \
  INSTANTIATE_XDIAG_BITS_BITSET(CHUNK_T, 2)                                    \
  INSTANTIATE_XDIAG_BITS_BITSET(CHUNK_T, 4)                                    \
  INSTANTIATE_XDIAG_BITS_BITSET(CHUNK_T, 8)

INSTANTIATE_XDIAG_BITS_BITSET_FOR_NCHUNKS(uint8_t)
INSTANTIATE_XDIAG_BITS_BITSET_FOR_NCHUNKS(uint16_t)
INSTANTIATE_XDIAG_BITS_BITSET_FOR_NCHUNKS(uint32_t)
INSTANTIATE_XDIAG_BITS_BITSET_FOR_NCHUNKS(uint64_t)

#undef INSTANTIATE_XDIAG_BITS_BITSET_FOR_NCHUNKS
#undef INSTANTIATE_XDIAG_BITS_BITSET

} // namespace xdiag::bits
