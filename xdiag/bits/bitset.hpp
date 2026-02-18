// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <array>
#include <cstdint>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include <xdiag/bits/bitmask.hpp>
#include <xdiag/bits/log2.hpp>
#include <xdiag/common.hpp>

namespace xdiag::bits {

// Multi-precision bit storage with arithmetic operations.
//
// Bitset stores a sequence of bits using chunks (uint8_t/16/32/64) and supports
// both dynamic sizing (nchunks=0) and static sizing (nchunks>0). Provides bitwise
// operations, shifts, arithmetic (+,-,*,/), and comparisons, modeling unsigned integers.
//
// Template parameters:
//   chunk_tt: Chunk type (uint8_t, uint16_t, uint32_t, uint64_t)
//   nchunks: Number of chunks (0=dynamic using std::vector, >0=static using std::array)
//
// Example:
//   Bitset<uint64_t, 2> bits;      // Static: 128 bits (2 × 64)
//   bits.set(65);                  // Set bit 65 to 1
//   bits = bits + bits;            // Arithmetic operations
//   Bitset<uint64_t, 0> dynamic(200);  // Dynamic: 200 bits
template <typename chunk_tt = uint64_t, int64_t nchunks = 0> class Bitset {
public:
  using chunk_t = chunk_tt;

  // storage_t depends template parameter nchunks:
  // nchunks == 0 -> std::vector<chunk_t>          (dynamic)
  // nchunks != 0 -> std::array<chunk_t, nchunk>   (static)
  using storage_t =
      typename std::conditional<(bool)nchunks, std::array<chunk_t, nchunks>,
                                std::vector<chunk_t>>::type;

  Bitset() = default;
  explicit Bitset(int64_t nbits);

  // Bit-level access
  bool test(int64_t pos) const noexcept;
  void set(int64_t pos) noexcept;             // fast: always sets to 1
  void set(int64_t pos, bool value) noexcept; // conditional set
  void reset(int64_t pos) noexcept;
  void flip(int64_t pos) noexcept;

  // Bit-level access (ranged), assuming length <= nchunkbits_
  void set_range(int64_t start, int64_t length, chunk_t bits) noexcept;
  chunk_t get_range(int64_t start, int64_t length) const noexcept;

  // Bitwise operations
  Bitset operator&(Bitset const &rhs) const;
  Bitset operator|(Bitset const &rhs) const;
  Bitset operator^(Bitset const &rhs) const;
  Bitset operator~() const;
  Bitset &operator&=(Bitset const &rhs) noexcept;
  Bitset &operator|=(Bitset const &rhs) noexcept;
  Bitset &operator^=(Bitset const &rhs) noexcept;

  // Shift operations
  Bitset operator<<(int64_t shift) const;
  Bitset operator>>(int64_t shift) const;
  Bitset &operator<<=(int64_t shift) noexcept;
  Bitset &operator>>=(int64_t shift) noexcept;

  // Arithmetic operations
  Bitset operator+(Bitset const &rhs) const;
  Bitset operator-(Bitset const &rhs) const;
  Bitset operator-() const;
  Bitset operator*(Bitset const &rhs) const;
  Bitset operator/(Bitset const &rhs) const;
  Bitset &operator+=(Bitset const &rhs) noexcept;
  Bitset &operator-=(Bitset const &rhs) noexcept;
  Bitset &operator*=(Bitset const &rhs) noexcept;
  Bitset &operator/=(Bitset const &rhs) noexcept;

  // Predicates
  bool all() const noexcept;
  bool any() const noexcept;
  bool none() const noexcept;
  int64_t count() const noexcept;

  bool operator==(Bitset<chunk_t, nchunks> const &rhs) const noexcept;
  bool operator!=(Bitset<chunk_t, nchunks> const &rhs) const noexcept;
  bool operator<(Bitset<chunk_t, nchunks> const &rhs) const noexcept;
  bool operator<=(Bitset<chunk_t, nchunks> const &rhs) const noexcept;
  bool operator>(Bitset<chunk_t, nchunks> const &rhs) const noexcept;
  bool operator>=(Bitset<chunk_t, nchunks> const &rhs) const noexcept;

  storage_t const &chunks() const noexcept;
  static constexpr size_t nchunkbits() noexcept { return nchunkbits_; }

private:
  // Optimized shift-by-1 for division algorithm
  void shift_left_by_1() noexcept;

  static constexpr size_t nchunkbits_ = std::numeric_limits<chunk_t>::digits;
  static constexpr size_t chunkshift_ = floorlog2(nchunkbits_);
  static constexpr chunk_t chunkmask_ = bitmask<chunk_t>(chunkshift_);
  storage_t chunks_;
};

template <typename chunk_t, int64_t nchunks>
std::string to_string(Bitset<chunk_t, nchunks> const &bits);

template <typename chunk_t, int64_t nchunks>
std::ostream &operator<<(std::ostream &out,
                         Bitset<chunk_t, nchunks> const &bits);

// Conversion functions between Bitset and uint64_t (for testing/interop)
template <typename chunk_t, int64_t nchunks>
Bitset<chunk_t, nchunks> make_bitset(uint64_t value);

template <typename chunk_t, int64_t nchunks>
uint64_t to_uint64(Bitset<chunk_t, nchunks> const &bits);

using BitsetDynamic = Bitset<uint64_t, 0>;
using BitsetStatic1 = Bitset<uint64_t, 1>;
using BitsetStatic2 = Bitset<uint64_t, 2>;
using BitsetStatic4 = Bitset<uint64_t, 4>;
using BitsetStatic8 = Bitset<uint64_t, 8>;

} // namespace xdiag::bits

// Specialization of numeric_limits for Bitset (needed for BitArray)
namespace std {
template <> class numeric_limits<xdiag::bits::Bitset<uint64_t, 0>> {
public:
  static constexpr int digits = numeric_limits<int>::max();
};
template <> class numeric_limits<xdiag::bits::Bitset<uint64_t, 1>> {
public:
  static constexpr int digits = numeric_limits<uint64_t>::digits;
};
template <> class numeric_limits<xdiag::bits::Bitset<uint64_t, 2>> {
public:
  static constexpr int digits = numeric_limits<uint64_t>::digits * 2;
};
template <> class numeric_limits<xdiag::bits::Bitset<uint64_t, 4>> {
public:
  static constexpr int digits = numeric_limits<uint64_t>::digits * 4;
};
template <> class numeric_limits<xdiag::bits::Bitset<uint64_t, 8>> {
public:
  static constexpr int digits = numeric_limits<uint64_t>::digits * 8;
};

} // namespace std
