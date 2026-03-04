// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

namespace xdiag::bits {

template <typename T> struct basic_bit_type {
  using type = typename T::chunk_t;
};
template <> struct basic_bit_type<uint8_t> {
  using type = uint8_t;
};
template <> struct basic_bit_type<uint16_t> {
  using type = uint16_t;
};
template <> struct basic_bit_type<uint32_t> {
  using type = uint32_t;
};
template <> struct basic_bit_type<uint64_t> {
  using type = uint64_t;
};

} // namespace xdiag::bits
