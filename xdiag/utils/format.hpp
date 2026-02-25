// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <complex>

#define FMT_HEADER_ONLY
#include <xdiag/extern/fmt/format.hpp>

namespace fmt {

// Formatter for complex numbers
template <> struct formatter<std::complex<double>> {
  template <typename ParseContext> constexpr auto parse(ParseContext &ctx) {
    return ctx.begin();
  }

  template <typename FormatContext>
  inline auto format(std::complex<double> const &number, FormatContext &ctx) {
    if (std::imag(number) < 0.) {
      return fmt::format_to(ctx.out(), "{0}-i{1}", std::real(number),
                            -std::imag(number));
    } else {
      return fmt::format_to(ctx.out(), "{0}+i{1}", std::real(number),
                            std::imag(number));
    }
  }
};

} // namespace fmt
