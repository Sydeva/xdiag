// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>
#include <variant>

#include <xdiag/math/complex.hpp>
#include <xdiag/math/scalar.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

// Coeff represents an operator coefficient that is either a numeric Scalar
// (double or complex) or a named string constant to be resolved later via
// OpSum::plain().
class Coeff {
public:
  using value_t = std::variant<Scalar, std::string>;

  XDIAG_API Coeff() = default;
  XDIAG_API explicit Coeff(std::string value);
  XDIAG_API explicit Coeff(const char *value);
  XDIAG_API explicit Coeff(double value);
  XDIAG_API explicit Coeff(complex value);
  XDIAG_API explicit Coeff(Scalar value);

  XDIAG_API bool operator==(Coeff const &rhs) const;
  XDIAG_API bool operator!=(Coeff const &rhs) const;

  XDIAG_API bool isscalar() const;
  XDIAG_API bool isstring() const;
  XDIAG_API Scalar scalar() const;
  XDIAG_API std::string string() const;

private:
  value_t value_;
};

XDIAG_API bool isscalar(Coeff const &c);
XDIAG_API bool isstring(Coeff const &c);
XDIAG_API Scalar scalar(Coeff const &c);
XDIAG_API std::string string(Coeff const &c);

// Multiply two Coeffs. Both must be scalar (throws if either is a string).
// Call OpSum::plain() first to resolve named coefficients before multiplying.
XDIAG_API Coeff operator*(Coeff const &lhs, Coeff const &rhs);
XDIAG_API Coeff operator*(Coeff const &lhs, Scalar const &rhs);
XDIAG_API Coeff operator*(Scalar const &lhs, Coeff const &rhs);

XDIAG_API std::ostream &operator<<(std::ostream &out, Coeff const &c);
XDIAG_API std::string to_string(Coeff const &c);

} // namespace xdiag
