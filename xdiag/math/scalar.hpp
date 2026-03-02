// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>
#include <variant>
#include <xdiag/math/complex.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

// Scalar is a type-erased numeric value: either double (real) or complex.
// It stays real until complex arithmetic forces a promotion; the variant
// then widens to complex and never narrows back.
//
// Use is<double>() / is<complex>() to query the active type.
// Use as<double>() to extract (throws if complex); as<complex>() always works.
// Use to_real() to narrow back to double if imag() is zero (throws otherwise).
class Scalar {
public:
  using value_t = std::variant<double, complex>;

  XDIAG_API Scalar() = default;
  XDIAG_API Scalar(double value);
  XDIAG_API Scalar(complex value);
  XDIAG_API
  Scalar(int64_t value); // prevents ambiguous integer literal overload
  XDIAG_API Scalar(int value);

  XDIAG_API bool operator==(Scalar const &rhs) const;
  XDIAG_API bool operator!=(Scalar const &rhs) const;

  // Field operations. Mixed real/complex widens the result to complex.
  XDIAG_API Scalar &operator+=(Scalar const &rhs);
  XDIAG_API Scalar &operator-=(Scalar const &rhs);
  XDIAG_API Scalar &operator*=(Scalar const &rhs);
  XDIAG_API Scalar &operator/=(Scalar const &rhs);

  XDIAG_API Scalar operator-() const;
  XDIAG_API Scalar operator+(Scalar const &b) const;
  XDIAG_API Scalar operator-(Scalar const &b) const;
  XDIAG_API Scalar operator*(Scalar const &b) const;
  XDIAG_API Scalar operator/(Scalar const &b) const;

  template <typename T> bool is() const; // is<double>() or is<complex>()
  template <typename T> T as() const;    // as<double>() throws if complex

  XDIAG_API bool isreal() const;
  XDIAG_API double real() const;
  XDIAG_API double imag() const; // 0 for real Scalars
  XDIAG_API double abs() const;
  XDIAG_API Scalar conj() const; // identity for real Scalars
  XDIAG_API Scalar to_real(
      double tol = 1e-12) const; // narrows to double; throws if |imag| > tol
  XDIAG_API bool isapprox(Scalar const &y, double rtol, double atol) const;

private:
  value_t value_;
};

// Returns a zero of the same type (double or complex) as s.
XDIAG_API Scalar zero(Scalar s);

XDIAG_API bool isreal(Scalar const &s);
XDIAG_API double real(Scalar const &s);
XDIAG_API double imag(Scalar const &s);
XDIAG_API double abs(Scalar const &s);
XDIAG_API Scalar conj(Scalar const &s);
XDIAG_API bool isapprox(Scalar const &a, Scalar const &b, double rtol = 1e-12,
                        double atol = 1e-12);

XDIAG_API std::ostream &operator<<(std::ostream &out, Scalar const &v);
XDIAG_API std::string to_string(Scalar const &v);

} // namespace xdiag
