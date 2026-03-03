// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>
#include <variant>

#include <xdiag/armadillo.hpp>
#include <xdiag/math/scalar.hpp>
#include <xdiag/math/vector.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

// Matrix is a type-erased dense matrix: either arma::mat (real) or
// arma::cx_mat (complex). Mixed-type arithmetic widens to complex.
//
// Use is<arma::mat>() / is<arma::cx_mat>() to query; as<T>() to extract
// (as<arma::cx_mat>() always succeeds; as<arma::mat>() throws if complex).
// hc() returns the conjugate transpose (Hermitian conjugate).
// to_real() narrows back to arma::mat if the imaginary part is zero.
class Matrix {
public:
  using value_t = std::variant<arma::mat, arma::cx_mat>;

  XDIAG_API Matrix() = default;
  XDIAG_API Matrix(arma::mat const &mat);
  XDIAG_API Matrix(arma::cx_mat const &mat);

  XDIAG_API bool operator==(Matrix const &rhs) const;
  XDIAG_API bool operator!=(Matrix const &rhs) const;
  XDIAG_API bool operator<(Matrix const &rhs) const;

  // Linear combination (scalar multiply/add). Mixed real/complex widens to
  // complex.
  XDIAG_API Matrix &operator+=(Matrix const &rhs);
  XDIAG_API Matrix &operator-=(Matrix const &rhs);
  XDIAG_API Matrix &operator*=(Scalar const &rhs);
  XDIAG_API Matrix &operator/=(Scalar const &rhs);

  XDIAG_API Matrix operator-() const;
  XDIAG_API Matrix operator+(Matrix const &b) const;
  XDIAG_API Matrix operator-(Matrix const &b) const;
  XDIAG_API Matrix operator*(Scalar const &rhs) const;
  XDIAG_API Matrix operator/(Scalar const &rhs) const;

  // Matrix-matrix and matrix-vector products. Widens to complex if either side
  // is complex.
  XDIAG_API Matrix operator*(Matrix const &rhs) const;
  XDIAG_API Vector operator*(Vector const &rhs) const;

  template <typename T> XDIAG_API bool is() const;
  template <typename T> XDIAG_API T as() const;

  XDIAG_API int64_t n_rows() const;
  XDIAG_API int64_t n_cols() const;
  XDIAG_API bool isreal() const;
  XDIAG_API arma::mat real() const;
  XDIAG_API arma::mat imag() const; // zero matrix for real Matrices
  XDIAG_API Matrix hc() const;      // conjugate transpose
  XDIAG_API Matrix to_real(double tol = 1e-12)
      const; // narrows to arma::mat; throws if any |imag| > tol
  XDIAG_API bool isapprox(Matrix const &y, double rtol = 1e-12,
                          double atol = 1e-12) const;

private:
  value_t mat_;
};

XDIAG_API bool isreal(Matrix const &m);
XDIAG_API arma::mat real(Matrix const &m);
XDIAG_API arma::mat imag(Matrix const &m);
XDIAG_API Matrix hc(Matrix const &m);
XDIAG_API bool isapprox(Matrix const &a, Matrix const &b, double rtol = 1e-12,
                        double atol = 1e-12);

XDIAG_API std::ostream &operator<<(std::ostream &out, Matrix const &mat);
XDIAG_API std::string to_string(Matrix const &mat);

} // namespace xdiag
