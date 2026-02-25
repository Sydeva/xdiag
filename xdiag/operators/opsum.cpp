// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "opsum.hpp"

#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>
#include <xdiag/utils/to_string_generic.hpp>

namespace xdiag {

// --- Constructors ---

OpSum::OpSum(Op const &op) : terms_({{Coeff(1.0), Monomial(op)}}) {}
OpSum::OpSum(Monomial const &mono) : terms_({{Coeff(1.0), mono}}) {}
OpSum::OpSum(Coeff const &coeff, Op const &op)
    : terms_({{coeff, Monomial(op)}}) {}
OpSum::OpSum(Coeff const &coeff, Monomial const &mono)
    : terms_({{coeff, mono}}) {}
OpSum::OpSum(std::string const &coeff, Op const &op)
    : OpSum(Coeff(coeff), op) {}
OpSum::OpSum(double coeff, Op const &op) : OpSum(Coeff(coeff), op) {}
OpSum::OpSum(complex coeff, Op const &op) : OpSum(Coeff(coeff), op) {}
OpSum::OpSum(std::string const &coeff, Monomial const &mono)
    : OpSum(Coeff(coeff), mono) {}
OpSum::OpSum(double coeff, Monomial const &mono) : OpSum(Coeff(coeff), mono) {}
OpSum::OpSum(complex coeff, Monomial const &mono) : OpSum(Coeff(coeff), mono) {}

// --- Scalar scaling ---

OpSum &OpSum::operator*=(double scalar) try {
  return operator*=(Scalar(scalar));
}
XDIAG_CATCH

OpSum &OpSum::operator*=(complex scalar) try {
  return operator*=(Scalar(scalar));
}
XDIAG_CATCH

OpSum &OpSum::operator*=(Scalar const &scalar) try {
  for (auto &[c, m] : terms_) {
    if (c.isscalar()) {
      c = Coeff(c.scalar() * scalar);
    } else {
      auto it = params_.find(c.string());
      if (it != params_.end()) {
        c = Coeff(it->second * scalar);
      } else {
        XDIAG_THROW(fmt::format(
            "Cannot scale OpSum: coefficient \"{}\" has not been defined. "
            "Call plain() first to resolve named coefficients.",
            c.string()));
      }
    }
  }
  return *this;
}
XDIAG_CATCH

OpSum &OpSum::operator/=(double scalar) try {
  return operator/=(Scalar(scalar));
}
XDIAG_CATCH

OpSum &OpSum::operator/=(complex scalar) try {
  return operator/=(Scalar(scalar));
}
XDIAG_CATCH

OpSum &OpSum::operator/=(Scalar const &scalar) try {
  return operator*=(Scalar(1.0) / scalar);
}
XDIAG_CATCH

// --- Addition / Subtraction ---

void OpSum::merge_params(std::map<std::string, Scalar> const &other) try {
  for (auto const &[key, val] : other) {
    auto it = params_.find(key);
    if (it == params_.end()) {
      params_[key] = val;
    } else if (it->second != val) {
      XDIAG_THROW(fmt::format(
          "Conflicting values for named coefficient \"{}\": {} vs {}", key,
          to_string(it->second), to_string(val)));
    }
  }
}
XDIAG_CATCH

OpSum &OpSum::operator+=(OpSum const &ops) try {
  terms_.insert(terms_.end(), ops.terms_.begin(), ops.terms_.end());
  merge_params(ops.params_);
  return *this;
}
XDIAG_CATCH

OpSum &OpSum::operator+=(Op const &op) try { return operator+=(OpSum(op)); }
XDIAG_CATCH

OpSum OpSum::operator+(OpSum const &ops) const try {
  OpSum result = *this;
  result += ops;
  return result;
}
XDIAG_CATCH

OpSum OpSum::operator+(Op const &op) const try { return operator+(OpSum(op)); }
XDIAG_CATCH

OpSum &OpSum::operator-=(OpSum const &ops) try {
  OpSum neg = ops;
  neg *= Scalar(-1.0);
  return operator+=(neg);
}
XDIAG_CATCH

OpSum &OpSum::operator-=(Op const &op) try { return operator-=(OpSum(op)); }
XDIAG_CATCH

OpSum OpSum::operator-(OpSum const &ops) const try {
  OpSum result = *this;
  result -= ops;
  return result;
}
XDIAG_CATCH

OpSum OpSum::operator-(Op const &op) const try { return operator-(OpSum(op)); }
XDIAG_CATCH

OpSum OpSum::operator-() const try {
  OpSum result = *this;
  result *= Scalar(-1.0);
  return result;
}
XDIAG_CATCH

// --- Algebra product ---

OpSum OpSum::operator*(OpSum const &rhs) const try {
  OpSum result;
  for (auto const &[cl, ml] : terms_) {
    for (auto const &[cr, mr] : rhs.terms_) {
      result.terms_.push_back({cl * cr, ml * mr});
    }
  }
  result.merge_params(params_);
  result.merge_params(rhs.params_);
  return result;
}
XDIAG_CATCH

OpSum &OpSum::operator*=(OpSum const &rhs) try {
  *this = operator*(rhs);
  return *this;
}
XDIAG_CATCH

OpSum OpSum::operator*(Op const &rhs) const try {
  return operator*(OpSum(rhs));
}
XDIAG_CATCH

OpSum &OpSum::operator*=(Op const &rhs) try { return operator*=(OpSum(rhs)); }
XDIAG_CATCH

OpSum OpSum::operator*(Monomial const &rhs) const try {
  return operator*(OpSum(rhs));
}
XDIAG_CATCH

OpSum &OpSum::operator*=(Monomial const &rhs) try {
  return operator*=(OpSum(rhs));
}
XDIAG_CATCH

// --- Named parameters ---

Scalar &OpSum::operator[](std::string const &name) { return params_[name]; }
Scalar const &OpSum::operator[](std::string const &name) const try {
  return params_.at(name);
}
XDIAG_CATCH

OpSum OpSum::plain() const try {
  OpSum result;
  for (auto const &[coeff, mono] : terms_) {
    if (coeff.isscalar()) {
      result.terms_.push_back({coeff, mono});
    } else {
      auto it = params_.find(coeff.string());
      if (it != params_.end()) {
        result.terms_.push_back({Coeff(it->second), mono});
      } else {
        XDIAG_THROW(fmt::format(
            "Cannot make OpSum plain: coefficient \"{}\" has not been defined.",
            coeff.string()));
      }
    }
  }
  return result;
}
XDIAG_CATCH

// --- Access ---

std::vector<Term> const &OpSum::terms() const noexcept { return terms_; }
std::map<std::string, Scalar> const &OpSum::params() const noexcept {
  return params_;
}
int64_t OpSum::size() const noexcept { return (int64_t)terms_.size(); }
OpSum::iterator_t OpSum::begin() const noexcept { return terms_.begin(); }
OpSum::iterator_t OpSum::end() const noexcept { return terms_.end(); }

bool OpSum::operator==(OpSum const &rhs) const {
  return (terms_ == rhs.terms_) && (params_ == rhs.params_);
}
bool OpSum::operator!=(OpSum const &rhs) const { return !operator==(rhs); }

// --- Free operators ---

// Coeff/scalar * Op
OpSum operator*(double coeff, Op const &op) { return OpSum(Coeff(coeff), op); }
OpSum operator*(complex coeff, Op const &op) { return OpSum(Coeff(coeff), op); }
OpSum operator*(std::string const &coeff, Op const &op) {
  return OpSum(Coeff(coeff), op);
}
OpSum operator*(Coeff const &coeff, Op const &op) { return OpSum(coeff, op); }
OpSum operator*(Op const &op, double coeff) { return coeff * op; }
OpSum operator*(Op const &op, complex coeff) { return coeff * op; }
OpSum operator*(Op const &op, std::string const &coeff) { return coeff * op; }
OpSum operator*(Op const &op, Coeff const &coeff) { return coeff * op; }

// Coeff/scalar * Monomial
OpSum operator*(double coeff, Monomial const &mono) {
  return OpSum(Coeff(coeff), mono);
}
OpSum operator*(complex coeff, Monomial const &mono) {
  return OpSum(Coeff(coeff), mono);
}
OpSum operator*(std::string const &coeff, Monomial const &mono) {
  return OpSum(Coeff(coeff), mono);
}
OpSum operator*(Coeff const &coeff, Monomial const &mono) {
  return OpSum(coeff, mono);
}
OpSum operator*(Monomial const &mono, double coeff) { return coeff * mono; }
OpSum operator*(Monomial const &mono, complex coeff) { return coeff * mono; }
OpSum operator*(Monomial const &mono, std::string const &coeff) {
  return coeff * mono;
}
OpSum operator*(Monomial const &mono, Coeff const &coeff) {
  return coeff * mono;
}

// scalar * OpSum (scalar scaling)
OpSum operator*(double scalar, OpSum const &ops) {
  return Scalar(scalar) * ops;
}
OpSum operator*(complex scalar, OpSum const &ops) {
  return Scalar(scalar) * ops;
}
OpSum operator*(Scalar const &scalar, OpSum const &ops) {
  OpSum result = ops;
  result *= scalar;
  return result;
}
OpSum operator*(OpSum const &ops, double scalar) { return scalar * ops; }
OpSum operator*(OpSum const &ops, complex scalar) { return scalar * ops; }
OpSum operator*(OpSum const &ops, Scalar const &scalar) { return scalar * ops; }

OpSum operator/(OpSum const &ops, double scalar) {
  return ops / Scalar(scalar);
}
OpSum operator/(OpSum const &ops, complex scalar) {
  return ops / Scalar(scalar);
}
OpSum operator/(OpSum const &ops, Scalar const &scalar) {
  return ops * (Scalar(1.0) / scalar);
}

// Algebra products: Op * OpSum, Monomial * OpSum
OpSum operator*(Op const &lhs, OpSum const &rhs) { return OpSum(lhs) * rhs; }
OpSum operator*(Monomial const &lhs, OpSum const &rhs) {
  return OpSum(lhs) * rhs;
}

// --- I/O ---

std::ostream &operator<<(std::ostream &out, OpSum const &ops) {
  out << "Terms:\n";
  out << "------\n";
  for (auto const &[coeff, mono] : ops) {
    out << coeff << " * " << mono << "\n";
  }
  if (!ops.params().empty()) {
    out << "\nParameters:\n";
    out << "-----------\n";
    for (auto const &[name, val] : ops.params()) {
      out << name << ": " << val << "\n";
    }
  }
  return out;
}

std::string to_string(OpSum const &ops) { return to_string_generic(ops); }

} // namespace xdiag
