// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <utility>
#include <vector>

#include <xdiag/armadillo.hpp>
#include <xdiag/bits/bit_traits.hpp>
#include <xdiag/operators/coeff.hpp>
#include <xdiag/operators/op.hpp>

namespace xdiag::matrices {

template <typename bit_t, typename coeff_t> class NonBranchingOp {
public:
  using cbit_t = typename xdiag::bits::basic_bit_type<bit_t>::type;

  explicit NonBranchingOp(std::vector<int64_t> const &sites,
                          arma::Mat<coeff_t> const &mat,
                          double precision = 1e-12);
  bool isdiagonal() const;
  coeff_t coeff(cbit_t local_state) const;
  std::pair<cbit_t, coeff_t> state_coeff(cbit_t local_state) const;
  bool non_zero_term(cbit_t local_state) const;

  cbit_t extract(bit_t state) const;
  bit_t deposit(cbit_t local_state, bit_t state) const;

private:
  std::vector<int64_t> sites_;
  bit_t mask_;
  bool diagonal_;
  std::vector<bool> non_zero_term_;
  std::vector<cbit_t> state_applied_;
  std::vector<coeff_t> coeff_;
};

template <typename bit_t, typename coeff_t>
std::vector<NonBranchingOp<bit_t, coeff_t>>
non_branching_ops(Coeff const &c, Op const &op, double precision = 1e-12);

} // namespace xdiag::matrices
