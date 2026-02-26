// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "apply.hpp"

#include <xdiag/armadillo.hpp>
#include <xdiag/basis/dispatcher.hpp>
#include <xdiag/basis/plain/apply.hpp>
#include <xdiag/basis/plain/basis_onthefly.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::basis {

template <typename mat_t>
void apply(OpSum const &ops, std::shared_ptr<Basis> const &basis_in,
           mat_t const &mat_in, std::shared_ptr<Basis> const &basis_out,
           mat_t &mat_out) try {
  using namespace combinatorics;

  Dispatcher d;
#define ADD_DISPATCH(BASIS, FUNCTION)                                          \
  d.add<BASIS>([&](BASIS const &in, BASIS const &out) {                        \
    FUNCTION(ops, basis_in, mat_in, basis_out, mat_out);                       \
  });

  ADD_DISPATCH(BasisOnTheFly<Subsets<uint32_t>>, plain::apply);
  ADD_DISPATCH(BasisOnTheFly<Subsets<uint64_t>>, plain::apply);
  ADD_DISPATCH(BasisOnTheFly<Combinations<uint32_t>>, plain::apply);
  ADD_DISPATCH(BasisOnTheFly<Combinations<uint64_t>>, plain::apply);
  ADD_DISPATCH(BasisOnTheFly<LinTable<uint32_t>>, plain::apply);
  ADD_DISPATCH(BasisOnTheFly<LinTable<uint64_t>>, plain::apply);
#undef ADD_DISPATCH

  d.dispatch(basis_in, basis_out);
}
XDIAG_CATCH

template void apply(OpSum const &, std::shared_ptr<Basis> const &,
                    arma::vec const &, std::shared_ptr<Basis> const &,
                    arma::vec &);
template void apply(OpSum const &, std::shared_ptr<Basis> const &,
                    arma::cx_vec const &, std::shared_ptr<Basis> const &,
                    arma::cx_vec &);
template void apply(OpSum const &, std::shared_ptr<Basis> const &,
                    arma::mat const &, std::shared_ptr<Basis> const &,
                    arma::mat &);
template void apply(OpSum const &, std::shared_ptr<Basis> const &,
                    arma::cx_mat const &, std::shared_ptr<Basis> const &,
                    arma::cx_mat &);
} // namespace xdiag::basis
