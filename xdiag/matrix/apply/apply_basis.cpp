// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "apply.hpp"

#include <xdiag/armadillo.hpp>
#include <xdiag/basis/plain/apply.hpp>
#include <xdiag/basis/plain/basis_onthefly.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/matrix/implementation/dispatcher.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::matrix {

template <typename mat_t>
void apply_basis(OpSum const &ops, std::shared_ptr<Basis> const &basis_in,
                 mat_t const &mat_in, std::shared_ptr<Basis> const &basis_out,
                 mat_t &mat_out) try {
  using namespace combinatorics;
  using namespace bits;

  Dispatcher d;
#define ADD_DISPATCH(BASIS, FUNCTION)                                          \
  d.add<BASIS>([&](BASIS const &in, BASIS const &out) {                        \
    FUNCTION(ops, in, mat_in, out, mat_out);                                   \
  });

  ADD_DISPATCH(BasisOnTheFly<Subsets<uint32_t>>, plain::apply);
  ADD_DISPATCH(BasisOnTheFly<Subsets<uint64_t>>, plain::apply);
  ADD_DISPATCH(BasisOnTheFly<Combinations<uint32_t>>, plain::apply);
  ADD_DISPATCH(BasisOnTheFly<Combinations<uint64_t>>, plain::apply);
  ADD_DISPATCH(BasisOnTheFly<Combinations<BitsetDynamic>>, plain::apply);
  ADD_DISPATCH(BasisOnTheFly<Combinations<BitsetStatic2>>, plain::apply);
  ADD_DISPATCH(BasisOnTheFly<Combinations<BitsetStatic4>>, plain::apply);
  ADD_DISPATCH(BasisOnTheFly<Combinations<BitsetStatic8>>, plain::apply);
  ADD_DISPATCH(BasisOnTheFly<LinTable<uint32_t>>, plain::apply);
  ADD_DISPATCH(BasisOnTheFly<LinTable<uint64_t>>, plain::apply);
#undef ADD_DISPATCH

  d.dispatch(basis_in, basis_out);
}
XDIAG_CATCH

#define INSTANTIATE_XDIAG_BASIS_APPLY(MAT_TYPE)                                \
  template void apply(OpSum const &, std::shared_ptr<Basis> const &,           \
                      MAT_TYPE const &, std::shared_ptr<Basis> const &,        \
                      MAT_TYPE &);

INSTANTIATE_XDIAG_BASIS_APPLY(arma::vec);
INSTANTIATE_XDIAG_BASIS_APPLY(arma::cx_vec);
INSTANTIATE_XDIAG_BASIS_APPLY(arma::mat);
INSTANTIATE_XDIAG_BASIS_APPLY(arma::cx_mat);

} // namespace xdiag::matrix
