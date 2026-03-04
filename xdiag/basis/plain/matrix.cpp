// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "matrix.hpp"

#include <xdiag/basis/fill_functions.hpp>
#include <xdiag/basis/plain/basis_onthefly.hpp>
#include <xdiag/basis/plain/implementation/apply_generic.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/math/complex.hpp>

namespace xdiag::basis::plain {

template <typename basis_t, typename coeff_t>
void matrix(OpSum const &ops, basis_t const &basis_in, basis_t const &basis_out,
            coeff_t *mat) try {
  int64_t m = basis_out.size();
  apply_generic<coeff_t>(ops, basis_in, basis_out,
                         [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
                           fill_matrix(mat, m, idx_in, idx_out, val);
                         });
}
XDIAG_CATCH

} // namespace xdiag::basis::plain
using namespace arma;
using namespace xdiag::combinatorics;
using namespace xdiag::basis;
using namespace xdiag::bits;

// Macro for declaring a single apply instantiation
#define INSTANITATE_XDIAG_BASIS_PLAIN_MATRIX(BASIS_TYPE, ENUM_TYPE, INT_TYPE,  \
                                             COEFF_TYPE)                       \
  template void xdiag::basis::plain::matrix(                                   \
      OpSum const &, BASIS_TYPE<ENUM_TYPE<INT_TYPE>> const &,                  \
      BASIS_TYPE<ENUM_TYPE<INT_TYPE>> const &, COEFF_TYPE *);

// BEGIN_INSTANTIATION_GROUP(onthefly_subsets_uint32_t)
INSTANITATE_XDIAG_BASIS_PLAIN_MATRIX(BasisOnTheFly, Subsets, uint32_t, double);
INSTANITATE_XDIAG_BASIS_PLAIN_MATRIX(BasisOnTheFly, Subsets, uint32_t, complex);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_subsets_uint64_t)
INSTANITATE_XDIAG_BASIS_PLAIN_MATRIX(BasisOnTheFly, Subsets, uint64_t, double);
INSTANITATE_XDIAG_BASIS_PLAIN_MATRIX(BasisOnTheFly, Subsets, uint64_t, complex);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_uint32_t)
INSTANITATE_XDIAG_BASIS_PLAIN_MATRIX(BasisOnTheFly, Combinations, uint32_t,
                                     double);
INSTANITATE_XDIAG_BASIS_PLAIN_MATRIX(BasisOnTheFly, Combinations, uint32_t,
                                     complex);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_uint64_t)
INSTANITATE_XDIAG_BASIS_PLAIN_MATRIX(BasisOnTheFly, Combinations, uint64_t,
                                     double);
INSTANITATE_XDIAG_BASIS_PLAIN_MATRIX(BasisOnTheFly, Combinations, uint64_t,
                                     complex);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_bitset_dynamic)
INSTANITATE_XDIAG_BASIS_PLAIN_MATRIX(BasisOnTheFly, Combinations, BitsetDynamic,
                                     double);
INSTANITATE_XDIAG_BASIS_PLAIN_MATRIX(BasisOnTheFly, Combinations, BitsetDynamic,
                                     complex);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_bitset_static_2)
INSTANITATE_XDIAG_BASIS_PLAIN_MATRIX(BasisOnTheFly, Combinations, BitsetStatic2,
                                     double);
INSTANITATE_XDIAG_BASIS_PLAIN_MATRIX(BasisOnTheFly, Combinations, BitsetStatic2,
                                     complex);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_bitset_static_4)
INSTANITATE_XDIAG_BASIS_PLAIN_MATRIX(BasisOnTheFly, Combinations, BitsetStatic4,
                                     double);
INSTANITATE_XDIAG_BASIS_PLAIN_MATRIX(BasisOnTheFly, Combinations, BitsetStatic4,
                                     complex);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_bitset_static_8)
INSTANITATE_XDIAG_BASIS_PLAIN_MATRIX(BasisOnTheFly, Combinations, BitsetStatic8,
                                     double);
INSTANITATE_XDIAG_BASIS_PLAIN_MATRIX(BasisOnTheFly, Combinations, BitsetStatic8,
                                     complex);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_lintable_uint32_t)
INSTANITATE_XDIAG_BASIS_PLAIN_MATRIX(BasisOnTheFly, Combinations, uint32_t,
                                     double);
INSTANITATE_XDIAG_BASIS_PLAIN_MATRIX(BasisOnTheFly, Combinations, uint32_t,
                                     complex);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_lintable_uint64_t)
INSTANITATE_XDIAG_BASIS_PLAIN_MATRIX(BasisOnTheFly, LinTable, uint64_t, double);
INSTANITATE_XDIAG_BASIS_PLAIN_MATRIX(BasisOnTheFly, LinTable, uint64_t,
                                     complex);
// END_INSTANTIATION_GROUP
