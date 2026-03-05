// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "apply.hpp"

#include <xdiag/matrix/utils/fill_functions.hpp>

namespace xdiag::matrix::spinhalf {

template <typename basis_t, typename mat_t>
void apply(OpSum const &ops, basis_t const &basis_in, mat_t const &mat_in,
           basis_t const &basis_out, mat_t &mat_out) try {
  using coeff_t = typename mat_t::elem_type;
  spinhalf::matrix_generic(ops, basis_in, basis_out,
                           [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
                             fill_apply(mat_in, mat_out, idx_in, idx_out, val);
                           });
}
XDIAG_CATCH

} // namespace xdiag::matrix::spinhalf

#define INSTANTIATE_XDIAG_MATRIX_SPINHALF_APPLY(BASIS_TYPE, ENUM_TYPE,         \
                                                INT_TYPE, MAT_TYPE)            \
  template void xdiag::matrix::spinhalf::apply(                                \
      OpSum const &, BASIS_TYPE<ENUM_TYPE<INT_TYPE>> const &,                  \
      MAT_TYPE const &, BASIS_TYPE<ENUM_TYPE<INT_TYPE>> const &, MAT_TYPE &);

#define INSTANTIATE_XDIAG_MATRIX_SPINHALF_APPLY_MATRICES(BASIS_TYPE,           \
                                                         ENUM_TYPE, INT_TYPE)  \
  INSTANTIATE_XDIAG_MATRIX_SPINHALF_APPLY(BASIS_TYPE, ENUM_TYPE, INT_TYPE,     \
                                          vec);                                \
  INSTANTIATE_XDIAG_MATRIX_SPINHALF_APPLY(BASIS_TYPE, ENUM_TYPE, INT_TYPE,     \
                                          cx_vec);                             \
  INSTANTIATE_XDIAG_MATRIX_SPINHALF_APPLY(BASIS_TYPE, ENUM_TYPE, INT_TYPE,     \
                                          mat);                                \
  INSTANTIATE_XDIAG_MATRIX_SPINHALF_APPLY(BASIS_TYPE, ENUM_TYPE, INT_TYPE,     \
                                          cx_mat);

// BEGIN_INSTANTIATION_GROUP(onthefly_subsets_uint32_t)
INSTANTIATE_XDIAG_BASIS_PLAIN_APPLY_MATRICES(BasisOnTheFly, Subsets, uint32_t);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_subsets_uint64_t)
INSTANTIATE_XDIAG_BASIS_PLAIN_APPLY_MATRICES(BasisOnTheFly, Subsets, uint64_t);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_uint32_t)
INSTANTIATE_XDIAG_BASIS_PLAIN_APPLY_MATRICES(BasisOnTheFly, Combinations,
                                             uint32_t);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_uint64_t)
INSTANTIATE_XDIAG_BASIS_PLAIN_APPLY_MATRICES(BasisOnTheFly, Combinations,
                                             uint64_t);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_lintable_uint32_t)
INSTANTIATE_XDIAG_BASIS_PLAIN_APPLY_MATRICES(BasisOnTheFly, LinTable, uint32_t);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_lintable_uint64_t)
INSTANTIATE_XDIAG_BASIS_PLAIN_APPLY_MATRICES(BasisOnTheFly, LinTable, uint64_t);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_bitset_dynamic)
INSTANTIATE_XDIAG_BASIS_PLAIN_APPLY_MATRICES(BasisOnTheFly, Combinations,
                                             BitsetDynamic);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_bitset_static_2)
INSTANTIATE_XDIAG_BASIS_PLAIN_APPLY_MATRICES(BasisOnTheFly, Combinations,
                                             BitsetStatic2);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_bitset_static_4)
INSTANTIATE_XDIAG_BASIS_PLAIN_APPLY_MATRICES(BasisOnTheFly, Combinations,
                                             BitsetStatic4);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_bitset_static_8)
INSTANTIATE_XDIAG_BASIS_PLAIN_APPLY_MATRICES(BasisOnTheFly, Combinations,
                                             BitsetStatic8);
// END_INSTANTIATION_GROUP
