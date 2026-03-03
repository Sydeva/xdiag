// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "apply.hpp"

#include <cstdint>

#include <xdiag/armadillo.hpp>
#include <xdiag/basis/fill_functions.hpp>
#include <xdiag/basis/plain/basis_onthefly.hpp>
#include <xdiag/basis/plain/implementation/apply_generic.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::basis::plain {

template <typename mat_t, typename basis_t>
void apply(OpSum const &ops, basis_t const &basis_in, mat_t const &mat_in,
           basis_t const &basis_out, mat_t &mat_out) try {
  using coeff_t = typename mat_t::elem_type;
  apply_generic<coeff_t>(ops, basis_in, basis_out,
                         [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
                           fill_apply(mat_in, mat_out, idx_in, idx_out, val);
                         });
}
XDIAG_CATCH

} // namespace xdiag::basis::plain

using namespace arma;
using namespace xdiag::combinatorics;
using namespace xdiag::basis;
using namespace xdiag::bits;

// Macro for declaring a single apply instantiation
#define INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BASIS_TYPE, ENUM_TYPE, INT_TYPE,   \
                                            MAT_TYPE)                          \
  template void xdiag::basis::plain::apply(                                    \
      OpSum const &, BASIS_TYPE<ENUM_TYPE<INT_TYPE>> const &,                  \
      MAT_TYPE const &, BASIS_TYPE<ENUM_TYPE<INT_TYPE>> const &, MAT_TYPE &);

// BEGIN_INSTANTIATION_GROUP(onthefly_subsets_uint32_t)
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Subsets, uint32_t, vec);
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Subsets, uint32_t, cx_vec);
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Subsets, uint32_t, mat);
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Subsets, uint32_t, cx_mat);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_subsets_uint64_t)
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Subsets, uint64_t, vec);
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Subsets, uint64_t, cx_vec);
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Subsets, uint64_t, mat);
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Subsets, uint64_t, cx_mat);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_uint32_t)
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Combinations, uint32_t, vec);
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Combinations, uint32_t,
                                    cx_vec);
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Combinations, uint32_t, mat);
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Combinations, uint32_t,
                                    cx_mat);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_uint64_t)
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Combinations, uint64_t, vec);
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Combinations, uint64_t,
                                    cx_vec);
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Combinations, uint64_t, mat);
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Combinations, uint64_t,
                                    cx_mat);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_bitset_dynamic)
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Combinations, BitsetDynamic,
                                    vec);
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Combinations, BitsetDynamic,
                                    cx_vec);
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Combinations, BitsetDynamic,
                                    mat);
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Combinations, BitsetDynamic,
                                    cx_mat);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_bitset_static_2)
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Combinations, BitsetStatic2,
                                    vec);
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Combinations, BitsetStatic2,
                                    cx_vec);
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Combinations, BitsetStatic2,
                                    mat);
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Combinations, BitsetStatic2,
                                    cx_mat);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_bitset_static_4)
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Combinations, BitsetStatic4,
                                    vec);
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Combinations, BitsetStatic4,
                                    cx_vec);
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Combinations, BitsetStatic4,
                                    mat);
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Combinations, BitsetStatic4,
                                    cx_mat);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_bitset_static_8)
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Combinations, BitsetStatic8,
                                    vec);
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Combinations, BitsetStatic8,
                                    cx_vec);
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Combinations, BitsetStatic8,
                                    mat);
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, Combinations, BitsetStatic8,
                                    cx_mat);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_lintable_uint32_t)
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, LinTable, uint32_t, vec);
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, LinTable, uint32_t, cx_vec);
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, LinTable, uint32_t, mat);
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, LinTable, uint32_t, cx_mat);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_lintable_uint64_t)
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, LinTable, uint64_t, vec);
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, LinTable, uint64_t, cx_vec);
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, LinTable, uint64_t, mat);
INSTANITATE_XDIAG_BASIS_PLAIN_APPLY(BasisOnTheFly, LinTable, uint64_t, cx_mat);
// END_INSTANTIATION_GROUP

#undef INSTANITATE_XDIAG_BASIS_PLAIN_APPLY
