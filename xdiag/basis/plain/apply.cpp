// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "apply.hpp"

#include <cstdint>

#include <xdiag/armadillo.hpp>
#include <xdiag/basis/fill_functions.hpp>
#include <xdiag/basis/plain/basis_onthefly.hpp>
#include <xdiag/basis/plain/implementation/apply_generic.hpp>
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


// Macro for declaring a single apply instantiation
#define INSTANTIATE_APPLY(BASIS_TYPE, ENUM_TYPE, INT_TYPE, MAT_TYPE)           \
  template void xdiag::basis::plain::apply(                                    \
      OpSum const &, BASIS_TYPE<ENUM_TYPE<INT_TYPE>> const &,                  \
      MAT_TYPE const &, BASIS_TYPE<ENUM_TYPE<INT_TYPE>> const &, MAT_TYPE &);

// BEGIN_INSTANTIATION_GROUP(onthefly_subsets_uint32_t)
INSTANTIATE_APPLY(BasisOnTheFly, Subsets, uint32_t, vec);
INSTANTIATE_APPLY(BasisOnTheFly, Subsets, uint32_t, cx_vec);
INSTANTIATE_APPLY(BasisOnTheFly, Subsets, uint32_t, mat);
INSTANTIATE_APPLY(BasisOnTheFly, Subsets, uint32_t, cx_mat);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_subsets_uint64_t)
INSTANTIATE_APPLY(BasisOnTheFly, Subsets, uint64_t, vec);
INSTANTIATE_APPLY(BasisOnTheFly, Subsets, uint64_t, cx_vec);
INSTANTIATE_APPLY(BasisOnTheFly, Subsets, uint64_t, mat);
INSTANTIATE_APPLY(BasisOnTheFly, Subsets, uint64_t, cx_mat);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_uint32_t)
INSTANTIATE_APPLY(BasisOnTheFly, Combinations, uint32_t, vec);
INSTANTIATE_APPLY(BasisOnTheFly, Combinations, uint32_t, cx_vec);
INSTANTIATE_APPLY(BasisOnTheFly, Combinations, uint32_t, mat);
INSTANTIATE_APPLY(BasisOnTheFly, Combinations, uint32_t, cx_mat);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_combinations_uint64_t)
INSTANTIATE_APPLY(BasisOnTheFly, Combinations, uint64_t, vec);
INSTANTIATE_APPLY(BasisOnTheFly, Combinations, uint64_t, cx_vec);
INSTANTIATE_APPLY(BasisOnTheFly, Combinations, uint64_t, mat);
INSTANTIATE_APPLY(BasisOnTheFly, Combinations, uint64_t, cx_mat);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_lintable_uint32_t)
INSTANTIATE_APPLY(BasisOnTheFly, LinTable, uint32_t, vec);
INSTANTIATE_APPLY(BasisOnTheFly, LinTable, uint32_t, cx_vec);
INSTANTIATE_APPLY(BasisOnTheFly, LinTable, uint32_t, mat);
INSTANTIATE_APPLY(BasisOnTheFly, LinTable, uint32_t, cx_mat);
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(onthefly_lintable_uint64_t)
INSTANTIATE_APPLY(BasisOnTheFly, LinTable, uint64_t, vec);
INSTANTIATE_APPLY(BasisOnTheFly, LinTable, uint64_t, cx_vec);
INSTANTIATE_APPLY(BasisOnTheFly, LinTable, uint64_t, mat);
INSTANTIATE_APPLY(BasisOnTheFly, LinTable, uint64_t, cx_mat);
// END_INSTANTIATION_GROUP

#undef INSTANTIATE_APPLY
