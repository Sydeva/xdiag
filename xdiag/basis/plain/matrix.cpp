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

using namespace combinatorics;
template void matrix(OpSum const &, BasisOnTheFly<Subsets<uint32_t>> const &,
                     BasisOnTheFly<Subsets<uint32_t>> const &, double *);
template void matrix(OpSum const &, BasisOnTheFly<Subsets<uint32_t>> const &,
                     BasisOnTheFly<Subsets<uint32_t>> const &, complex *);
template void matrix(OpSum const &, BasisOnTheFly<Subsets<uint64_t>> const &,
                     BasisOnTheFly<Subsets<uint64_t>> const &, double *);
template void matrix(OpSum const &, BasisOnTheFly<Subsets<uint64_t>> const &,
                     BasisOnTheFly<Subsets<uint64_t>> const &, complex *);

template void matrix(OpSum const &,
                     BasisOnTheFly<Combinations<uint32_t>> const &,
                     BasisOnTheFly<Combinations<uint32_t>> const &, double *);
template void matrix(OpSum const &,
                     BasisOnTheFly<Combinations<uint32_t>> const &,
                     BasisOnTheFly<Combinations<uint32_t>> const &, complex *);
template void matrix(OpSum const &,
                     BasisOnTheFly<Combinations<uint64_t>> const &,
                     BasisOnTheFly<Combinations<uint64_t>> const &, double *);
template void matrix(OpSum const &,
                     BasisOnTheFly<Combinations<uint64_t>> const &,
                     BasisOnTheFly<Combinations<uint64_t>> const &, complex *);

template void matrix(OpSum const &, BasisOnTheFly<LinTable<uint32_t>> const &,
                     BasisOnTheFly<LinTable<uint32_t>> const &, double *);
template void matrix(OpSum const &, BasisOnTheFly<LinTable<uint32_t>> const &,
                     BasisOnTheFly<LinTable<uint32_t>> const &, complex *);
template void matrix(OpSum const &, BasisOnTheFly<LinTable<uint64_t>> const &,
                     BasisOnTheFly<LinTable<uint64_t>> const &, double *);
template void matrix(OpSum const &, BasisOnTheFly<LinTable<uint64_t>> const &,
                     BasisOnTheFly<LinTable<uint64_t>> const &, complex *);

} // namespace xdiag::basis::plain
