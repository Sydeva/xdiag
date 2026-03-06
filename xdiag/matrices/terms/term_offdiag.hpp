// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <xdiag/basis/basis_onthefly.hpp>
#include <xdiag/matrices/fill_functions.hpp>

namespace xdiag::matrices {

template <typename enumeration_t, typename non_zero_term_f,
          typename term_action_f, typename fill_f>
void term_offdiag(basis::BasisOnTheFly<enumeration_t> const &basis_in,
                  basis::BasisOnTheFly<enumeration_t> const &basis_out,
                  non_zero_term_f non_zero_term, term_action_f term_action,
                  fill_f fill) {

  // OpenMP parallel implementation
#ifdef _OPENMP
  int64_t size = basis_in.size();

#pragma omp parallel
  {
    int num_thread = omp_get_thread_num();
    int nthreads = omp_get_num_threads();
    int64_t idx_in = num_thread * (size / nthreads);
    auto begin = basis_in.begin() + idx_in;
    auto end = (num_thread == nthreads - 1)
                   ? basis_in.end()
                   : basis_in.begin() + (num_thread + 1) * (size / nthreads);
    for (auto it = begin; it != end; ++it, ++idx_in) {
      auto spins_in = *it;
      if (non_zero_term(spins_in)) {
        auto [spins_out, coeff] = term_action(spins_in);
        auto idx_out = basis_out.index(spins_out);
        XDIAG_FILL(idx_in, idx_out, coeff);
      }
    }
  }

  // Serial implementation
#else
  int64_t idx_in = 0;
  for (auto const &spins_in : basis_in) {
    if (non_zero_term(spins_in)) {
      auto [spins_out, coeff] = term_action(spins_in);
      auto idx_out = basis_out.index(spins_out);
      fill(idx_in, idx_out, coeff);
    }
    ++idx_in;
  }
#endif
}

} // namespace xdiag::matrices
