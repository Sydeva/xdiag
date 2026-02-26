// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <xdiag/basis/plain/basis_onthefly.hpp>

namespace xdiag::basis {

template <typename enumeration_t, typename term_coeff_f, typename fill_f>
void apply_diag(BasisOnTheFly<enumeration_t> const &basis,
                term_coeff_f term_coeff, fill_f fill) {

  // OpenMP parallel implementation
#ifdef _OPENMP
  int64_t size = basis.size();

#pragma omp parallel
  {
    int num_thread = omp_get_thread_num();
    int nthreads = omp_get_num_threads();
    int64_t idx = num_thread * (size / nthreads);
    auto begin = basis.begin() + idx;
    auto end = (num_thread == nthreads - 1)
                   ? basis.end()
                   : basis.begin() + (num_thread + 1) * (size / nthreads);

    for (auto it = begin; it != end; ++it, ++idx) {
      auto coeff = term_coeff(*it);
      fill(idx, idx, coeff);
    }
  }

  // Serial implementation
#else
  int64_t idx = 0;
  for (auto spins : basis) {
    auto coeff = term_coeff(spins);
    fill(idx, idx, coeff);
    ++idx;
  }
#endif
}

} // namespace xdiag::basis
