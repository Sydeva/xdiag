// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/algorithms/lanczos/tmatrix.hpp>
#include <xdiag/armadillo.hpp>

namespace xdiag::lanczos {

bool converged_eigenvalues(Tmatrix const &tmat, int n_eigenvalue,
                           double precision);
bool converged_time_evolution(Tmatrix const &tmat, complex tau,
                              double precision, double nrm);

} // namespace xdiag::lanczos
