// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/armadillo.hpp>

namespace xdiag::math {

// Promote real armadillo vectors/matrices to complex equivalents.
// If the argument is already complex, it is returned unchanged (identity).
arma::cx_vec to_cx_vec(arma::vec const &A);
arma::cx_vec to_cx_vec(arma::cx_vec const &A);
arma::cx_mat to_cx_mat(arma::mat const &A);
arma::cx_mat to_cx_mat(arma::cx_mat const &A);

} // namespace xdiag::math
