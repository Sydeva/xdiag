// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/math/complex.hpp>
#include <xdiag/blocks/blocks.hpp>

namespace xdiag::math {

// Internal routines (used to differentiate dot product for distributed blocks)
double dot(Block const &block, arma::vec const &v, arma::vec const &w);
complex dot(Block const &block, arma::cx_vec const &v, arma::cx_vec const &w);

} // namespace xdiag::math
