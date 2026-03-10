// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/symmetries/action/site_permutation.hpp>

namespace xdiag::symmetries {

// determines norm of a state in an orbit
template <typename bit_t, typename coeff_t>
double norm(bit_t state, SitePermutation const &action,
            arma::Col<coeff_t> const &characters);

} // namespace xdiag::symmetries
