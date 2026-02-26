// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/operators/opsum.hpp>

namespace xdiag::basis::plain {

template <typename basis_t, typename coeff_t>
void matrix(OpSum const &ops, basis_t const &basis_in, basis_t const &basis_out,
            coeff_t *mat);

}
