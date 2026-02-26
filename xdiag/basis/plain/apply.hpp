// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/operators/opsum.hpp>

namespace xdiag::basis::plain {

template <typename mat_t, typename basis_t>
void apply(OpSum const &ops, basis_t const &basis_in, mat_t const &mat_in,
           basis_t const &basis_out, mat_t &mat_out);
 
} // namespace xdiag::basis::plain
