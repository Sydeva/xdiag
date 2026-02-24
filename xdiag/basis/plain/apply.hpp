// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/operators/opsum.hpp>

namespace xdiag::basis::plain {

template <typename basis_t>
void apply(OpSum const &ops, basis_t const &in, basis_t const &out);

} // namespace xdiag::basis::plain
