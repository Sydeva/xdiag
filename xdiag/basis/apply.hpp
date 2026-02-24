// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/basis/basis.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::basis {

void apply(OpSum const &ops, Basis *in, Basis *out);

} // namespace xdiag::basis
