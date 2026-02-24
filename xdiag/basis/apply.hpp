// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <memory>
#include <xdiag/basis/basis.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::basis {

void apply(OpSum const &ops, std::shared_ptr<Basis> const &in,
           std::shared_ptr<Basis> const &out);

} // namespace xdiag::basis
