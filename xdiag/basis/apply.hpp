// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <memory>
#include <xdiag/basis/basis.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::basis {

template <typename mat_t>
void apply(OpSum const &ops, std::shared_ptr<Basis> const &basis_in,
           mat_t const &mat_in, std::shared_ptr<Basis> const &basis_out,
           mat_t &mat_out);

} // namespace xdiag::basis
