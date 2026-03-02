// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <memory>
#include <xdiag/basis/basis.hpp>

namespace xdiag::basis {

template <typename op_t, typename mat_t>
void apply(op_t const &ops, std::shared_ptr<Basis> const &basis_in,
           mat_t const &mat_in, std::shared_ptr<Basis> const &basis_out,
           mat_t &mat_out);

} // namespace xdiag::basis
