// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/states/state.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

template <typename op_t> XDIAG_API State apply(op_t const &op, State const &v);
template <typename op_t>
XDIAG_API void apply(op_t const &op, State const &v, State &w);

} // namespace xdiag
