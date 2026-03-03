// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

template <typename op_t, typename mat_t, typename block_t>
XDIAG_API void apply_block(op_t const &ops, block_t const &block_in,
                           mat_t const &mat_in, block_t const &block_out,
                           mat_t &mat_out);

} // namespace xdiag
