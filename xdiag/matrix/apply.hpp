// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

template <typename op_t, typename mat_t>
XDIAG_API void apply(op_t const &ops, Block const &block_in,
                     mat_t const &vec_in, Block const &block_out,
                     mat_t &vec_out);

template <typename mat_t>
XDIAG_API void apply(OpSum const &ops, Spinhalf const &block_in,
                     mat_t const &vec_in, Spinhalf const &block_out,
                     mat_t &vec_out);

} // namespace xdiag
