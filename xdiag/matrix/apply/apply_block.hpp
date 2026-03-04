// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/operators/opsum.hpp>

namespace xdiag::matrix {

template <typename mat_t, typename block_t>
void apply_block(OpSum const &ops, block_t const &block_in, mat_t const &mat_in,
                 block_t const &block_out, mat_t &mat_out);

} // namespace xdiag::matrix
