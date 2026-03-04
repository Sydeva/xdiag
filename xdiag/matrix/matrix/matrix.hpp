// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

template <typename op_t>
XDIAG_API arma::mat matrix(op_t const &op, Block const &block);

template <typename op_t>
XDIAG_API arma::cx_mat matrixC(op_t const &op, Block const &block);

template <typename op_t>
XDIAG_API arma::mat matrix(op_t const &op, Block const &block_in,
                           Block const &block_out);

template <typename op_t>
XDIAG_API arma::cx_mat matrixC(op_t const &op, Block const &block_in,
                               Block const &block_out);

} // namespace xdiag
