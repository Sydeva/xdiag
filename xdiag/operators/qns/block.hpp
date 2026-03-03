// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

XDIAG_API Block block(OpSum const &ops, Block const &block);
XDIAG_API Spinhalf block(OpSum const &ops, Spinhalf const &block);

} // namespace xdiag
