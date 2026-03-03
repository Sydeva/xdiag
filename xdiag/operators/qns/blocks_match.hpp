// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag {
bool blocks_match(OpSum const &ops, Block const &b1, Block const &b2);
bool blocks_match(OpSum const &ops, Spinhalf const &b1, Spinhalf const &b2);
} // namespace xdiag
