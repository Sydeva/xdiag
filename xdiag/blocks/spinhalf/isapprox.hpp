// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <optional>

#include <xdiag/blocks/spinhalf/spinhalf.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::blocks {

bool isapprox(Spinhalf const &block, OpSum const &o1, OpSum const &o2);
std::optional<Scalar> isapprox_multiple(Spinhalf const &block, OpSum const &o1,
                                        OpSum const &o)

} // namespace xdiag::blocks
