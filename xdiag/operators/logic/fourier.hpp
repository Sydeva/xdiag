// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <vector>

#include <xdiag/common.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/symmetries/representation.hpp>

namespace xdiag {

XDIAG_API OpSum fourier(Op const &op,
                        std::vector<Representation> const &irreps);

XDIAG_API OpSum fourier(OpSum const &ops,
                        std::vector<Representation> const &irreps);

} // namespace xdiag

