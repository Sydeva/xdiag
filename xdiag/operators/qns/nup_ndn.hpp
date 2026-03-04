// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <optional>

#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::operators {

// these functions work for Spinhalf, Electron, and tJ Blocks only
std::optional<int64_t> nup(Op const &op);
std::optional<int64_t> nup(Monomial const &ops);
std::optional<int64_t> nup(OpSum const &ops);
std::optional<int64_t> ndn(Op const &op);
std::optional<int64_t> ndn(Monomial const &ops);
std::optional<int64_t> ndn(OpSum const &ops);

} // namespace xdiag::operators
