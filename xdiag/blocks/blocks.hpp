// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <cstdint>
#include <ostream>
#include <string>
#include <variant>

#include <xdiag/blocks/spinhalf.hpp>

namespace xdiag {
using Block = std::variant<Spinhalf>;

int64_t dim(Block const &block);
int64_t size(Block const &block);
int64_t nsites(Block const &block);
bool isreal(Block const &block);

std::ostream &operator<<(std::ostream &out, Block const &block);
std::string to_string(Block const &block);

} // namespace xdiag
