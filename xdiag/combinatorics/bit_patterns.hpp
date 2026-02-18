// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/bits/bitops.hpp>
#include <xdiag/common.hpp>

namespace xdiag::combinatorics {

template <typename bit_t>
bit_t get_next_pattern(bit_t v) noexcept;

template <typename bit_t>
bit_t get_nth_pattern(int64_t n, int64_t nsites, int64_t nupspins);

template <typename bit_t>
int64_t get_n_for_pattern(bit_t pattern, int64_t nsites, int64_t nupspins);

} // namespace xdiag::combinatorics
