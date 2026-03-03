// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>

#include <xdiag/blocks/spinhalf/spinhalf.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::blocks {

OpSum compile(Spinhalf const &block, OpSum const &ops, std::string mode);

// mode can be one of:
//
// 1: "implementation": compiles to a form understood by the basis
// implementation
//
// 2: "matrix"        : compiles to a "Matrix"-only type, used
// for symmetry detection

} // namespace xdiag::blocks
