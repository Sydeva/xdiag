// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <xdiag/operators/opsum.hpp>
#include <xdiag/blocks/spinhalf.hpp>

namespace xdiag::blocks {

// functions to compile an OpSum into a form which is understood by the basis
// implementation

OpSum compile(OpSum const &ops, Spinhalf const &block);

} // namespace xdiag::blocks
