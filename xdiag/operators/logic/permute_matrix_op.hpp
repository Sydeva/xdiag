// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/operators/op.hpp>

namespace xdiag::operators {

// Takes a "Matrix" Op with (potentially unsorted) sites and returns an
// equivalent Op with sites in ascending order, with the matrix permuted
// accordingly so that the operator is unchanged.
Op permute_matrix_op(Op const &op);

} // namespace xdiag::operators
