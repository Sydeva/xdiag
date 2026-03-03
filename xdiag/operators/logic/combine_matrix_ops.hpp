// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <vector>

#include <xdiag/operators/op.hpp>

namespace xdiag::operators {
Op combine_matrix_ops(std::vector<Op> const &ops);
}
