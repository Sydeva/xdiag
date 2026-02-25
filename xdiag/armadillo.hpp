// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

// We stick to armadillo version 12.8.4 to minimize compiler constraints

#define ARMA_DONT_USE_WRAPPER

// Redefinition of "Op" to avoid name colission with xdiag::Op
#define Op XDIAG_OP_SUBSTITUTE
#include <xdiag/extern/armadillo/armadillo>
#undef Op
