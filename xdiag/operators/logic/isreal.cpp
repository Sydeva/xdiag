// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "isreal.hpp"

#include <string>
#include <vector>

#include <xdiag/operators/logic/types.hpp>
#include <xdiag/operators/logic/valid.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag {

bool isreal(Op const &op) try {
  check_valid(op);

  std::string type = op.type();
  if (op.hasmatrix()) {
    return op.matrix().isreal();
  } else {
    return is_real_type(type);
  }
}
XDIAG_CATCH

bool isreal(OpSum const &ops) try {
  for (auto const &[coeff, mono] : ops.plain()) {
    if (!isreal(coeff.scalar())) {
      return false;
    }
    for (auto const &op : mono) {
      if (!isreal(op)) {
        return false;
      }
    }
  }
  return true;
}
XDIAG_CATCH

} // namespace xdiag
