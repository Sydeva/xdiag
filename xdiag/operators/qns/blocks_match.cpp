// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "blocks_match.hpp"

#include <xdiag/operators/qns/nup_ndn.hpp>
#include <xdiag/utils/variants.hpp>

namespace xdiag {

bool blocks_match(OpSum const &ops, Block const &block1,
                  Block const &block2) try {
  return std::visit(utils::overload{
                        [&](Spinhalf const &b1, Spinhalf const &b2) {
                          return blocks_match(ops, b1, b2);
                        },
                        [&](auto const &b1, auto const &b2) { return false; },
                    },
                    block1, block2);
}
XDIAG_CATCH

bool blocks_match(OpSum const &ops, Spinhalf const &b1,
                  Spinhalf const &b2) try {
  if (b1.nup() && b2.nup()) {
    return (*b2.nup()) == (*b1.nup()) + (*operators::nup(ops));
  } else if (!b1.nup() && !b2.nup()) {
    return true;
  } else {
    return false;
  }
}
XDIAG_CATCH

} // namespace xdiag
