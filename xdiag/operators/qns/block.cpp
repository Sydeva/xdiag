// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "block.hpp"

#include <xdiag/operators/qns/nup_ndn.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag {

Spinhalf block(OpSum const &ops, Spinhalf const &block) try {
  auto nup_block = block.nup();
  auto nup_ops = operators::nup(ops);
  if (nup_block && nup_ops) {
    return ((*nup_ops) == 0)
               ? block
               : Spinhalf(block.nsites(), (*nup_block) + (*nup_ops));
  } else { // block has no nup defined
    return block;
  }
}
XDIAG_CATCH

Block block(OpSum const &ops, Block const &blocki) try {
  return std::visit([&](auto &&b) { return Block(block(ops, b)); }, blocki);
}
XDIAG_CATCH

} // namespace xdiag
