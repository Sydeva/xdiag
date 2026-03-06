// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "dispatcher.hpp"

#include <xdiag/utils/error.hpp>

namespace xdiag::matrices {

void Dispatcher::dispatch(typename Dispatcher::pointer_t const &a,
                          typename Dispatcher::pointer_t const &b) const try {
  if (a->type() != b->type()) {
    XDIAG_THROW("Type mismatch for Basis: \n  1 -> " + std::string(a->name()) +
                "\n  2 -> " + std::string(b->name()));
  }
  auto it = table_.find(a->type());
  if (it == table_.end()) {
    XDIAG_THROW("No handler registered for type: " + std::string(a->name()));
  }
  it->second(a, b);
}
XDIAG_CATCH

} // namespace xdiag::matrices
