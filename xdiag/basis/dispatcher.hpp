// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <functional>
#include <xdiag/basis/basis.hpp>

namespace xdiag::basis {

class Dispatcher {
public:
  using func_t = std::function<void(Basis *, Basis *)>;
  Dispatcher() = default;

  // Register a function for a specific derived type
  template <typename T>
  void add(std::function<void(T const &, T const &)> fn) {
    table_[T::static_type()] = [fn](Basis *a, Basis *b) {
      fn(*static_cast<T *>(a), *static_cast<T *>(b));
    };
  }

  // Call the function for two shapes, throws if types differ or unregistered
  void dispatch(Basis *a, Basis *b) const;

private:
  std::unordered_map<std::size_t, func_t> table_;
};

} // namespace xdiag::basis
