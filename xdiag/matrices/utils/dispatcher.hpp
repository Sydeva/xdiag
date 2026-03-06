// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <functional>
#include <memory>
#include <xdiag/basis/basis.hpp>

namespace xdiag::matrices {

class Dispatcher {
public:
  using pointer_t = std::shared_ptr<basis::Basis>;
  using func_t = std::function<void(pointer_t const &, pointer_t const &)>;
  Dispatcher() = default;

  // Register a function for a specific derived type
  template <typename T> void add(std::function<void(T const &, T const &)> fn) {
    table_[T::static_type()] = [fn](pointer_t const &a, pointer_t const &b) {
      fn(static_cast<T const &>(*a), static_cast<T const &>(*b));
    };
  }

  // Call the function for two shapes, throws if types differ or unregistered
  void dispatch(pointer_t const &a, pointer_t const &b) const;

private:
  std::unordered_map<std::size_t, func_t> table_;
};

} // namespace xdiag::matrices
