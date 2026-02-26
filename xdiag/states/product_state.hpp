// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>
#include <vector>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

class ProductState {
public:
  using iterator_t = std::vector<std::string>::const_iterator;

  XDIAG_API ProductState() = default;
  XDIAG_API explicit ProductState(int64_t nsites);
  XDIAG_API explicit ProductState(std::vector<std::string> const &local_states);

  XDIAG_API std::string const &operator[](int64_t i) const;
  XDIAG_API std::string &operator[](int64_t i);

  XDIAG_API int64_t size() const;
  XDIAG_API int64_t nsites() const;
  XDIAG_API void push_back(std::string l);

  XDIAG_API iterator_t begin() const;
  XDIAG_API iterator_t end() const;

  XDIAG_API bool operator==(ProductState const &rhs) const;
  XDIAG_API bool operator!=(ProductState const &rhs) const;

private:
  std::vector<std::string> local_states_;
};

XDIAG_API int64_t size(ProductState const &p);
XDIAG_API int64_t nsites(ProductState const &p);
XDIAG_API std::ostream &operator<<(std::ostream &out,
                                   ProductState const &state);
XDIAG_API std::string to_string(ProductState const &state,
                                std::string format = "fancy");

} // namespace xdiag
