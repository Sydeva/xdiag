// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "valid.hpp"

#include <algorithm>

#include <xdiag/math/ipow.hpp>
#include <xdiag/operators/logic/valid.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

namespace xdiag::blocks {

void check_valid(Spinhalf const &block, OpSum const &ops) try {

  const std::vector<std::string> known_types = {
      "Sx",    "Sy",       "Sz",
      "Sp",    "Sm",       "SdotS",
      "SzSz",  "Exchange", "ScalarChirality",
      "Matrix"};

  std::string str;
  for (auto type : known_types) {
    str += fmt::format("\"{}\", ", type);
  }
  std::string known_types_string = str.substr(0, str.size() - 2);

  const std::map<std::string, int64_t> sites_for_type = {
      {"Sx", 1},   {"Sy", 1},       {"Sz", 1},
      {"Sp", 1},   {"Sm", 1},       {"SdotS", 2},
      {"SzSz", 2}, {"Exchange", 2}, {"ScalarChirality", 3}};

  for (auto const &[c, monomial] : ops) {
    for (auto const &op : monomial) {

      std::string type = op.type();

      // check if known type
      if (std::find(known_types.begin(), known_types.end(), type) ==
          known_types.end()) {
        XDIAG_THROW(
            fmt::format("Encountered invalid Op type \"{}\" for Spinhalf "
                        "block. Known types are:\n{}",
                        type, known_types_string));
      }

      // check whether sites of operators are valid
      if (op.hassites()) {
        must_have_sites_in_range(op, 0, block.nsites());
      }

      // check if Matrix type has proper dimension
      if (type == "Matrix") {
        must_have_matrix(op);
        must_have_disjoint_sites(op);
        Matrix mat = op.matrix();
        int64_t op_nsites = op.sites().size();
        int64_t expected_dim = math::ipow(2, op_nsites);
        if ((mat.n_rows() != expected_dim) || (mat.n_cols()) != expected_dim) {
          XDIAG_THROW(
              fmt::format("Encountered a \"Matrix\"-type Op on {} sites so "
                          "expected a 2^{} x 2^{} = {} x {} matrix. However, I "
                          "got a {} x {} matrix instead.",
                          op_nsites, op_nsites, op_nsites, expected_dim,
                          expected_dim, mat.n_rows(), mat.n_cols()));
        }
      } else {
        must_not_have_matrix(op);
        must_have_sites(op);
        must_have_disjoint_sites(op);
        must_have_nsites(op, sites_for_type.at(type));
      }
    }
  }
}
XDIAG_CATCH

} // namespace xdiag::blocks
