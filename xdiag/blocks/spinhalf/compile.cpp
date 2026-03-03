// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "compile.hpp"

#include <xdiag/operators/logic/rewrite.hpp>
#include <xdiag/spinhalf/valid.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

namespace xdiag::blocks {
OpSum compile(Spinhalf const &block, OpSum const &ops, std::string mode) try {
  check_valid(block, ops);

  // convert Sx and Sy operators to S+ and S-
  // Sx -> 0.5(Sp + Sm)
  // Sy -> -i*0.5(Sp - Sm)
  OpRule sxsy_rule = [](Op const &op) -> std::optional<OpSum> {
    if (op.type() == "Sx") {
      return 0.5 * Op("S+", op[0]) + 0.5 * Op("S-", op[0]);
    } else if (op.type() == "Sy") {
      return complex(0.0, -0.5) * Op("S+", op[0]) +
             complex(0.0, 0.5) * Op("S-", op[0]);
    }
    return std::nullopt;
  };

  OpRule sdots_rule = [](Op const &op) -> std::optional<OpSum> {
    if (op.type() == "SdotS") {
      return Op("SzSz", {op[0], op[1]}) + Op("Exchange", {op[0], op[1]});
    }
    return std::nullopt;
  };

  // turn every operator in to a matrix operator
  OpRule matrix_rule = [](Op const &op) -> std::optional<OpSum> {
    if (op.type() == "S+") {
      arma::mat M = {{0, 0}, {1, 0}};
      return Op("Matrix", op[0], mat);
    } else if (op.type() == "S-") {
      arma::mat M = {{0, 1}, {0, 0}};
      return Op("Matrix", op[0], mat);
    } else if (op.type() == "Sz") {
      arma::mat M = {{-0.5, 0}, {0, 0.5}};
      return Op("Matrix", op[0], mat);
    } else if (op.type() == "SzSz") {
      arma::mat M = {{0.25, 0.0, 0.0, 0.0},
                     {0.0, -0.25, 0.0, 0.0},
                     {0.0, 0.0, -0.25, 0.0},
                     {0.0, 0.0, 0.0, 0.25}};
      return Op("Matrix", {op[0], op[1]}, mat);
    } else if (op.type() == "Exchange") {
      arma::mat M = {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.5, 0.0},
                     {0.0, 0.5, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0};
      return Op("Matrix", {op[0], op[1]}, mat);
    } else if (op.type() == "ScalarChirality") {
      complex iq(0, 0.25);
      arma::cx_mat M = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                        {0.0, 0.0, iq, 0.0, -iq, 0.0, 0.0, 0.0},
                        {0.0, -iq, 0.0, 0.0, iq, 0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0, 0.0, 0.0, -iq, iq, 0.0},
                        {0.0, iq, -iq, 0.0, 0.0, 0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0, iq, 0.0, 0.0, -iq, 0.0},
                        {0.0, 0.0, 0.0, -iq, 0.0, iq, 0.0, 0.0},
                        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
      return Op("Matrix", {op[0], op[1], op[2]}, mat);
    }
    return std::nullopt;
  };

  if (mode == "implementation") {
    std::vector<OpRule> orules = {[](Op const &op) -> std::optional<OpSum> {
      if (op.type() == "Sx") {
        return OpSum(Op("Y", op.sites()));
      }
      return std::nullopt;
    }};
  } else if (mode == "matrix") {
  } else {
    XDIAG_THROW(fmt::format("Unknown compilation mode \"{}\". must be either "
                            "\"implementation\" or \"matrix\".",
                            mode));
  }
}
XDIAG_CATCH
} // namespace xdiag::blocks
