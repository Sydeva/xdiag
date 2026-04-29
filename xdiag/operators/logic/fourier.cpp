// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

// Usage notes for fourier(Op/OpSum, irreps):
// - Performs slot-wise Fourier averaging over a common PermutationGroup.
// - Requires exactly one Representation per operator slot.
// - The resulting operator must transform trivially under the group so that it
//   maps an irrep-restricted block to itself.
// - Character weights are currently fermion-focused: annihilation slots use the
//   complex-conjugated character. At present this is implemented for Cup, Cdn,
//   CdagupCdagupCupCup, and CdagupCdagupCupCupHC.
// - For translation irreps this means admissibility is given by the signed
//   momentum sum implied by the creation/annihilation slot structure.
// - Current limitations: arity must be at least 2, OpSum inputs cannot mix
//   arities, all irreps must belong to the same group, and support is currently
//   intended for fermionic operators.

#include "fourier.hpp"

#include <optional>

#include <xdiag/operators/logic/order.hpp>
#include <xdiag/operators/logic/qns.hpp>
#include <xdiag/operators/logic/real.hpp>
#include <xdiag/operators/logic/valid.hpp>

namespace xdiag {

namespace {

enum class SlotCharacterMode { Plain, Conjugate };

int64_t arity(Op const &op) {
  if (!op.hassites()) {
    return 0;
  }
  return (int64_t)op.sites().size();
}

SlotCharacterMode slot_character_mode(Op const &op, int64_t slot) {
  auto const &type = op.type();

  if ((type == "Cup") || (type == "Cdn")) {
    return SlotCharacterMode::Conjugate;
  }

  if (type == "CdagupCdagupCupCup") {
    if ((slot == 0) || (slot == 1)) {
      return SlotCharacterMode::Conjugate;
    }
    return SlotCharacterMode::Plain;
  }

  if (type == "CdagupCdagupCupCupHC") {
    if ((slot == 2) || (slot == 3)) {
      return SlotCharacterMode::Conjugate;
    }
    return SlotCharacterMode::Plain;
  }

  return SlotCharacterMode::Plain;
}

template <typename T>
T character_for_slot(Op const &op, int64_t slot, T const &character) {
  if (slot_character_mode(op, slot) == SlotCharacterMode::Conjugate) {
    return xdiag::conj(character);
  }
  return character;
}

void check_irreps_groups(std::vector<Representation> const &irreps) {
  if (irreps.empty()) {
    XDIAG_THROW("Fourier transform needs at least two irreducible "
                "representations");
  }
  auto const &group = irreps[0].group();
  for (auto const &irrep : irreps) {
    if (irrep.group() != group) {
      XDIAG_THROW("Fourier transform requires all irreducible "
                  "representations to be defined for the same group");
    }
  }
}

void check_op_sites(Op const &op, PermutationGroup const &group) {
  if (!op.hassites()) {
    return;
  }

  for (int64_t s : op.sites()) {
    if (s < 0) {
      XDIAG_THROW(fmt::format("Cannot apply fourier transform: found a "
                              "negative site index \"{}\"\nOp:\n{}",
                              s, to_string(op)));
    }
    if (s >= group.nsites()) {
      XDIAG_THROW(fmt::format("Cannot apply fourier transform: Op has a site "
                              "with number \"{}\" which exceeds the number "
                              "of sites in the PermutationGroup \"{}\"\n"
                              "Op:\n{}\nGroup:\n{}",
                              s, group.nsites(), to_string(op),
                              to_string(group)));
    }
  }
}

bool generated_term_is_valid(Op const &op) {
  try {
    check_valid(op);
    return true;
  } catch (Error const &) {
    return false;
  }
}

template <typename T>
OpSum fourier_impl(OpSum const &ops, std::vector<Representation> const &irreps,
                   std::vector<arma::Col<T>> const &characters) {
  auto const &group = irreps[0].group();
  int64_t Ng = group.size();
  int64_t m = (int64_t)irreps.size();
  T norm = (T)1;
  for (int64_t i = 0; i < m; ++i) {
    norm /= (T)Ng;
  }

  OpSum out;

  for (auto const &[cpl, op] : ops.plain()) {
    auto const &sites = op.sites();
    std::vector<int64_t> tuple(m, 0);
    bool done = false;

    while (!done) {
      std::vector<int64_t> sites_new(m);
      T weight = cpl.scalar().as<T>() * norm;

      for (int64_t slot = 0; slot < m; ++slot) {
        int64_t gi = tuple[slot];
        sites_new[slot] = group[gi][sites[slot]];
        auto ch = characters[(size_t)slot](gi);
        weight *= character_for_slot(op, slot, ch);
      }

      if (op.hasmatrix()) {
        Op op_new(op.type(), sites_new, op.matrix());
        if (generated_term_is_valid(op_new)) {
          out += Scalar(weight) * op_new;
        }
      } else {
        Op op_new(op.type(), sites_new);
        if (generated_term_is_valid(op_new)) {
          out += Scalar(weight) * op_new;
        }
      }

      for (int64_t slot = m - 1; slot >= 0; --slot) {
        tuple[slot] += 1;
        if (tuple[slot] < Ng) {
          break;
        }
        tuple[slot] = 0;
        if (slot == 0) {
          done = true;
        }
      }
    }
  }

  out = order(out);
  auto repr = representation(out, group);
  if (!repr || !isapprox(*repr, Representation(group))) {
    XDIAG_THROW("Fourier result does not transform trivially and therefore "
                "does not map an irrep-restricted block to itself");
  }

  return out;
}

} // namespace

OpSum fourier(Op const &op, std::vector<Representation> const &irreps) try {
  return fourier(OpSum(op), irreps);
}
XDIAG_CATCH

OpSum fourier(OpSum const &ops, std::vector<Representation> const &irreps) try {
  check_irreps_groups(irreps);

  int64_t common_arity = -1;
  for (auto const &[cpl, op] : ops.plain()) {
    (void)cpl;
    auto op_arity = arity(op);
    if (common_arity < 0) {
      common_arity = op_arity;
    } else if (common_arity != op_arity) {
      XDIAG_THROW("fourier(OpSum, irreps) does not support mixed operator "
                  "arities");
    }
  }

  if (common_arity < 0) {
    XDIAG_THROW("fourier(OpSum, irreps) requires a non-empty OpSum");
  }

  if (common_arity < 2) {
    XDIAG_THROW("Fourier transform requires operator arity >= 2");
  }

  if ((int64_t)irreps.size() != common_arity) {
    XDIAG_THROW(fmt::format("Fourier transform requires exactly one irrep per "
                            "operator slot: got {} irreps for arity {}",
                            irreps.size(), common_arity));
  }

  auto const &group = irreps[0].group();
  for (auto const &[cpl, op] : ops.plain()) {
    (void)cpl;
    check_valid(op);
    check_op_sites(op, group);
  }

  bool all_real = isreal(ops);
  for (auto const &irrep : irreps) {
    all_real = all_real && isreal(irrep);
  }

  if (all_real) {
    std::vector<arma::vec> chars;
    chars.reserve(irreps.size());
    for (auto const &irrep : irreps) {
      chars.push_back(irrep.characters().as<arma::vec>());
    }
    return fourier_impl(ops, irreps, chars);
  }

  std::vector<arma::cx_vec> chars;
  chars.reserve(irreps.size());
  for (auto const &irrep : irreps) {
    chars.push_back(irrep.characters().as<arma::cx_vec>());
  }
  return fourier_impl(ops, irreps, chars);
}
XDIAG_CATCH

} // namespace xdiag

