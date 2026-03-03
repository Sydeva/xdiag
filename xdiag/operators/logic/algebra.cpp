// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "algebra.hpp"

#include <xdiag/armadillo.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/operators/logic/combine_matrix_ops.hpp>
#include <xdiag/operators/logic/permute_matrix_op.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

namespace xdiag::operators {

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// Given a monomial and a position k of an adjacent pair (k, k+1) that has
// been matched by a rule, builds the full replacement OpSum:
//   prefix * repl_pair * suffix
// where prefix = mono[0..k-1], suffix = mono[k+2..end].
static OpSum replace_pair(Monomial const &mono, int64_t k,
                          OpSum const &repl_pair) {
  std::vector<Op> pre(mono.ops().begin(), mono.ops().begin() + k);
  std::vector<Op> suf(mono.ops().begin() + k + 2, mono.ops().end());
  Monomial prefix(pre), suffix(suf);

  OpSum result;
  for (auto const &[c, m] : repl_pair) {
    result += OpSum(c, prefix * m * suffix);
  }
  return result;
}

// Swap the pair at position (k, k+1) with a scalar sign.
static OpSum swap_pair(Monomial const &mono, int64_t k, double sign) {
  std::vector<Op> pre(mono.ops().begin(), mono.ops().begin() + k);
  std::vector<Op> suf(mono.ops().begin() + k + 2, mono.ops().end());
  Monomial prefix(pre), suffix(suf);
  Monomial swapped({mono[k + 1], mono[k]});
  return sign * (prefix * swapped * suffix);
}

// ---------------------------------------------------------------------------
// Common rules
// ---------------------------------------------------------------------------

// Remove a site-free "Id" op from a monomial when other operators surround it.
// A standalone {Id} monomial (length 1) is left unchanged — it represents the
// scalar identity and must not be collapsed to an empty monomial.
static MonomialRule id_absorption_rule() {
  return [](Monomial const &mono) -> std::optional<OpSum> {
    if (mono.size() <= 1)
      return std::nullopt;
    for (int64_t k = 0; k < mono.size(); ++k) {
      if (mono[k].type() == "Id" && !mono[k].hassites()) {
        std::vector<Op> ops;
        ops.reserve(mono.size() - 1);
        for (int64_t j = 0; j < mono.size(); ++j) {
          if (j != k)
            ops.push_back(mono[j]);
        }
        return OpSum(Monomial(ops));
      }
    }
    return std::nullopt;
  };
}

// Sort adjacent single-site operators by site index.
// sign_if_swap: -1 for fermionic (anticommuting) pairs, +1 for bosonic.
static MonomialRule sort_rule(std::set<std::string> const &fermionic_types) {
  return [fermionic_types](Monomial const &mono) -> std::optional<OpSum> {
    int64_t n = mono.size();
    for (int64_t k = 0; k + 1 < n; ++k) {
      Op const &a = mono[k];
      Op const &b = mono[k + 1];
      if (!a.hassites() || !b.hassites())
        continue;
      if (a.size() != 1 || b.size() != 1)
        continue;
      if (a[0] == b[0])
        continue; // same-site: handled by algebra rules, not here
      if (!(b < a))
        continue; // already in order

      bool a_fermi = fermionic_types.count(a.type()) > 0;
      bool b_fermi = fermionic_types.count(b.type()) > 0;
      double sign = (a_fermi && b_fermi) ? -1.0 : 1.0;
      return swap_pair(mono, k, sign);
    }
    return std::nullopt;
  };
}

// ---------------------------------------------------------------------------
// Spin-1/2 algebra
// ---------------------------------------------------------------------------

// Spin-1/2 same-site relations (adjacent pair, same site):
//   S+*S+ = 0         S-*S- = 0
//   S+*Sz = -1/2 S+   Sz*S+ = +1/2 S+
//   S-*Sz = +1/2 S-   Sz*S- = -1/2 S-
//   S+*S- = 1/2 Id + Sz
//   S-*S+ = 1/2 Id - Sz
//   Sz*Sz = 1/4 Id
static MonomialRule spin_same_site_rule() {
  return [](Monomial const &mono) -> std::optional<OpSum> {
    int64_t n = mono.size();
    for (int64_t k = 0; k + 1 < n; ++k) {
      Op const &a = mono[k];
      Op const &b = mono[k + 1];
      if (!a.hassites() || !b.hassites())
        continue;
      if (a.size() != 1 || b.size() != 1)
        continue;
      if (a[0] != b[0])
        continue;

      int64_t i = a[0];
      std::string ta = a.type(), tb = b.type();

      OpSum repl;
      if (ta == "S+" && tb == "S+") {
        repl = OpSum{}; // 0
      } else if (ta == "S-" && tb == "S-") {
        repl = OpSum{};
      } else if (ta == "S+" && tb == "Sz") {
        repl = -0.5 * Op("S+", i);
      } else if (ta == "Sz" && tb == "S+") {
        repl = 0.5 * Op("S+", i);
      } else if (ta == "S-" && tb == "Sz") {
        repl = 0.5 * Op("S-", i);
      } else if (ta == "Sz" && tb == "S-") {
        repl = -0.5 * Op("S-", i);
      } else if (ta == "S+" && tb == "S-") {
        repl = 0.5 * Op("Id") + OpSum(Op("Sz", i));
      } else if (ta == "S-" && tb == "S+") {
        repl = 0.5 * Op("Id") - OpSum(Op("Sz", i));
      } else if (ta == "Sz" && tb == "Sz") {
        repl = 0.25 * Op("Id");
      } else {
        continue;
      }
      return replace_pair(mono, k, repl);
    }
    return std::nullopt;
  };
}

// Spin expansion rules (compound -> products of S+, S-, Sz)
static std::vector<OpRule> spin_expansion_rules() {
  std::vector<OpRule> rules;

  // SdotS{i,j} -> Sz{i}*Sz{j} + 1/2 S+{i}*S-{j} + 1/2 S-{i}*S+{j}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "SdotS")
      return std::nullopt;
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += Op("Sz", i) * Op("Sz", j);
    r += 0.5 * (Op("S+", i) * Op("S-", j));
    r += 0.5 * (Op("S-", i) * Op("S+", j));
    return r;
  });

  // SzSz{i,j} -> Sz{i}*Sz{j}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "SzSz")
      return std::nullopt;
    int64_t i = op[0], j = op[1];
    return OpSum(Op("Sz", i) * Op("Sz", j));
  });

  // Exchange{i,j} -> 1/2 S+{i}*S-{j} + 1/2 S-{i}*S+{j}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Exchange")
      return std::nullopt;
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += 0.5 * (Op("S+", i) * Op("S-", j));
    r += 0.5 * (Op("S-", i) * Op("S+", j));
    return r;
  });

  // ScalarChirality{i,j,k} ->
  //   (i/2)*[ S+{i}*S-{j}*Sz{k} - S-{i}*S+{j}*Sz{k}
  //         + Sz{i}*S+{j}*S-{k} - Sz{i}*S-{j}*S+{k}
  //         + S-{i}*Sz{j}*S+{k} - S+{i}*Sz{j}*S-{k} ]
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "ScalarChirality")
      return std::nullopt;
    int64_t i = op[0], j = op[1], k = op[2];
    complex I(0.0, 1.0);
    OpSum r;
    r += (I * 0.5) * (Op("S+", i) * Op("S-", j) * Op("Sz", k));
    r += (-I * 0.5) * (Op("S-", i) * Op("S+", j) * Op("Sz", k));
    r += (I * 0.5) * (Op("Sz", i) * Op("S+", j) * Op("S-", k));
    r += (-I * 0.5) * (Op("Sz", i) * Op("S-", j) * Op("S+", k));
    r += (I * 0.5) * (Op("S-", i) * Op("Sz", j) * Op("S+", k));
    r += (-I * 0.5) * (Op("S+", i) * Op("Sz", j) * Op("S-", k));
    return r;
  });

  // Sx{i} -> 1/2 S+{i} + 1/2 S-{i}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Sx")
      return std::nullopt;
    int64_t i = op[0];
    return 0.5 * Op("S+", i) + 0.5 * Op("S-", i);
  });

  // Sy{i} -> -i/2 S+{i} + i/2 S-{i}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Sy")
      return std::nullopt;
    int64_t i = op[0];
    complex I(0.0, 1.0);
    return (-I * 0.5) * Op("S+", i) + (I * 0.5) * Op("S-", i);
  });

  return rules;
}

// ---------------------------------------------------------------------------
// Electron algebra
// ---------------------------------------------------------------------------

// Electron same-site CAR (adjacent pair, same site):
//   Cdagup*Cdagup = 0    Cup*Cup = 0
//   Cdagdn*Cdagdn = 0    Cdn*Cdn = 0
//   Cup*Cdagup    = Id - Cdagup*Cup
//   Cdn*Cdagdn    = Id - Cdagdn*Cdn
//   Cdagup*Cdagdn = -Cdagdn*Cdagup
//   Cdagdn*Cdagup = -Cdagup*Cdagdn   (already covered by above at same site)
//   Cup*Cdn       = -Cdn*Cup
//   Cdn*Cup       = -Cup*Cdn
//   Cdagup*Cdn    = -Cdn*Cdagup
//   Cdn*Cdagup    = -Cdagup*Cdn
//   Cup*Cdagdn    = -Cdagdn*Cup
//   Cdagdn*Cup    = -Cup*Cdagdn
static MonomialRule electron_same_site_rule() {
  return [](Monomial const &mono) -> std::optional<OpSum> {
    int64_t n = mono.size();
    for (int64_t k = 0; k + 1 < n; ++k) {
      Op const &a = mono[k];
      Op const &b = mono[k + 1];
      if (!a.hassites() || !b.hassites())
        continue;
      if (a.size() != 1 || b.size() != 1)
        continue;
      if (a[0] != b[0])
        continue;

      int64_t i = a[0];
      std::string ta = a.type(), tb = b.type();

      OpSum repl;
      if (ta == "Cdagup" && tb == "Cdagup") {
        repl = OpSum{};
      } else if (ta == "Cup" && tb == "Cup") {
        repl = OpSum{};
      } else if (ta == "Cdagdn" && tb == "Cdagdn") {
        repl = OpSum{};
      } else if (ta == "Cdn" && tb == "Cdn") {
        repl = OpSum{};
      } else if (ta == "Cup" && tb == "Cdagup") {
        // {Cup, Cdagup} = Id  =>  Cup*Cdagup = Id - Cdagup*Cup
        repl = OpSum(Op("Id")) - OpSum(Op("Cdagup", i) * Op("Cup", i));
      } else if (ta == "Cdn" && tb == "Cdagdn") {
        repl = OpSum(Op("Id")) - OpSum(Op("Cdagdn", i) * Op("Cdn", i));
        // Anticommutation rules: only the direction that moves toward canonical
        // alphabetical order (Cdagdn < Cdagup < Cdn < Cup). The reverse
        // directions are already canonical and must NOT fire (they would
        // cycle).
      } else if (ta == "Cdagup" && tb == "Cdagdn") {
        repl = -1.0 * (Op("Cdagdn", i) * Op("Cdagup", i)); // Cdagup>Cdagdn
      } else if (ta == "Cup" && tb == "Cdn") {
        repl = -1.0 * (Op("Cdn", i) * Op("Cup", i)); // Cup>Cdn
      } else if (ta == "Cdn" && tb == "Cdagup") {
        repl = -1.0 * (Op("Cdagup", i) * Op("Cdn", i)); // Cdn>Cdagup
      } else if (ta == "Cup" && tb == "Cdagdn") {
        repl = -1.0 * (Op("Cdagdn", i) * Op("Cup", i)); // Cup>Cdagdn
      } else {
        continue;
      }
      return replace_pair(mono, k, repl);
    }
    return std::nullopt;
  };
}

// Shared electron/tJ expansion rules
static std::vector<OpRule> fermionic_expansion_rules() {
  std::vector<OpRule> rules;

  // Hopup{i,j} -> -Cdagup{i}*Cup{j} - Cdagup{j}*Cup{i}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Hopup")
      return std::nullopt;
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += -1.0 * (Op("Cdagup", i) * Op("Cup", j));
    r += -1.0 * (Op("Cdagup", j) * Op("Cup", i));
    return r;
  });

  // Hopdn{i,j} -> -Cdagdn{i}*Cdn{j} - Cdagdn{j}*Cdn{i}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Hopdn")
      return std::nullopt;
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += -1.0 * (Op("Cdagdn", i) * Op("Cdn", j));
    r += -1.0 * (Op("Cdagdn", j) * Op("Cdn", i));
    return r;
  });

  // Hop{i,j} -> -Cdagup{i}*Cup{j} - Cdagup{j}*Cup{i}
  //             -Cdagdn{i}*Cdn{j} - Cdagdn{j}*Cdn{i}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Hop")
      return std::nullopt;
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += -1.0 * (Op("Cdagup", i) * Op("Cup", j));
    r += -1.0 * (Op("Cdagup", j) * Op("Cup", i));
    r += -1.0 * (Op("Cdagdn", i) * Op("Cdn", j));
    r += -1.0 * (Op("Cdagdn", j) * Op("Cdn", i));
    return r;
  });

  // Nup{i} -> Cdagup{i}*Cup{i}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Nup")
      return std::nullopt;
    int64_t i = op[0];
    return OpSum(Op("Cdagup", i) * Op("Cup", i));
  });

  // Ndn{i} -> Cdagdn{i}*Cdn{i}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Ndn")
      return std::nullopt;
    int64_t i = op[0];
    return OpSum(Op("Cdagdn", i) * Op("Cdn", i));
  });

  // Ntot{i} -> Cdagup{i}*Cup{i} + Cdagdn{i}*Cdn{i}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Ntot")
      return std::nullopt;
    int64_t i = op[0];
    OpSum r;
    r += Op("Cdagup", i) * Op("Cup", i);
    r += Op("Cdagdn", i) * Op("Cdn", i);
    return r;
  });

  // Nupdn{i} -> Cdagup{i}*Cup{i}*Cdagdn{i}*Cdn{i}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Nupdn")
      return std::nullopt;
    int64_t i = op[0];
    return OpSum(Op("Cdagup", i) * Op("Cup", i) * Op("Cdagdn", i) *
                 Op("Cdn", i));
  });

  // NtotNtot{i,j} -> Ntot{i}*Ntot{j}  (further expanded by Ntot rule)
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "NtotNtot")
      return std::nullopt;
    int64_t i = op[0], j = op[1];
    return OpSum(Op("Ntot", i) * Op("Ntot", j));
  });

  // NupdnNupdn{i,j} -> Nupdn{i}*Nupdn{j}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "NupdnNupdn")
      return std::nullopt;
    int64_t i = op[0], j = op[1];
    return OpSum(Op("Nupdn", i) * Op("Nupdn", j));
  });

  // NupNdn{i,j} -> Nup{i}*Ndn{j}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "NupNdn")
      return std::nullopt;
    int64_t i = op[0], j = op[1];
    return OpSum(Op("Nup", i) * Op("Ndn", j));
  });

  // NupNup{i,j} -> Nup{i}*Nup{j}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "NupNup")
      return std::nullopt;
    int64_t i = op[0], j = op[1];
    return OpSum(Op("Nup", i) * Op("Nup", j));
  });

  // NdnNdn{i,j} -> Ndn{i}*Ndn{j}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "NdnNdn")
      return std::nullopt;
    int64_t i = op[0], j = op[1];
    return OpSum(Op("Ndn", i) * Op("Ndn", j));
  });

  // NdnNup{i,j} -> Ndn{i}*Nup{j}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "NdnNup")
      return std::nullopt;
    int64_t i = op[0], j = op[1];
    return OpSum(Op("Ndn", i) * Op("Nup", j));
  });

  // Sz{i} -> 1/2 Nup{i} - 1/2 Ndn{i}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Sz")
      return std::nullopt;
    int64_t i = op[0];
    return 0.5 * Op("Nup", i) - 0.5 * Op("Ndn", i);
  });

  // S+{i} -> Cdagup{i} * Cdn{i}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "S+")
      return std::nullopt;
    int64_t i = op[0];
    return OpSum(Op("Cdagup", i) * Op("Cdn", i));
  });

  // S-{i} -> Cdagdn{i} * Cup{i}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "S-")
      return std::nullopt;
    int64_t i = op[0];
    return OpSum(Op("Cdagdn", i) * Op("Cup", i));
  });

  // Sx{i} -> 1/2 S+{i} + 1/2 S-{i}  (further expanded by S+/S- rules above)
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Sx")
      return std::nullopt;
    int64_t i = op[0];
    return 0.5 * Op("S+", i) + 0.5 * Op("S-", i);
  });

  // Sy{i} -> -i/2 S+{i} + i/2 S-{i}  (further expanded by S+/S- rules above)
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Sy")
      return std::nullopt;
    int64_t i = op[0];
    complex I(0.0, 1.0);
    return (-I * 0.5) * Op("S+", i) + (I * 0.5) * Op("S-", i);
  });

  return rules;
}

// ---------------------------------------------------------------------------
// tJ algebra
// ---------------------------------------------------------------------------

// tJ same-site CAR (projected, no double occupancy):
//   Cdagup*Cdagup = 0    Cup*Cup = 0
//   Cdagdn*Cdagdn = 0    Cdn*Cdn = 0
//   Cup*Cdagup  = Id - Cdagup*Cup - Cdagdn*Cdn   ({Cdagup,Cup} = Id - Ndn)
//   Cdn*Cdagdn  = Id - Cdagdn*Cdn - Cdagup*Cup   ({Cdagdn,Cdn} = Id - Nup)
//   Cdn*Cdagup  = 0   (no double occupancy: would create doublon)
//   Cup*Cdagdn  = 0   (no double occupancy)
//   Cdagup*Cdagdn = -Cdagdn*Cdagup
//   Cdagdn*Cdagup = -Cdagup*Cdagdn
//   Cup*Cdn = -Cdn*Cup
//   Cdn*Cup = -Cup*Cdn
//   Cdagup*Cdn = -Cdn*Cdagup
//   Cup*Cdagdn = -Cdagdn*Cup
static MonomialRule tj_same_site_rule() {
  return [](Monomial const &mono) -> std::optional<OpSum> {
    int64_t n = mono.size();
    for (int64_t k = 0; k + 1 < n; ++k) {
      Op const &a = mono[k];
      Op const &b = mono[k + 1];
      if (!a.hassites() || !b.hassites())
        continue;
      if (a.size() != 1 || b.size() != 1)
        continue;
      if (a[0] != b[0])
        continue;

      int64_t i = a[0];
      std::string ta = a.type(), tb = b.type();

      OpSum repl;
      if (ta == "Cdagup" && tb == "Cdagup") {
        repl = OpSum{};
      } else if (ta == "Cup" && tb == "Cup") {
        repl = OpSum{};
      } else if (ta == "Cdagdn" && tb == "Cdagdn") {
        repl = OpSum{};
      } else if (ta == "Cdn" && tb == "Cdn") {
        repl = OpSum{};
      } else if (ta == "Cup" && tb == "Cdagup") {
        // {Cdagup, Cup} = Id - Ndn  =>  Cup*Cdagup = Id - Cdagup*Cup -
        // Cdagdn*Cdn
        repl = OpSum(Op("Id")) - OpSum(Op("Cdagup", i) * Op("Cup", i)) -
               OpSum(Op("Cdagdn", i) * Op("Cdn", i));
      } else if (ta == "Cdn" && tb == "Cdagdn") {
        // {Cdagdn, Cdn} = Id - Nup  =>  Cdn*Cdagdn = Id - Cdagdn*Cdn -
        // Cdagup*Cup
        repl = OpSum(Op("Id")) - OpSum(Op("Cdagdn", i) * Op("Cdn", i)) -
               OpSum(Op("Cdagup", i) * Op("Cup", i));
      } else if (ta == "Cdn" && tb == "Cdagup") {
        repl = OpSum{}; // 0 (no double occupancy)
      } else if (ta == "Cup" && tb == "Cdagdn") {
        repl = OpSum{}; // 0 (no double occupancy)
        // Anticommutation rules: only the out-of-canonical-order direction.
        // Canonical order: Cdagdn < Cdagup < Cdn < Cup (alphabetical).
        // Cdn*Cdagup and Cup*Cdagdn are zero in tJ (handled above).
      } else if (ta == "Cdagup" && tb == "Cdagdn") {
        repl = -1.0 * (Op("Cdagdn", i) * Op("Cdagup", i)); // Cdagup>Cdagdn
      } else if (ta == "Cup" && tb == "Cdn") {
        repl = -1.0 * (Op("Cdn", i) * Op("Cup", i)); // Cup>Cdn
      } else {
        continue;
      }
      return replace_pair(mono, k, repl);
    }
    return std::nullopt;
  };
}

// tJ also needs spin exchange expansion (in terms of fermionic operators)
// Exchange{i,j} = S+{i}*S-{j} + S-{i}*S+{j} ... but for tJ,
// S+{i} = Cdagup{i}*Cdn{i}  and  S-{i} = Cdagdn{i}*Cup{i}
static std::vector<OpRule> tj_expansion_rules() {
  auto rules = fermionic_expansion_rules();

  // tJSzSz{i,j} -> (Nup{i}-Ndn{i})/2 * (Nup{j}-Ndn{j})/2
  //  = 1/4*(Nup{i}*Nup{j} - Nup{i}*Ndn{j} - Ndn{i}*Nup{j} + Ndn{i}*Ndn{j})
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "tJSzSz")
      return std::nullopt;
    int64_t i = op[0], j = op[1];
    OpSum r;
    r += 0.25 * (Op("Nup", i) * Op("Nup", j));
    r += -0.25 * (Op("Nup", i) * Op("Ndn", j));
    r += -0.25 * (Op("Ndn", i) * Op("Nup", j));
    r += 0.25 * (Op("Ndn", i) * Op("Ndn", j));
    return r;
  });

  // tJSdotS{i,j} -> tJSzSz{i,j} + Exchange{i,j}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "tJSdotS")
      return std::nullopt;
    int64_t i = op[0], j = op[1];
    return OpSum(Op("tJSzSz", {i, j})) + OpSum(Op("Exchange", {i, j}));
  });

  // Exchange{i,j} in tJ: S+{i}*S-{j} + S-{i}*S+{j}
  //   S+{i} = Cdagup{i}*Cdn{i},  S-{i} = Cdagdn{i}*Cup{i}
  rules.push_back([](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "Exchange")
      return std::nullopt;
    int64_t i = op[0], j = op[1];
    OpSum r;
    // S+{i}*S-{j} = Cdagup{i}*Cdn{i}*Cdagdn{j}*Cup{j}
    r += Op("Cdagup", i) * Op("Cdn", i) * Op("Cdagdn", j) * Op("Cup", j);
    // S-{i}*S+{j} = Cdagdn{i}*Cup{i}*Cdagup{j}*Cdn{j}
    r += Op("Cdagdn", i) * Op("Cup", i) * Op("Cdagup", j) * Op("Cdn", j);
    return r;
  });

  return rules;
}

// ---------------------------------------------------------------------------
// Public factory functions
// ---------------------------------------------------------------------------

Algebra spin_algebra() {
  std::set<std::string> fermionic{}; // spin operators commute

  auto algebra_rules_vec = std::vector<MonomialRule>{
      spin_same_site_rule(),
      sort_rule(fermionic),
      id_absorption_rule(),
  };

  return Algebra{
      .name = "spin-1/2",
      .elementary_types = {"S+", "S-", "Sz"},
      .fermionic_types = fermionic,
      .allowed_types = {"Exchange", "Id", "S+", "S-", "ScalarChirality",
                        "SdotS", "Sx", "Sy", "Sz", "SzSz"},
      .expansion_rules = spin_expansion_rules(),
      .algebra_rules = algebra_rules_vec,
  };
}

Algebra electron_algebra() {
  std::set<std::string> fermionic{"Cdagup", "Cup", "Cdagdn", "Cdn"};

  auto algebra_rules_vec = std::vector<MonomialRule>{
      electron_same_site_rule(),
      sort_rule(fermionic),
      id_absorption_rule(),
  };

  return Algebra{
      .name = "electron",
      .elementary_types = {"Cdagup", "Cup", "Cdagdn", "Cdn"},
      .fermionic_types = fermionic,
      .allowed_types = {"Cdagdn", "Cdagup", "Cdn",      "Cup",        "Hop",
                        "Hopdn",  "Hopup",  "HubbardU", "Id",         "Ndn",
                        "NdnNdn", "NdnNup", "Ntot",     "NtotNtot",   "Nup",
                        "NupNdn", "NupNup", "Nupdn",    "NupdnNupdn", "S+",
                        "S-",     "Sx",     "Sy",       "Sz"},
      .expansion_rules = fermionic_expansion_rules(),
      .algebra_rules = algebra_rules_vec,
  };
}

Algebra tj_algebra() {
  std::set<std::string> fermionic{"Cdagup", "Cup", "Cdagdn", "Cdn"};

  auto algebra_rules_vec = std::vector<MonomialRule>{
      tj_same_site_rule(),
      sort_rule(fermionic),
      id_absorption_rule(),
  };

  return Algebra{
      .name = "tJ",
      .elementary_types = {"Cdagup", "Cup", "Cdagdn", "Cdn"},
      .fermionic_types = fermionic,
      .allowed_types = {"Cdagdn", "Cdagup", "Cdn",   "Cup",        "Exchange",
                        "Hop",    "Hopdn",  "Hopup", "Id",         "Ndn",
                        "NdnNdn", "NdnNup", "Ntot",  "NtotNtot",   "Nup",
                        "NupNdn", "NupNup", "Nupdn", "NupdnNupdn", "S+",
                        "S-",     "Sx",     "Sy",    "Sz",         "tJSdotS",
                        "tJSzSz"},
      .expansion_rules = tj_expansion_rules(),
      .algebra_rules = algebra_rules_vec,
  };
}

// ---------------------------------------------------------------------------
// Matrix algebra  (and shared matrix-building MonomialRules)
// ---------------------------------------------------------------------------

// MonomialRule: combine the first pair of adjacent "Matrix" ops (with sites)
// into a single "Matrix" op using combine_matrix_ops.
static MonomialRule combine_matrix_rule() {
  return [](Monomial const &mono) -> std::optional<OpSum> {
    for (int64_t k = 0; k + 1 < mono.size(); ++k) {
      if (mono[k].type() == "Matrix" && mono[k].hassites() &&
          mono[k + 1].type() == "Matrix" && mono[k + 1].hassites()) {
        Op combined = combine_matrix_ops({mono[k], mono[k + 1]});
        std::vector<Op> ops;
        for (int64_t j = 0; j < k; ++j)
          ops.push_back(mono[j]);
        ops.push_back(combined);
        for (int64_t j = k + 2; j < mono.size(); ++j)
          ops.push_back(mono[j]);
        return OpSum(Monomial(ops));
      }
    }
    return std::nullopt;
  };
}

// MonomialRule: if a "Matrix" op has unsorted sites, sort them and permute
// the matrix accordingly using permute_matrix_op.
static MonomialRule sort_matrix_sites_rule() {
  return [](Monomial const &mono) -> std::optional<OpSum> {
    for (int64_t k = 0; k < mono.size(); ++k) {
      if (mono[k].type() != "Matrix" || !mono[k].hassites())
        continue;
      auto const &sites = mono[k].sites();
      bool sorted = true;
      for (int64_t j = 0; j + 1 < (int64_t)sites.size(); ++j) {
        if (sites[j] > sites[j + 1]) {
          sorted = false;
          break;
        }
      }
      if (!sorted) {
        Op permuted = permute_matrix_op(mono[k]);
        std::vector<Op> ops;
        for (int64_t j = 0; j < mono.size(); ++j)
          ops.push_back(j == k ? permuted : mono[j]);
        return OpSum(Monomial(ops));
      }
    }
    return std::nullopt;
  };
}

// ---------------------------------------------------------------------------
// SpinHalf implementation algebra
// ---------------------------------------------------------------------------

// Converts a single Op to an equivalent "Matrix" Op carrying its explicit
// matrix.  Used by spinhalf_to_matrix_rule when a non-Matrix op appears
// inside a multi-op monomial.
//
// Site/kron convention (mirrors combine_matrix_ops / embed_op):
//   For sites [s0, s1, ...], bit j of the local state index encodes the spin
//   at sites[j].  In armadillo kron(A, B):  A acts on the high (outer) bit,
//   B on the low (inner) bit.  So for sites [i, j]:
//     kron(outer=j_op, inner=i_op)
//   and for sites [i, j, k]:
//     kron(k_op, kron(j_op, i_op))
static Op op_to_matrix_op(Op const &op) try {
  static const arma::mat sp = {{0., 1.}, {0., 0.}};
  static const arma::mat sm = {{0., 0.}, {1., 0.}};
  static const arma::mat sz = {{0.5, 0.}, {0., -0.5}};

  auto const &type = op.type();

  if (type == "Matrix")
    return op;
  if (type == "S+")
    return Op("Matrix", op.sites(), sp);
  if (type == "S-")
    return Op("Matrix", op.sites(), sm);
  if (type == "Sz")
    return Op("Matrix", op.sites(), sz);

  // Sx{i} -> [[0,0.5],[0.5,0]]
  if (type == "Sx") {
    arma::mat m = {{0.0, 0.5}, {0.5, 0.0}};
    return Op("Matrix", op.sites(), m);
  }

  // Sy{i} -> [[0,-i/2],[i/2,0]]
  if (type == "Sy") {
    complex I(0.0, 1.0);
    arma::cx_mat m = {{0.0 + 0.0 * I, -I * 0.5}, {I * 0.5, 0.0 + 0.0 * I}};
    return Op("Matrix", op.sites(), m);
  }

  // Two-site compound ops.
  // When i == j (same site) the matrix is a plain 2×2 product; otherwise it
  // is a 4×4 Kronecker product.  Convention for different sites:
  //   kron(outer=site1, inner=site0) — matches combine_matrix_ops bit-index.

  // SdotS{i,j} = Sz_i Sz_j + 0.5*S+_i S-_j + 0.5*S-_i S+_j
  if (type == "SdotS") {
    int64_t i = op.sites()[0], j = op.sites()[1];
    if (i == j) {
      arma::mat m = arma::mat(sz * sz) + 0.5 * arma::mat(sp * sm) +
                    0.5 * arma::mat(sm * sp);
      return Op("Matrix", std::vector<int64_t>{i}, m);
    }
    arma::mat m = arma::mat(arma::kron(sz, sz)) +
                  0.5 * arma::mat(arma::kron(sm, sp)) +
                  0.5 * arma::mat(arma::kron(sp, sm));
    return Op("Matrix", op.sites(), m);
  }

  if (type == "SzSz") {
    int64_t i = op.sites()[0], j = op.sites()[1];
    if (i == j)
      return Op("Matrix", std::vector<int64_t>{i}, arma::mat(sz * sz));
    return Op("Matrix", op.sites(), arma::mat(arma::kron(sz, sz)));
  }
  if (type == "Exchange") {
    int64_t i = op.sites()[0], j = op.sites()[1];
    if (i == j) {
      arma::mat m = 0.5 * arma::mat(sp * sm) + 0.5 * arma::mat(sm * sp);
      return Op("Matrix", std::vector<int64_t>{i}, m);
    }
    arma::mat m = 0.5 * arma::mat(arma::kron(sm, sp)) +
                  0.5 * arma::mat(arma::kron(sp, sm));
    return Op("Matrix", op.sites(), m);
  }

  // Three-site ScalarChirality → 8×8 complex matrix.
  // S_i·(S_j×S_k) = (i/2)*[ S+_i S-_j Sz_k - S-_i S+_j Sz_k
  //                        + Sz_i S+_j S-_k - Sz_i S-_j S+_k
  //                        + S-_i Sz_j S+_k - S+_i Sz_j S-_k ]
  // Each product A_i B_j C_k maps to kron(C_k, kron(B_j, A_i)).
  if (type == "ScalarChirality") {
    complex I(0., 1.);
    arma::cx_mat sp_cx(sp, arma::zeros<arma::mat>(2, 2));
    arma::cx_mat sm_cx(sm, arma::zeros<arma::mat>(2, 2));
    arma::cx_mat sz_cx(sz, arma::zeros<arma::mat>(2, 2));
    // t(k_op, j_op, i_op) = kron(k_op, kron(j_op, i_op))
    auto t = [](arma::cx_mat const &a, arma::cx_mat const &b,
                arma::cx_mat const &c) {
      return arma::cx_mat(arma::kron(a, arma::cx_mat(arma::kron(b, c))));
    };
    arma::cx_mat m = (I * 0.5) * t(sz_cx, sm_cx, sp_cx)     // S+_i S-_j Sz_k
                     + (-I * 0.5) * t(sz_cx, sp_cx, sm_cx)  // S-_i S+_j Sz_k
                     + (I * 0.5) * t(sm_cx, sp_cx, sz_cx)   // Sz_i S+_j S-_k
                     + (-I * 0.5) * t(sp_cx, sm_cx, sz_cx)  // Sz_i S-_j S+_k
                     + (I * 0.5) * t(sp_cx, sz_cx, sm_cx)   // S-_i Sz_j S+_k
                     + (-I * 0.5) * t(sm_cx, sz_cx, sp_cx); // S+_i Sz_j S-_k
    return Op("Matrix", op.sites(), m);
  }

  XDIAG_THROW(
      fmt::format("Cannot convert Op of type \"{}\" to a Matrix op", type));
}
XDIAG_CATCH

// MonomialRule: a size-1 SdotS{i,j} monomial expands to Exchange{i,j} +
// SzSz{i,j}.  Used by spinhalf_implementation to keep those two types as
// named operators rather than collapsing them to a generic Matrix.
static MonomialRule spinhalf_sdots_rule() {
  return [](Monomial const &mono) -> std::optional<OpSum> {
    if (mono.size() != 1 || mono[0].type() != "SdotS")
      return std::nullopt;
    int64_t i = mono[0][0], j = mono[0][1];
    return OpSum(Op("Exchange", std::vector<int64_t>{i, j})) +
           OpSum(Op("SzSz", std::vector<int64_t>{i, j}));
  };
}

// MonomialRule: convert the first non-Matrix, non-Id op in a monomial to its
// Matrix form via op_to_matrix_op.
//
// protected_single_op_types — op types that are left as-is when they are the
//   sole op in a size-1 monomial (used by spinhalf_implementation to keep
//   "Sz", "S+", "S-", "SzSz", "Exchange", "ScalarChirality" in named form).
//   Pass an empty set for matrix_algebra, where every op becomes a Matrix.
//
// Iterated to fixed point together with combine_matrix_rule so all ops in
// multi-op monomials eventually collapse to one Matrix op.
static MonomialRule
convert_to_matrix_rule(std::set<std::string> const &protected_single_op_types) {
  return [protected_single_op_types](
             Monomial const &mono) -> std::optional<OpSum> {
    for (int64_t k = 0; k < mono.size(); ++k) {
      auto const &type = mono[k].type();
      if (type == "Matrix")
        continue;
      if (type == "Id")
        continue; // handled by id_absorption_rule
      // In a size-1 monomial, protected types stay as-is.
      if (mono.size() == 1 && protected_single_op_types.count(type) > 0)
        return std::nullopt;
      Op mat_op = op_to_matrix_op(mono[k]);
      std::vector<Op> ops;
      for (int64_t j = 0; j < mono.size(); ++j)
        ops.push_back(j == k ? mat_op : mono[j]);
      return OpSum(Monomial(ops));
    }
    return std::nullopt;
  };
}

Algebra spinhalf_implementation_algebra() {
  std::set<std::string> fermionic{};

  // Types that remain in named form when they appear as a size-1 monomial.
  // All other allowed types (Sx, Sy, SdotS after sdots_rule, multi-op
  // combinations) are converted to Matrix.
  static const std::set<std::string> protected_types = {
      "Exchange", "S+", "S-", "ScalarChirality", "Sz", "SzSz"};

  auto algebra_rules_vec = std::vector<MonomialRule>{
      id_absorption_rule(),
      spinhalf_sdots_rule(),
      convert_to_matrix_rule(protected_types),
      combine_matrix_rule(),
      sort_matrix_sites_rule(),
  };

  return Algebra{
      .name = "spinhalf_implementation_algebra",
      .elementary_types = {"Exchange", "Matrix", "S+", "S-", "ScalarChirality",
                           "Sz", "SzSz"},
      .fermionic_types = fermionic,
      .allowed_types = {"Exchange", "Id", "Matrix", "S+", "S-",
                        "ScalarChirality", "SdotS", "Sx", "Sy", "Sz", "SzSz"},
      .expansion_rules = {},
      .algebra_rules = algebra_rules_vec,
  };
}

Algebra matrix_algebra() {
  std::set<std::string> fermionic{}; // spin operators are bosonic

  // No types are protected: every op (including size-1 Sz, S+, etc.) is
  // converted to an explicit Matrix.  op_to_matrix_op handles all spin-1/2
  // types directly, so no OpRules are needed.
  auto algebra_rules_vec = std::vector<MonomialRule>{
      id_absorption_rule(),
      convert_to_matrix_rule({}),
      combine_matrix_rule(),
      sort_matrix_sites_rule(),
  };

  return Algebra{
      .name = "matrix",
      .elementary_types = {"Matrix"},
      .fermionic_types = fermionic,
      .allowed_types = {"Exchange", "Id", "Matrix", "S+", "S-",
                        "ScalarChirality", "SdotS", "Sx", "Sy", "Sz", "SzSz"},
      .expansion_rules = {},
      .algebra_rules = algebra_rules_vec,
  };
}

} // namespace xdiag::operators
