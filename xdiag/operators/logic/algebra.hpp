// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <set>
#include <string>
#include <vector>

#include <xdiag/operators/logic/rewrite.hpp>

namespace xdiag::operators {

// Algebra describes how to bring an OpSum into normal order for a specific
// physical system. It contains:
//
//   elementary_types  — the canonical generators (e.g. {S+, S-, Sz} for spin)
//   fermionic_types   — subset of elementary_types that anticommute at
//                       different sites (used to determine swap signs)
//   expansion_rules   — OpRules that expand compound operators into products
//                       of elementary ones (applied first, to fixed point)
//   algebra_rules     — MonomialRules that simplify same-site products and
//                       sort operators into canonical order (applied after
//                       expansion, to fixed point)
//
// The three pre-built algebras cover spin-1/2, electron, and tJ systems.
// Pass one to normal_order() to bring an OpSum into normal form.

struct Algebra {
  std::string name; // human-readable name for errors
  std::vector<std::string> elementary_types;
  std::set<std::string> fermionic_types;
  std::set<std::string>
      allowed_types; // all valid input types (checked before expansion)
  std::vector<OpRule> expansion_rules;
  std::vector<MonomialRule> algebra_rules;
};

Algebra spin_algebra();
Algebra electron_algebra();
Algebra tj_algebra();
Algebra matrix_algebra();
Algebra spinhalf_implementation_algebra();

} // namespace xdiag::operators
