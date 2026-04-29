// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/apply.hpp>
#include <xdiag/algebra/isapprox.hpp>
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/operators/logic/hc.hpp>
#include <xdiag/operators/logic/qns.hpp>
#include <xdiag/states/create_state.hpp>
#include <xdiag/states/state.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;

// Returns apply(Op("S+", si), apply(Op("S+", sj), psi))
// using Spinhalf blocks without fixed nup so that the nup sector change is
// handled correctly.
static State apply_sp_sp(int64_t si, int64_t sj, State const &psi) {
  auto step1 = apply(Op("S+", sj), psi);
  auto step2 = apply(Op("S+", si), step1);
  return step2;
}

// Returns apply(Op("S-", si), apply(Op("S-", sj), psi))
static State apply_sm_sm(int64_t si, int64_t sj, State const &psi) {
  auto step1 = apply(Op("S-", sj), psi);
  auto step2 = apply(Op("S-", si), step1);
  return step2;
}

TEST_CASE("spinhalf_spsp_smsm", "[spinhalf]") try {
  Log("Test S+S+ and S-S- operators against successive single-site S+/S-");

  // ── quantum-number checks ──────────────────────────────────────────────────
  REQUIRE(*nup(Op("S+S+", {0, 1})) == 2);
  REQUIRE(*nup(Op("S-S-", {0, 1})) == -2);

  // ── hermitian conjugate checks ─────────────────────────────────────────────
  REQUIRE(hc(Op("S+S+", {0, 1})) == Op("S-S-", {0, 1}));
  REQUIRE(hc(Op("S-S-", {0, 1})) == Op("S+S+", {0, 1}));

  // ── apply checks ──────────────────────────────────────────────────────────
  // Use small lattices; no fixed nup so that the sector change is allowed.
  for (int64_t nsites = 2; nsites <= 5; ++nsites) {
    auto block = Spinhalf(nsites); // all Sz sectors

    // random real state in the unrestricted block
    int64_t seed = 42;
    auto psi = random_state(block, /*real=*/true, /*ncols=*/1, seed);

    for (int64_t i = 0; i < nsites; ++i) {
      for (int64_t j = 0; j < nsites; ++j) {
        if (i == j)
          continue;

        // ── S+S+ ─────────────────────────────────────────────────────────────
        // built-in operator
        auto spsp_psi = apply(Op("S+S+", {i, j}), psi);
        // reference: successive single-site S+
        auto ref_spsp = apply_sp_sp(i, j, psi);
        REQUIRE(isapprox(spsp_psi.vector(), ref_spsp.vector(), 1e-12, 1e-12));

        // ── S-S- ─────────────────────────────────────────────────────────────
        auto smsm_psi = apply(Op("S-S-", {i, j}), psi);
        auto ref_smsm = apply_sm_sm(i, j, psi);
        REQUIRE(isapprox(smsm_psi.vector(), ref_smsm.vector(), 1e-12, 1e-12));

        // ── S-S- is h.c. of S+S+:  <phi|S+S+|psi> = <psi|S-S-|phi>* ────────
        auto phi = random_state(block, /*real=*/true, /*ncols=*/1, seed + 1);
        double lhs = dot(phi, spsp_psi);           // <phi|S+S+|psi>
        auto smsm_phi = apply(Op("S-S-", {i, j}), phi);
        double rhs = dot(psi, smsm_phi);           // <psi|S-S-|phi>
        REQUIRE(isapprox(lhs, rhs, 1e-12, 1e-12));
      }
    }
  }
} catch (xdiag::Error const &e) {
  xdiag::error_trace(e);
}

