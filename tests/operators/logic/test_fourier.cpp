// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include <cmath>

#include <xdiag/algebra/matrix.hpp>
#include <xdiag/blocks/electron.hpp>
#include <xdiag/operators/logic/fourier.hpp>
#include <xdiag/operators/logic/qns.hpp>

using namespace xdiag;

namespace {

PermutationGroup cyclic_group(int64_t nsites) {
  std::vector<Permutation> perms;
  perms.reserve((size_t)nsites);
  for (int64_t shift = 0; shift < nsites; ++shift) {
    std::vector<int64_t> p;
    p.reserve((size_t)nsites);
    for (int64_t site = 0; site < nsites; ++site) {
      p.push_back((site + shift) % nsites);
    }
    perms.push_back(Permutation(p));
  }
  return PermutationGroup(perms);
}

std::vector<Representation> cyclic_irreps(PermutationGroup const &group) {
  int64_t n = group.size();
  std::vector<Representation> irreps;
  irreps.reserve((size_t)n);

  for (int64_t k = 0; k < n; ++k) {
    std::vector<complex> chis;
    chis.reserve((size_t)n);
    double two_pi = 2.0 * std::acos(-1.0);
    for (int64_t l = 0; l < n; ++l) {
      double angle = two_pi * (double)(l * k) / (double)n;
      chis.push_back({std::cos(angle), std::sin(angle)});
    }
    irreps.push_back(Representation(group, chis));
  }

  return irreps;
}

} // namespace

TEST_CASE("fourier", "[operators]") try {
  int64_t nsites = 6;
  auto group = cyclic_group(nsites);
  auto irreps = cyclic_irreps(group);

  SECTION("arity smaller than two is rejected") {
    REQUIRE_THROWS(fourier(Op("Ntot", 0), {irreps[0]}));
  }

  SECTION("mixed arities are rejected") {
    OpSum ops;
    ops += Op("NtotNtot", {0, 1});
    ops += Op("CdagupCdagupCupCup", {0, 1, 2, 3});
    REQUIRE_THROWS(fourier(ops, {irreps[0], irreps[0]}));
  }

  SECTION("NtotNtot admissible and non-admissible momentum choices") {
    Op op("NtotNtot", {0, 1});

    // Ntot slots use plain characters, hence admissibility is k0 + k1 = 0 mod N.
    auto F_ok = fourier(op, {irreps[1], irreps[5]});
    auto r_ok = representation(F_ok, group);
    REQUIRE(r_ok);
    REQUIRE(isapprox(*r_ok, Representation(group)));

    auto block = Electron(nsites, 2, 2, Representation(group));
    REQUIRE(matrixC(F_ok, block).n_rows == (uint64_t)block.size());

    REQUIRE_THROWS(fourier(op, {irreps[1], irreps[1]}));
  }

  SECTION("CdagupCdagupCupCup uses conjugated characters on annihilation slots") {
    Op op("CdagupCdagupCupCup", {0, 1, 2, 3});

    // Slot convention from the apply kernel:
    // slots 0,1 are annihilation (conjugated chars), slots 2,3 are creation.
    // Admissibility: -k0 - k1 + k2 + k3 = 0 mod N.
    auto F_ok = fourier(op, {irreps[1], irreps[2], irreps[3], irreps[0]});
    auto r_ok = representation(F_ok, group);
    REQUIRE(r_ok);
    REQUIRE(isapprox(*r_ok, Representation(group)));

    auto block = Electron(nsites, 3, 0, Representation(group));
    REQUIRE(matrixC(F_ok, block).n_rows == (uint64_t)block.size());

    REQUIRE_THROWS(fourier(op, {irreps[1], irreps[2], irreps[0], irreps[0]}));
  }

  SECTION("CdagupCdagupCupCupHC uses conjugated characters on annihilation slots") {
    Op op("CdagupCdagupCupCupHC", {0, 1, 2, 3});

    // Slot convention from slot_character_mode:
    // slots 2,3 are annihilation (conjugated chars), slots 0,1 are creation.
    // Admissibility: k0 + k1 - k2 - k3 = 0 mod N.
    auto F_ok = fourier(op, {irreps[1], irreps[2], irreps[0], irreps[3]});
    auto r_ok = representation(F_ok, group);
    REQUIRE(r_ok);
    REQUIRE(isapprox(*r_ok, Representation(group)));

    auto block = Electron(nsites, 3, 0, Representation(group));
    REQUIRE(matrixC(F_ok, block).n_rows == (uint64_t)block.size());

    REQUIRE_THROWS(fourier(op, {irreps[1], irreps[2], irreps[0], irreps[0]}));
  }

} catch (xdiag::Error e) {
  xdiag::error_trace(e);
}


