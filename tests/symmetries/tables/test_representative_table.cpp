// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include <cmath>
#include <vector>

#include <xdiag/bits/bitarray.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/combinatorics/bounded_multisets/bounded_multisets.hpp>
#include <xdiag/combinatorics/bounded_partitions/bounded_partitions.hpp>
#include <xdiag/combinatorics/bounded_partitions/schaefer_table.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/config.hpp>
#include <xdiag/io/file_toml.hpp>
#include <xdiag/io/read.hpp>
#include <xdiag/symmetries/action/isrepresentative.hpp>
#include <xdiag/symmetries/action/site_permutation.hpp>
#include <xdiag/symmetries/cyclic_group.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/tables/representative_table.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;
using namespace xdiag::bits;
using namespace xdiag::combinatorics;
using namespace xdiag::symmetries;

// -----------------------------------------------------------------------
// Generic invariant checker: verifies structural correctness for any
// enumeration type without depending on specific representative values.
//
// Checks:
//   1. size() in [0, enumeration.size()]
//   2. All iterated representatives satisfy isrepresentative()
//   3. Symmetry property: apply(symmetry(idx), state) == representative(idx)
//   4. index_of_representative is in-bounds
//   5. All norms are positive
//   6. Two identically-constructed tables compare equal
// -----------------------------------------------------------------------
template <typename enumeration_t>
void test_invariants(enumeration_t const &enumeration,
                     SitePermutation const &sp,
                     Representation const &irrep) {
  using bit_t = typename enumeration_t::bit_t;
  RepresentativeTable<enumeration_t> table(enumeration, sp, irrep);

  REQUIRE(table.size() >= 0);
  REQUIRE(table.size() <= enumeration.size());

  for (auto rep : table) {
    REQUIRE(isrepresentative(rep, sp));
  }

  // Symmetry property: only valid for states in the orbit of a nonzero-norm
  // representative. Zero-norm states are not included in the table and their
  // representative_index / representative_symmetry fields remain at their
  // default values, so we must not check them here.
  for (auto rep : table) {
    for (int64_t sym = 0; sym < sp.size(); ++sym) {
      bit_t state = sp.apply(sym, rep);
      int64_t idx = enumeration.index(state);
      REQUIRE(table.representative(idx) == rep);
      REQUIRE(sp.apply(table.symmetry(idx), state) == rep);
      int64_t rep_idx = table.index_of_representative(idx);
      REQUIRE(rep_idx >= 0);
      REQUIRE(rep_idx < table.size());
    }
  }

  for (int64_t ri = 0; ri < table.size(); ++ri) {
    REQUIRE(table.norm(ri) > 0.0);
  }

  RepresentativeTable<enumeration_t> table2(enumeration, sp, irrep);
  REQUIRE(table == table2);
  REQUIRE_FALSE(table != table2);
}

// -----------------------------------------------------------------------
// Size-check helper: constructs a table and asserts the number of reps.
// -----------------------------------------------------------------------
template <typename enumeration_t>
void check_size(enumeration_t const &enumeration, SitePermutation const &sp,
                Representation const &irrep, int64_t expected_size) {
  RepresentativeTable<enumeration_t> table(enumeration, sp, irrep);
  REQUIRE(table.size() == expected_size);
}

// -----------------------------------------------------------------------
// TEST CASE
// -----------------------------------------------------------------------

TEST_CASE("representative_table", "[symmetries]") try {
  Log("Test RepresentativeTable");

  // Groups and irreps used throughout
  auto group4 = cyclic_group(4);
  auto sp4 = SitePermutation(group4);
  auto irrep4_0 = cyclic_group_irrep(4, 0); // trivial: all chi=1
  auto irrep4_1 = cyclic_group_irrep(4, 1); // complex
  auto irrep4_2 = cyclic_group_irrep(4, 2); // real: chi={1,-1,1,-1}

  auto group3 = cyclic_group(3);
  auto sp3 = SitePermutation(group3);
  auto irrep3_0 = cyclic_group_irrep(3, 0);
  auto irrep3_1 = cyclic_group_irrep(3, 1); // complex

  // =====================================================================
  // Subsets — representatives of all n-bit strings under cyclic_group(4)
  //
  // cyclic_group(4) on Subsets(4) — 16 states, 6 orbits:
  //   {0}, {1,2,4,8}, {3,6,12,9}, {5,10}, {7,14,13,11}, {15}
  // k=0: 6 reps (all orbits have nonzero norm)
  // k=2: 4 reps (orbits {0} and {15} excluded: amplitude=0)
  // k=1: 3 reps (orbits {0}, {5,10}, {15} excluded)
  // =====================================================================

  SECTION("Subsets<uint32_t>") {
    auto subsets = Subsets<uint32_t>(4);
    test_invariants(subsets, sp4, irrep4_0);
    test_invariants(subsets, sp4, irrep4_2);
    test_invariants(subsets, sp4, irrep4_1);
    check_size(subsets, sp4, irrep4_0, 6);
    check_size(subsets, sp4, irrep4_2, 4);
    check_size(subsets, sp4, irrep4_1, 3);

    // Verify k=0 representatives appear in iteration order
    std::vector<uint32_t> expected = {0, 1, 3, 5, 7, 15};
    RepresentativeTable<Subsets<uint32_t>> table(subsets, sp4, irrep4_0);
    std::vector<uint32_t> got;
    for (auto r : table)
      got.push_back(r);
    REQUIRE(got == expected);

    // Verify norm values for k=0: norm = sqrt(|stabilizer|)
    // rep 0  (orbit size 1): norm = sqrt(4) = 2
    // rep 1  (orbit size 4): norm = 1
    // rep 3  (orbit size 4): norm = 1
    // rep 5  (orbit size 2): stabilizer={id,r2}, norm = sqrt(2)
    // rep 7  (orbit size 4): norm = 1
    // rep 15 (orbit size 1): norm = sqrt(4) = 2
    REQUIRE(table.norm(0) == Approx(2.0));
    REQUIRE(table.norm(1) == Approx(1.0));
    REQUIRE(table.norm(2) == Approx(1.0));
    REQUIRE(table.norm(3) == Approx(std::sqrt(2.0)));
    REQUIRE(table.norm(4) == Approx(1.0));
    REQUIRE(table.norm(5) == Approx(2.0));

    // operator!= with a different irrep
    RepresentativeTable<Subsets<uint32_t>> table2(subsets, sp4, irrep4_2);
    REQUIRE(table != table2);
    REQUIRE_FALSE(table == table2);
  }

  SECTION("Subsets<uint64_t>") {
    auto subsets = Subsets<uint64_t>(4);
    test_invariants(subsets, sp4, irrep4_0);
    check_size(subsets, sp4, irrep4_0, 6);
  }

  // =====================================================================
  // Combinations (native integer types)
  //
  // cyclic_group(4) on Combinations(4,2) — 6 states:
  //   {3(0011), 5(0101), 6(0110), 9(1001), 10(1010), 12(1100)}
  // Orbits: {3,6,12,9} and {5,10} → 2 representatives
  // =====================================================================

  SECTION("Combinations<uint32_t>") {
    auto combos = Combinations<uint32_t>(4, 2);
    test_invariants(combos, sp4, irrep4_0);
    test_invariants(combos, sp4, irrep4_2);
    check_size(combos, sp4, irrep4_0, 2);
    check_size(combos, sp4, irrep4_2, 2);

    // Verify representatives via iterator
    RepresentativeTable<Combinations<uint32_t>> table(combos, sp4, irrep4_0);
    auto it = table.begin();
    REQUIRE(*it == 3u);
    ++it;
    REQUIRE(*it == 5u);
    ++it;
    REQUIRE(it == table.end());

    // Orbit {3,6,12,9}: all four states map to rep=3
    for (uint32_t s : {3u, 6u, 9u, 12u})
      REQUIRE(table.representative(combos.index(s)) == 3u);
    // Orbit {5,10}: both map to rep=5
    for (uint32_t s : {5u, 10u})
      REQUIRE(table.representative(combos.index(s)) == 5u);
  }

  SECTION("Combinations<uint64_t>") {
    auto combos = Combinations<uint64_t>(4, 2);
    test_invariants(combos, sp4, irrep4_0);
    check_size(combos, sp4, irrep4_0, 2);
  }

  // =====================================================================
  // Combinations (Bitset types) — same 6 states as uint32_t but packed
  // into wider bitsets; orbit structure and rep count unchanged.
  // =====================================================================

  SECTION("Combinations<BitsetStatic2>") {
    auto combos = Combinations<BitsetStatic2>(4, 2);
    test_invariants(combos, sp4, irrep4_0);
    check_size(combos, sp4, irrep4_0, 2);
  }

  SECTION("Combinations<BitsetStatic4>") {
    auto combos = Combinations<BitsetStatic4>(4, 2);
    test_invariants(combos, sp4, irrep4_0);
    check_size(combos, sp4, irrep4_0, 2);
  }

  SECTION("Combinations<BitsetStatic8>") {
    auto combos = Combinations<BitsetStatic8>(4, 2);
    test_invariants(combos, sp4, irrep4_0);
    check_size(combos, sp4, irrep4_0, 2);
  }

  SECTION("Combinations<BitsetDynamic>") {
    auto combos = Combinations<BitsetDynamic>(4, 2);
    test_invariants(combos, sp4, irrep4_0);
    check_size(combos, sp4, irrep4_0, 2);
  }

  // =====================================================================
  // LinTable — drop-in Combinations with O(1) index(); same orbits.
  // =====================================================================

  SECTION("LinTable<uint32_t>") {
    auto lt = LinTable<uint32_t>(4, 2);
    test_invariants(lt, sp4, irrep4_0);
    test_invariants(lt, sp4, irrep4_2);
    check_size(lt, sp4, irrep4_0, 2);
  }

  SECTION("LinTable<uint64_t>") {
    auto lt = LinTable<uint64_t>(4, 2);
    test_invariants(lt, sp4, irrep4_0);
    check_size(lt, sp4, irrep4_0, 2);
  }

  // =====================================================================
  // BoundedMultisets
  //
  // With bound=2, BoundedMultisets(4, 2) enumerates the same 16 states
  // as Subsets(4): each slot is 0 or 1. Under cyclic_group(4):
  //   k=0: 6 reps, k=2: 4 reps (Burnside invariant of ordering)
  //
  // BitArrayN: uint64_t backing (1–8 bits per slot)
  // BitArrayLongN: BitsetDynamic backing (1–8 bits per slot)
  // =====================================================================

  SECTION("BoundedMultisets<BitArray1>") {
    auto ms = BoundedMultisets<BitArray1>(4, 2);
    test_invariants(ms, sp4, irrep4_0);
    test_invariants(ms, sp4, irrep4_2);
    check_size(ms, sp4, irrep4_0, 6);
    check_size(ms, sp4, irrep4_2, 4);
  }

  SECTION("BoundedMultisets<BitArray2>") {
    auto ms = BoundedMultisets<BitArray2>(4, 2);
    test_invariants(ms, sp4, irrep4_0);
    check_size(ms, sp4, irrep4_0, 6);
  }

  SECTION("BoundedMultisets<BitArray3>") {
    // bound=4 uses 3 bits; cyclic_group(4) on 4^4=256 states
    auto ms = BoundedMultisets<BitArray3>(4, 4);
    test_invariants(ms, sp4, irrep4_0);
    test_invariants(ms, sp4, irrep4_1);
  }

  SECTION("BoundedMultisets<BitArray4>") {
    auto ms = BoundedMultisets<BitArray4>(4, 2);
    test_invariants(ms, sp4, irrep4_0);
    check_size(ms, sp4, irrep4_0, 6);
  }

  SECTION("BoundedMultisets<BitArray5>") {
    auto ms = BoundedMultisets<BitArray5>(4, 2);
    test_invariants(ms, sp4, irrep4_0);
    check_size(ms, sp4, irrep4_0, 6);
  }

  SECTION("BoundedMultisets<BitArray6>") {
    auto ms = BoundedMultisets<BitArray6>(4, 2);
    test_invariants(ms, sp4, irrep4_0);
    check_size(ms, sp4, irrep4_0, 6);
  }

  SECTION("BoundedMultisets<BitArray7>") {
    auto ms = BoundedMultisets<BitArray7>(4, 2);
    test_invariants(ms, sp4, irrep4_0);
    check_size(ms, sp4, irrep4_0, 6);
  }

  SECTION("BoundedMultisets<BitArray8>") {
    auto ms = BoundedMultisets<BitArray8>(4, 2);
    test_invariants(ms, sp4, irrep4_0);
    check_size(ms, sp4, irrep4_0, 6);
  }

  SECTION("BoundedMultisets<BitArrayLong1>") {
    auto ms = BoundedMultisets<BitArrayLong1>(4, 2);
    test_invariants(ms, sp4, irrep4_0);
    check_size(ms, sp4, irrep4_0, 6);
    check_size(ms, sp4, irrep4_2, 4);
  }

  SECTION("BoundedMultisets<BitArrayLong2>") {
    auto ms = BoundedMultisets<BitArrayLong2>(4, 2);
    test_invariants(ms, sp4, irrep4_0);
    check_size(ms, sp4, irrep4_0, 6);
  }

  SECTION("BoundedMultisets<BitArrayLong3>") {
    auto ms = BoundedMultisets<BitArrayLong3>(4, 4);
    test_invariants(ms, sp4, irrep4_0);
    test_invariants(ms, sp4, irrep4_1);
  }

  SECTION("BoundedMultisets<BitArrayLong4>") {
    auto ms = BoundedMultisets<BitArrayLong4>(4, 2);
    test_invariants(ms, sp4, irrep4_0);
    check_size(ms, sp4, irrep4_0, 6);
  }

  SECTION("BoundedMultisets<BitArrayLong5>") {
    auto ms = BoundedMultisets<BitArrayLong5>(4, 2);
    test_invariants(ms, sp4, irrep4_0);
    check_size(ms, sp4, irrep4_0, 6);
  }

  SECTION("BoundedMultisets<BitArrayLong6>") {
    auto ms = BoundedMultisets<BitArrayLong6>(4, 2);
    test_invariants(ms, sp4, irrep4_0);
    check_size(ms, sp4, irrep4_0, 6);
  }

  SECTION("BoundedMultisets<BitArrayLong7>") {
    auto ms = BoundedMultisets<BitArrayLong7>(4, 2);
    test_invariants(ms, sp4, irrep4_0);
    check_size(ms, sp4, irrep4_0, 6);
  }

  SECTION("BoundedMultisets<BitArrayLong8>") {
    auto ms = BoundedMultisets<BitArrayLong8>(4, 2);
    test_invariants(ms, sp4, irrep4_0);
    check_size(ms, sp4, irrep4_0, 6);
  }

  // =====================================================================
  // BoundedPartitions
  //
  // cyclic_group(4), n=4, total=2:
  //   BitArray1 (bound=2): C(4,2)=6 states, same orbits as Combinations(4,2)
  //     → 2 reps for k=0
  //   BitArray2..8, BitArrayLong2..8 (bound=3): 10 states, 3 orbits:
  //     {(2,0,0,0) orbit}, {(1,1,0,0) orbit}, {(1,0,1,0) orbit}
  //     → 3 reps for k=0
  // =====================================================================

  SECTION("BoundedPartitions<BitArray1>") {
    // bound=2: sequences of 0s and 1s summing to 2 — same as Combinations(4,2)
    auto bp = BoundedPartitions<BitArray1>(4, 2, 2);
    test_invariants(bp, sp4, irrep4_0);
    check_size(bp, sp4, irrep4_0, 2);
  }

  SECTION("BoundedPartitions<BitArray2>") {
    auto bp = BoundedPartitions<BitArray2>(4, 2, 3);
    test_invariants(bp, sp4, irrep4_0);
    test_invariants(bp, sp4, irrep4_2);
    check_size(bp, sp4, irrep4_0, 3);
    check_size(bp, sp4, irrep4_2, 3);
  }

  SECTION("BoundedPartitions<BitArray3>") {
    auto bp = BoundedPartitions<BitArray3>(4, 2, 3);
    test_invariants(bp, sp4, irrep4_0);
    check_size(bp, sp4, irrep4_0, 3);
  }

  SECTION("BoundedPartitions<BitArray4>") {
    auto bp = BoundedPartitions<BitArray4>(4, 2, 3);
    test_invariants(bp, sp4, irrep4_0);
    check_size(bp, sp4, irrep4_0, 3);
  }

  SECTION("BoundedPartitions<BitArray5>") {
    auto bp = BoundedPartitions<BitArray5>(4, 2, 3);
    test_invariants(bp, sp4, irrep4_0);
    check_size(bp, sp4, irrep4_0, 3);
  }

  SECTION("BoundedPartitions<BitArray6>") {
    auto bp = BoundedPartitions<BitArray6>(4, 2, 3);
    test_invariants(bp, sp4, irrep4_0);
    check_size(bp, sp4, irrep4_0, 3);
  }

  SECTION("BoundedPartitions<BitArray7>") {
    auto bp = BoundedPartitions<BitArray7>(4, 2, 3);
    test_invariants(bp, sp4, irrep4_0);
    check_size(bp, sp4, irrep4_0, 3);
  }

  SECTION("BoundedPartitions<BitArray8>") {
    auto bp = BoundedPartitions<BitArray8>(4, 2, 3);
    test_invariants(bp, sp4, irrep4_0);
    check_size(bp, sp4, irrep4_0, 3);
  }

  SECTION("BoundedPartitions<BitArrayLong1>") {
    auto bp = BoundedPartitions<BitArrayLong1>(4, 2, 2);
    test_invariants(bp, sp4, irrep4_0);
    check_size(bp, sp4, irrep4_0, 2);
  }

  SECTION("BoundedPartitions<BitArrayLong2>") {
    auto bp = BoundedPartitions<BitArrayLong2>(4, 2, 3);
    test_invariants(bp, sp4, irrep4_0);
    check_size(bp, sp4, irrep4_0, 3);
  }

  SECTION("BoundedPartitions<BitArrayLong3>") {
    auto bp = BoundedPartitions<BitArrayLong3>(4, 2, 3);
    test_invariants(bp, sp4, irrep4_0);
    check_size(bp, sp4, irrep4_0, 3);
  }

  SECTION("BoundedPartitions<BitArrayLong4>") {
    auto bp = BoundedPartitions<BitArrayLong4>(4, 2, 3);
    test_invariants(bp, sp4, irrep4_0);
    check_size(bp, sp4, irrep4_0, 3);
  }

  SECTION("BoundedPartitions<BitArrayLong5>") {
    auto bp = BoundedPartitions<BitArrayLong5>(4, 2, 3);
    test_invariants(bp, sp4, irrep4_0);
    check_size(bp, sp4, irrep4_0, 3);
  }

  SECTION("BoundedPartitions<BitArrayLong6>") {
    auto bp = BoundedPartitions<BitArrayLong6>(4, 2, 3);
    test_invariants(bp, sp4, irrep4_0);
    check_size(bp, sp4, irrep4_0, 3);
  }

  SECTION("BoundedPartitions<BitArrayLong7>") {
    auto bp = BoundedPartitions<BitArrayLong7>(4, 2, 3);
    test_invariants(bp, sp4, irrep4_0);
    check_size(bp, sp4, irrep4_0, 3);
  }

  SECTION("BoundedPartitions<BitArrayLong8>") {
    auto bp = BoundedPartitions<BitArrayLong8>(4, 2, 3);
    test_invariants(bp, sp4, irrep4_0);
    check_size(bp, sp4, irrep4_0, 3);
  }

  // =====================================================================
  // SchaeferTable — same enumeration as BoundedPartitions, faster index()
  // Only instantiated for BitArray1..8 (not BitArrayLong)
  // =====================================================================

  SECTION("SchaeferTable<BitArray1>") {
    auto st = SchaeferTable<BitArray1>(4, 2, 2);
    test_invariants(st, sp4, irrep4_0);
    check_size(st, sp4, irrep4_0, 2);
  }

  SECTION("SchaeferTable<BitArray2>") {
    auto st = SchaeferTable<BitArray2>(4, 2, 3);
    test_invariants(st, sp4, irrep4_0);
    test_invariants(st, sp4, irrep4_2);
    check_size(st, sp4, irrep4_0, 3);
    check_size(st, sp4, irrep4_2, 3);
  }

  SECTION("SchaeferTable<BitArray3>") {
    auto st = SchaeferTable<BitArray3>(4, 2, 3);
    test_invariants(st, sp4, irrep4_0);
    check_size(st, sp4, irrep4_0, 3);
  }

  SECTION("SchaeferTable<BitArray4>") {
    auto st = SchaeferTable<BitArray4>(4, 2, 3);
    test_invariants(st, sp4, irrep4_0);
    check_size(st, sp4, irrep4_0, 3);
  }

  SECTION("SchaeferTable<BitArray5>") {
    auto st = SchaeferTable<BitArray5>(4, 2, 3);
    test_invariants(st, sp4, irrep4_0);
    check_size(st, sp4, irrep4_0, 3);
  }

  SECTION("SchaeferTable<BitArray6>") {
    auto st = SchaeferTable<BitArray6>(4, 2, 3);
    test_invariants(st, sp4, irrep4_0);
    check_size(st, sp4, irrep4_0, 3);
  }

  SECTION("SchaeferTable<BitArray7>") {
    auto st = SchaeferTable<BitArray7>(4, 2, 3);
    test_invariants(st, sp4, irrep4_0);
    check_size(st, sp4, irrep4_0, 3);
  }

  SECTION("SchaeferTable<BitArray8>") {
    auto st = SchaeferTable<BitArray8>(4, 2, 3);
    test_invariants(st, sp4, irrep4_0);
    check_size(st, sp4, irrep4_0, 3);
  }

  // =====================================================================
  // Cross-check using cyclic_group(3): verify rep counts via Burnside
  //
  // cyclic_group(3) on BoundedMultisets(3, bound=2):
  //   8 states, under rotation: {(0,0,0)}, {(1,0,0),(0,1,0),(0,0,1)},
  //                              {(1,1,0),(0,1,1),(1,0,1)}, {(1,1,1)}
  //   k=0: 4 reps
  //   k=1: (0,0,0) and (1,1,1) excluded → 2 reps
  // =====================================================================

  SECTION("BoundedMultisets<BitArray1> cyclic_group(3)") {
    auto ms = BoundedMultisets<BitArray1>(3, 2);
    test_invariants(ms, sp3, irrep3_0);
    test_invariants(ms, sp3, irrep3_1);
    check_size(ms, sp3, irrep3_0, 4);
    check_size(ms, sp3, irrep3_1, 2);
  }

  SECTION("BoundedMultisets<BitArrayLong1> cyclic_group(3)") {
    auto ms = BoundedMultisets<BitArrayLong1>(3, 2);
    test_invariants(ms, sp3, irrep3_0);
    test_invariants(ms, sp3, irrep3_1);
    check_size(ms, sp3, irrep3_0, 4);
    check_size(ms, sp3, irrep3_1, 2);
  }

  // =====================================================================
  // Triangular lattice 12-site space group (72 elements, non-abelian)
  //
  // This tests correctness on a physically-realistic group and various
  // types of irreps:
  //   Gamma.C6.A  – real, trivial (all chi=1)
  //   Gamma.C6.B  – real, non-trivial (chi in {+1,-1})
  //   K.C3.Ea     – complex (K-point, involves cube roots of unity)
  //
  // Enumerations used:
  //   Combinations<uint32_t>(12,2)       – C(12,2)=66 states, fast
  //   LinTable<uint32_t>(12,6)           – C(12,6)=924 states, half-filling
  //   BoundedPartitions<BitArray2>(12,2,3) – 78 states (sequences summing to 2)
  //   SchaeferTable<BitArray2>(12,2,3)     – same states, fast index()
  // =====================================================================

  SECTION("triangular_12_site_group") {
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/triangular.j1j2jch/"
                        "triangular.12.j1j2jch.sublattices.fsl.toml";
    auto fl = FileToml(lfile);

    // Load three qualitatively different irreps; each irrep's group may differ
    // (K.C3.Ea uses a 24-element little group, Gamma irreps use all 72).
    // SitePermutation is constructed from irrep.group() so that the group
    // check in RepresentativeTable is satisfied by construction.
    auto irrep_trivial = read_representation(fl, "Gamma.C6.A"); // real, chi=1
    auto irrep_B = read_representation(fl, "Gamma.C6.B");       // real, chi=±1
    auto irrep_Kea = read_representation(fl, "K.C3.Ea");        // complex

    auto sp_trivial = SitePermutation(irrep_trivial.group());
    auto sp_B = SitePermutation(irrep_B.group());
    auto sp_Kea = SitePermutation(irrep_Kea.group());

    // ------------------------------------------------------------------
    // Combinations<uint32_t>(12, 2): 66 2-particle states on 12 sites
    // ------------------------------------------------------------------
    {
      auto combos = Combinations<uint32_t>(12, 2);

      test_invariants(combos, sp_trivial, irrep_trivial);
      test_invariants(combos, sp_B, irrep_B);
      test_invariants(combos, sp_Kea, irrep_Kea);

      // For the trivial irrep all states have nonzero norm, so nreps =
      // number of orbits, which must satisfy 1 <= nreps < 66.
      RepresentativeTable<Combinations<uint32_t>> tbl(combos, sp_trivial,
                                                      irrep_trivial);
      REQUIRE(tbl.size() > 0);
      REQUIRE(tbl.size() < combos.size());
    }

    // ------------------------------------------------------------------
    // LinTable<uint32_t>(12, 6): 924 half-filling states
    // Tests a realistically-sized Hilbert-space sector.
    // ------------------------------------------------------------------
    {
      auto lt = LinTable<uint32_t>(12, 6);
      test_invariants(lt, sp_trivial, irrep_trivial);
      test_invariants(lt, sp_B, irrep_B);
      test_invariants(lt, sp_Kea, irrep_Kea);

      RepresentativeTable<LinTable<uint32_t>> tbl(lt, sp_trivial, irrep_trivial);
      REQUIRE(tbl.size() > 0);
      REQUIRE(tbl.size() < lt.size());
    }

    // ------------------------------------------------------------------
    // BoundedPartitions<BitArray2>(12,2,3): 78 states with integer
    // occupancies 0/1/2 per site summing to 2. Exercises multi-valued
    // site degrees of freedom.
    // ------------------------------------------------------------------
    {
      auto bp = BoundedPartitions<BitArray2>(12, 2, 3);
      test_invariants(bp, sp_trivial, irrep_trivial);
      test_invariants(bp, sp_B, irrep_B);
      test_invariants(bp, sp_Kea, irrep_Kea);
    }

    // ------------------------------------------------------------------
    // SchaeferTable<BitArray2>(12,2,3): same 78 states with fast index()
    // ------------------------------------------------------------------
    {
      auto st = SchaeferTable<BitArray2>(12, 2, 3);
      test_invariants(st, sp_trivial, irrep_trivial);
      test_invariants(st, sp_B, irrep_B);

      // Both BoundedPartitions and SchaeferTable for same params → same reps
      auto bp = BoundedPartitions<BitArray2>(12, 2, 3);
      RepresentativeTable<BoundedPartitions<BitArray2>> tbl_bp(bp, sp_trivial,
                                                               irrep_trivial);
      RepresentativeTable<SchaeferTable<BitArray2>> tbl_st(st, sp_trivial,
                                                           irrep_trivial);
      REQUIRE(tbl_bp.size() == tbl_st.size());
    }
  }

  Log("done");
} catch (xdiag::Error e) {
  xdiag::error_trace(e);
}
