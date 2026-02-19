#include <xdiag/bits/bitset.hpp>
#include <xdiag/bits/bitarray.hpp>
#include <xdiag/combinatorics/bounded_multisets/bounded_multisets.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/utils/timing.hpp>

using namespace xdiag;

int main() try {
  using namespace xdiag::bits;
  using namespace xdiag::combinatorics;

  // Combinations
  int n = 24;
  int k = 12;
  int64_t cnt = 0;
  // tic();
  // for (auto s : Combinations<uint64_t>(n, k)) {
  //   ++cnt;
  //   // Log("{}", to_string(make_bitset(s), n));
  // }
  // Log("Combinations short: {} ", cnt);
  // toc();

  // n = 80;
  // k = 4;
  // cnt = 0;
  // tic();
  // for (auto s : Combinations<Bitset<uint64_t, 2>>(n, k)) {
  //   ++cnt;
  //   // Log("{}", to_string(s, n));
  // }
  // Log("Combinations long:  {} ", cnt);
  // toc();

  n = 10;
  int64_t q = 5;
  cnt = 0;
  tic();
  for (auto s : BoundedMultisets<BitArray<uint64_t, 3>>(n, q)) {
    ++cnt;
    // Log("{}", to_string(s, n));
  }
  toc();
  Log("BoundedMultisets: {} {}", cnt, BoundedMultisets<BitArray<uint64_t, 3>>(n, q).size());

} catch (Error e) {
  error_trace(e);
}
