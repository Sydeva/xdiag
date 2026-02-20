#include <xdiag/bits/bitarray.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/combinatorics/bounded_multisets/bounded_multisets.hpp>
#include <xdiag/combinatorics/bounded_partitions/bounded_partitions.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/utils/timing.hpp>

using namespace xdiag;

int main() try {
  using namespace xdiag::bits;
  using namespace xdiag::combinatorics;

  // Combinations
  // int n = 24;
  // int k = 12;
  // int64_t cnt = 0;
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
  // {
  //   int64_t n = 12;
  //   int64_t q = 4;
  //   int64_t cnt = 0;
  //   tic();
  //   using A = BitArray<uint64_t, 2>;
  //   for (auto s : BoundedMultisets<A>(n, q)) {
  //     ++cnt;
  //     // Log("{}", to_string(s, n));
  //   }
  //   toc();
  //   Log("BoundedMultisets: {} {}", cnt, BoundedMultisets<A>(n, q).size());
  // }

  {
    int64_t n = 15;
    int64_t total = 15;
    int64_t q = 4;
    int64_t cnt = 0;
    tic();
    using A = BitArray<uint64_t, 2>;
    for (auto s : BoundedPartitions<A>(n, total, q)) {
      ++cnt;
      // Log("{}", to_string(s, n));
    }
    toc();
    Log("BoundedPartitions: {} {}", cnt,
        BoundedPartitions<A>(n, total, q).size());
  }

} catch (Error e) {
  error_trace(e);
}
