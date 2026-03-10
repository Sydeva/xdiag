#include "../catch.hpp"

#include <iostream>

#include <xdiag/symmetries/cyclic_group.hpp>
#include <xdiag/symmetries/representation.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;
using namespace arma;

TEST_CASE("representation", "[symmetries]") try {
  // irrep * conj(irrep) is always real
  for (int64_t n = 3; n < 8; ++n) {
    for (int64_t k = 0; k < n; ++k) {
      auto irrep = cyclic_group_irrep(n, k);
      arma::cx_vec chars_hc = arma::conj(irrep.characters().as<arma::cx_vec>());
      auto irrep_hc = Representation(irrep.group(), arma::cx_vec(chars_hc));
      REQUIRE((irrep * irrep_hc).isreal());
    }
  }

  // Trivial representation: all-ones real characters
  for (int64_t n = 2; n < 7; ++n) {
    auto group = cyclic_group(n);
    auto trivial = Representation(group);
    REQUIRE(trivial.isreal());
    REQUIRE(trivial.size() == n);
    for (auto c : trivial.characters().as<arma::vec>())
      REQUIRE(c == Approx(1.0));
  }

  // isreal distinguishes real (k=0) from complex (0 < k < n/2) irreps
  for (int64_t n = 5; n < 8; ++n) {
    REQUIRE(cyclic_group_irrep(n, 0).isreal());
    REQUIRE_FALSE(cyclic_group_irrep(n, 1).isreal());
  }

  // isapprox: same irrep compares equal, different irreps do not
  for (int64_t n = 3; n < 7; ++n) {
    auto r0 = cyclic_group_irrep(n, 0);
    REQUIRE(isapprox(r0, cyclic_group_irrep(n, 0)));
    REQUIRE_FALSE(isapprox(r0, cyclic_group_irrep(n, 1)));
  }

  // Characters that violate c(g)*c(h)=c(gh) should throw
  {
    auto group = cyclic_group(4);
    REQUIRE_THROWS(Representation(group, std::vector<double>{1.0, 2.0, 1.0, 2.0}));
  }

} catch (xdiag::Error e) {
  xdiag::error_trace(e);
}
