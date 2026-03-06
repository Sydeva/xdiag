// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include "testcases_spinhalf.hpp"

#include <xdiag/algorithms/sparse_diag.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/matrices/matrix.hpp>
#include <xdiag/matrices/sparse/coo_matrix.hpp>
#include <xdiag/matrices/sparse/csr_matrix.hpp>
#include <xdiag/matrices/sparse/csc_matrix.hpp>
#include <xdiag/operators/logic/algebra.hpp>
#include <xdiag/operators/logic/isapprox.hpp>
#include <xdiag/operators/logic/normal_order.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/xdiag_show.hpp>

TEST_CASE("spinhalf", "[blocks]") try {
  using namespace xdiag;

  // OpSum ops = complex(0.0, 1.0) * (Op("SdotS", {0, 1}) * Op("SdotS", {1, 2})
  // -
  //                                  Op("SdotS", {1, 2}) * Op("SdotS", {0,
  //                                  1}));
  // OpSum o1 = normal_order(ops, operators::spin_algebra(), 3);
  // OpSum o2 = normal_order(Op("ScalarChirality", {0, 1, 2}),
  //                         operators::spin_algebra(), 3);
  // XDIAG_SHOW(o1);
  // XDIAG_SHOW(o2);
  // XDIAG_SHOW(o2 +
  //            3.21 * Op("SdotS", {0, 1}) *
  //                Op("Matrix", 0, arma::mat{{0, 1}, {1, 0}}) +
  //            "U" * Op("HubbardU"));

  // XDIAG_SHOW(isapprox(o1, o2));

  int nsites = 6;
  auto block = Spinhalf(nsites, nsites / 2);
  auto ops = testcases::HBchain(nsites, 1.0);
  Log.set_verbosity(2);

  auto H = matrix(ops, block);
  arma::vec eigs;
  arma::eig_sym(eigs, H);

  double e0 = eigval0(ops, block);
  REQUIRE(isapprox(e0, eigs[0]));

  auto H_coo = coo_matrix(ops, block);
  auto H_coo_dense = to_dense(H_coo);

  arma::vec eigs_coo;
  arma::eig_sym(eigs_coo, H_coo_dense);
  REQUIRE(isapprox(H_coo_dense, H));
  REQUIRE(isapprox(e0, eigs_coo[0]));

  auto H_csr = csr_matrix(ops, block);
  auto H_csr_dense = to_dense(H_csr);

  arma::vec eigs_csr;
  arma::eig_sym(eigs_csr, H_csr_dense);
  REQUIRE(isapprox(H_csr_dense, H));
  REQUIRE(isapprox(e0, eigs_csr[0]));

  auto H_csc = csc_matrix(ops, block);
  auto H_csc_dense = to_dense(H_csc);

  arma::vec eigs_csc;
  arma::eig_sym(eigs_csc, H_csc_dense);
  REQUIRE(isapprox(H_csc_dense, H));
  REQUIRE(isapprox(e0, eigs_csc[0]));

  
  Log("e0: {}, eigs[0]: {}, eigs_coo[0]: {}, eigs_csc[0]: {}, eigs_csr[0]: {}", e0, eigs[0],
      eigs_coo[0], eigs_csc[0], eigs_csr[0]);

} catch (xdiag::Error e) {
  error_trace(e);
}
