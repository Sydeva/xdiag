//
// Created by ge45hub on 20/04/2026.
//
#include <complex.h>
#include <xdiag/all.hpp>
#include <xdiag/algorithms/lanczos/tmatrix.hpp>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace xdiag;

int main(int argc, char **argv) try {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0]
              << " <gs_symmetry_input.toml> <translation_symmetry_input.toml> [output.txt]\n";
    return 1;
  }

  std::string gs_input = argv[1];
  std::string translation_input = argv[2];
  std::string output_path =
      (argc >= 4) ? argv[3]
                  : "/home/t30/all/ge45hub/CLionProjects/xdiagAFM/cantedAFM/strucfac/teststrucfacwdmian02.outfile.txt";

  int nsites = 24;
  int max_iterations = 800;
  int numbeigs = 16;

  // Ground-state block is defined from the point-group symmetry file.
  std::string gs_irrep_label = "Gamma.C1.A";

  // Translation-only irreps used to construct momentum-resolved S_q operators.
  std::vector<std::string> momentum_irreps = {
    "Gamma.C1.A", "M1.C1.A", "X.C1.A"
  };
    //"X1.C1.A", "X2.C1.A", "Kp.C1.A","M2.C1.A",    "M3.C1.A", "K.C1.A",  "X3.C1.A", "X4.C1.A", "X5.C1.A"};

  std::cout << "Reading input files:\n"
            << "  GS symmetry file: " << gs_input << "\n"
            << "  Translation symmetry file: " << translation_input << "\n";

  auto gs_file = FileToml(gs_input);
  auto tr_file = FileToml(translation_input);
  OpSum ops = read_opsum(gs_file, "Interactions");

  OpSum Sz_ops;
  for (int i=0; i<nsites; ++i) {
    Sz_ops += "B" * Op("Sz", {i});
  }
  Sz_ops["B"] = 0.75;
  ops += Sz_ops;

  ops["J1"] = 1.0;
  ops["D"] = complex(0,0.1);
  ops["Jpm"] = 0.02;
  ops["B"] = 0.75;


  auto gs_irrep = read_representation(tr_file, gs_irrep_label, "Symmetries");
  auto gs_block = Spinhalf(nsites, gs_irrep);

  std::cout << "Computing ground state in irrep " << gs_irrep_label << " ...\n";
  auto [e0, gs] = eig0(ops, gs_block, 1e-12, max_iterations);
  gs.make_complex();
  std::cout << "Ground-state energy: " << std::fixed << std::setprecision(12) << e0
            << "\n\n";

  std::ofstream out(output_path);
  if (!out.is_open()) {
    std::cerr << "Could not open output file for writing: " << output_path << "\n";
    return 1;
  }

  out << "momentum,energy,weight\n";

  for (auto const &irrep_label : momentum_irreps) {
    std::cout << "Processing momentum irrep " << irrep_label << " ...\n";

    auto q_irrep = read_representation(tr_file, irrep_label, "Symmetries");
    OpSum SmSm;
    SmSm += Op("S+S+", {0,1});
    SmSm += Op("S+S+", {0,23});
    SmSm += Op("S+S+", {0,7});
    auto S_q = symmetrize(SmSm, q_irrep);

    State Av = apply(S_q, gs);
    double nrm = norm(Av);
    if (nrm < 1e-14) {
      std::cout << "  skipped (norm below threshold)\n";
      continue;
    }
    Av /= nrm;

    auto res = eigvals_lanczos_inplace(ops, Av, numbeigs, 1e-12, max_iterations, 1e-7);

    std::cout << " rescrit " << res.criterion << ")\n";

    std::vector<double> alphas(res.alphas.begin(), res.alphas.end());
    std::vector<double> betas(res.betas.begin(), res.betas.end());
    Tmatrix tmat(alphas, betas);
    auto [energies_t, evecs_t] = tmat.eigen();

    for (arma::uword i = 0; i < energies_t.n_elem; ++i) {
      double energy = energies_t(i) - e0;
      double amp = std::abs(evecs_t(0, i));
      double weight = (nrm * nrm) * (amp * amp);
      out << irrep_label << "," << std::setprecision(16) << energy << ","
          << weight << "\n";
    }

    std::cout << "  wrote " << energies_t.n_elem << " poles\n";
  }

  out.close();
  std::cout << "\nDone. Wrote long-form rows to: " << output_path << "\n";
  return 0;
} catch (Error e) {
  error_trace(e);
}
