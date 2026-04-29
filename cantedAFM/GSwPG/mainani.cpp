//
// Created by ge45hub on 20/04/2026.
//
#include <complex.h>
#include <xdiag/all.hpp>
#include <xdiag/algorithms/lanczos/tmatrix.hpp>
#include <iostream>
#include <iomanip>

using namespace xdiag;



int main() try {

    std::vector<double> energies;
    std::vector<double> SzExpect;

    std::vector<std::string> IrrepList;
    int Nsites = 24;
    int numbeig = 16; // number of Eigenvalues to converge.

    /*std::vector<std::string> Irreps = {
        "Gamma.C6.A", "Gamma.C6.B", "Gamma.C6.E1a", "Gamma.C6.E1b", "Gamma.C6.E2a",
        "Gamma.C6.E2b",  "M.C2.A", "M.C2.B", "K.C3.A", "K.C3.Ea", "K.C3.Ea"}; */

    std::vector<std::string> Irreps = {
        "Gamma.C6.A", "Gamma.C6.B", "Gamma.C6.E1a", "Gamma.C6.E1b"
    };
        //, "Gamma.C6.E1a", "Gamma.C6.E1b", "Gamma.C6.E2a",
        //"Gamma.C6.E2b"};

    std::cout << "Starting diagonalization for " << Nsites << " sites\n";
    std::cout << "Number of eigenvalues per sector: " << numbeig << "\n";

    auto fl = FileToml("/home/t30/all/ge45hub/CLionProjects/xdiagAFM/cantedAFM/strucfac/2a4latwithanianddmi.toml");

    OpSum ops = read_opsum(fl, "Interactions");

    OpSum Sz_ops;
    for (int i=0; i<Nsites; ++i) {
        Sz_ops += "B" * Op("Sz", {i});
    }
    Sz_ops["B"] = 0.75;
    ops += Sz_ops;

    ops["J1"] = 1.0;
    ops["D"] = complex(0,0.1);
    ops["Jpm"] = 0.02;

    std::cout << "Hamiltonian configured. Starting symmetry sector loop...\n\n";

    int irrep_count = 0;
    for (auto irrep : Irreps) {
        irrep_count++;
        auto irrep2 = read_representation(fl, irrep, "Symmetries");

        std::cout << "  [" << irrep_count << "/" << Irreps.size() << "] Processing irrep: " << irrep << "\n";


        auto block = Spinhalf(Nsites, irrep2);
        int64_t block_dim = block.dim();

        std::cout << "block dimension = " << block_dim << "\n";

        auto res = eigs_lanczos(ops, block, numbeig, 1e-12,1000,1e-7);
        arma::vec eig0 = res.eigenvalues;

        std::cout << " rescrit " << res.criterion << ")\n";

        std::cout << "      → converged " << eig0.n_elem << " eigenvalues (E_min = "
                  << std::fixed << std::setprecision(6) << eig0(0) << ")\n";

        for (int i = 0; i < numbeig; i++) {
            energies.push_back(eig0[i]);
            IrrepList.push_back(irrep);
            auto szc = innerC(Sz_ops, res.eigenvectors.col(i));
            if (std::abs(std::imag(szc)) > 1e-10) {
                std::cerr << "Warning: <Sz> has non-negligible imag part = "
                          << std::imag(szc) << "\n";
            }
            SzExpect.push_back(1.25 * std::real(szc));


        }
    }

    std::cout << "\nDiagonalization complete. Total eigenvalues collected: " << energies.size() << "\n\n";

    // Construct the filename
    std::string flstring = "/home/t30/all/ge45hub/CLionProjects/xdiagAFM/cantedAFM/strucfac/HeisenbergZeemanDMJani02" +
                           std::to_string(Nsites) + ".outfile.txt";
    std::ofstream outfile(flstring);
    for (int i = 0; i < energies.size(); i++) {
        outfile << energies[i] << "," << SzExpect[i] << "," << IrrepList[i] << "\n";
    }
    outfile.close();
    return 0;
} catch (Error e) {
    error_trace(e);
}
