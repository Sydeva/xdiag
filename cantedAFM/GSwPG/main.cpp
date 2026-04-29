//
// Created by ge45hub on 20/04/2026.
//
#include <complex.h>
#include <xdiag/all.hpp>
#include <iostream>
#include <iomanip>

using namespace xdiag;



int main() try {

    std::vector<double> energies;
    std::vector<std::string> IrrepList;
    std::vector<int> Szsector;
    int Nsites = 24;
    int numbeig = 12; // number of Eigenvalues to converge.

    /*std::vector<std::string> Irreps = {
        "Gamma.C6.A", "Gamma.C6.B", "Gamma.C6.E1a", "Gamma.C6.E1b", "Gamma.C6.E2a",
        "Gamma.C6.E2b",  "M.C2.A", "M.C2.B", "K.C3.A", "K.C3.Ea", "K.C3.Ea"}; */

    std::vector<std::string> Irreps = {
        "Gamma.C6.A", "Gamma.C6.B", "Gamma.C6.E1a", "Gamma.C6.E1b", "Gamma.C6.E2a",
        "Gamma.C6.E2b"};

    std::cout << "Starting diagonalization for " << Nsites << " sites\n";
    std::cout << "Number of eigenvalues per sector: " << numbeig << "\n";

    auto fl = FileToml("/home/t30/all/ge45hub/CLionProjects/xdiagAFM/cantedAFM/GSwPG/honeycomb.24.no_pairflip.toml");

    OpSum ops = read_opsum(fl, "Interactions");

    for (int i=0; i<Nsites; ++i) {
        ops += "B" * Op("Sz", {i});
    }

    ops["J1"] = 1.0;
    ops["D"] = complex(0,0.001);
    ops["B"] = 0.75;

    std::cout << "Hamiltonian configured. Starting symmetry sector loop...\n\n";

    int irrep_count = 0;
    for (auto irrep : Irreps) {
        irrep_count++;
        auto irrep2 = read_representation(fl, irrep, "Symmetries");

        std::cout << "  [" << irrep_count << "/" << Irreps.size() << "] Processing irrep: " << irrep << "\n";

        int sz_count = 0;
        for (int nup = Nsites/2-2-6; nup <= Nsites/2-2+6; nup++) {
            sz_count++;
            auto block = Spinhalf(Nsites, nup, irrep2);
            int64_t block_dim = block.dim();

            std::cout << "    Sz sector " << std::setw(2) << sz_count << " (nup=" << std::setw(2) << nup
                      << "): block dimension = " << block_dim << "\n";

            auto res = eigvals_lanczos(ops, block, numbeig);
            arma::vec eig0 = res.eigenvalues;
            std::cout << " rescrit " << res.criterion << ")\n";
            std::cout << "      → converged " << eig0.n_elem << " eigenvalues (E_min = "
                      << std::fixed << std::setprecision(6) << eig0(0) << ")\n";

            for (int i = 0; i < eig0.n_elem; i++) {
                energies.push_back(eig0[i]);
                IrrepList.push_back(irrep);
                Szsector.push_back(nup);
            }
        }
    }

    std::cout << "\nDiagonalization complete. Total eigenvalues collected: " << energies.size() << "\n\n";

    // Construct the filename
    std::string flstring = "HeisenbergZeemanDMSzvar" +
                           std::to_string(Nsites) + ".outfile.txt";
    std::ofstream outfile(flstring);
    for (int i = 0; i < energies.size(); i++) {
        outfile << energies[i] << "," << Szsector[i] << "," << IrrepList[i] << "\n";
    }
    outfile.close();
    return 0;
} catch (Error e) {
    error_trace(e);
}
