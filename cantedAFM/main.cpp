//
// Created by ge45hub on 20/04/2026.
//
#include <xdiag/all.hpp>

using namespace xdiag;

int main() try {

    std::vector<double> energies;
    std::vector<std::string> IrrepList;
    std::vector<int> Szsector;
    int Nsites = 24;
    int numbeig = 12; // number of Eigenvalues to converge.

    std::vector<std::string> Irreps = {
        "Gamma.C6.A", "Gamma.C6.B", "Gamma.C6.E1a", "Gamma.C6.E1b", "Gamma.C6.E2a",
        "Gamma.C6.E2b",  "M.C2.A", "M.C2.B", "K.C3.A", "K.C3.Ea", "K.C3.Ea"};

    auto fl = FileToml("honeycomb.24.toml");

    OpSum ops = read_opsum(fl, "Interactions");

    for (int i=0; i<Nsites; ++i) {
        ops += "B" * Op("Sz", {i});
    }

    ops["J1"] = 1.0;
    ops["B"] = 0.75;


    for (auto irrep : Irreps) {

        auto irrep2 = read_representation(fl, irrep, "Symmetries");

        for (int nup = Nsites/2-2-6; nup <= Nsites/2-2+6; nup++) {
            auto block = Spinhalf(Nsites, nup, irrep2);

            auto res = eigvals_lanczos(ops, block, numbeig);
            arma::vec eig0 = res.eigenvalues;

            for (int i = 0; i < eig0.n_elem; i++) {
                energies.push_back(eig0[i]);
                IrrepList.push_back(irrep);
                Szsector.push_back(nup);
            }
        }
    }

    // Construct the filename
    std::string flstring = "HeisenbergZeemanSzvar" +
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
