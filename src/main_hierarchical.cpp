#include <filesystem>

#include "coeff_class.hpp"
#include "io_functions.hpp"
#include "index_functions.hpp"
#include "integration_parameters.hpp"
#include "integrators.hpp"
#include "matrix.hpp"
#include "tree_class.hpp"

int main()
{
    double t = 0.0;
    blas_ops blas;
    cme_lr_tree tree;
    const Index kNsteps = ceil(kTstar / Tau); // Number of time steps

    tree.Read("input/input.nc");
    cout << tree;
    tree.Orthogonalize(blas);
    double norm = tree.Normalize();
    cout << "Norm: " << norm << endl; 

    // Check if folder in ../output/ exists, otherwise create folder
    std::stringstream fname;
    fname << "output/" << kFilename;
    std::filesystem::create_directory(fname.str());

    // Store initial values
    fname.str("");
    fname << "output/" << kFilename << "/output_t0.nc";
    tree.Write(fname.str(), t, Tau);

    auto start_time(std::chrono::high_resolution_clock::now());

    for (Index ts = 0; ts < kNsteps; ++ts)
    {
        if (kTstar - t < Tau)
            Tau = kTstar - t;

        TTNIntegrator(tree.root, blas, Tau);
        norm = tree.Normalize();

        t += Tau;

        // Print progress bar
        PrintProgressBar(ts, kNsteps, start_time, norm);

        // Write snapshot
        if ((ts + 1) % kSnapshot == 0 || (ts + 1) == kNsteps)
        {
            fname.str("");
            fname << "output/" << kFilename << "/output_t" << ts + 1 << ".nc";
            tree.Write(fname.str(), t, Tau);
        }
    }

    return 0;
}