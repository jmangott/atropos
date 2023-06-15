#include <chrono>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>
#include <lr/coefficients.hpp>
#include <lr/lr.hpp>

#include "grid_class.hpp"
#include "index_functions.hpp"
#include "integrators.hpp"
#include "io_functions.hpp"
#include "kl_step_functions.hpp"
#include "parameters.hpp"
#include "partition_class.hpp"
#include "print_functions.hpp"
#include "s_step_functions.hpp"
#include "timer_class.hpp"
#include "weight_functions.hpp"


using std::cout;
using std::endl;
using std::string;
using std::vector;


int main()
{
    // get_time::start("main");
    /////////////////////////////////////////////
    /////////////////// SETUP ///////////////////
    /////////////////////////////////////////////

    grid_info grid(kM1, kM2, kR, kN1, kN2, kBinsize1, kBinsize2, kLiml1, kLiml2);

    // Perform partition of the network
    partition_info<1> partition1(grid, mysystem);
    partition_info<2> partition2(grid, mysystem);

    // Integration weights
    vector<multi_array<double, 2>> w_x_dep(mysystem.mu());

    // Coefficients
    vector<multi_array<double, 3>> c_coeff1(mysystem.mu()), d_coeff1(mysystem.mu());
    vector<multi_array<double, 3>> c_coeff2(mysystem.mu()), d_coeff2(mysystem.mu());
    multi_array<double, 5> e_coeff({mysystem.mu(), grid.r, grid.r, grid.r, grid.r}), f_coeff({mysystem.mu(), grid.r, grid.r, grid.r, grid.r});

    // Allocate memory
    for (Index mu = 0; mu < mysystem.mu(); mu++)
    {
        w_x_dep[mu].resize({partition1.dx_dep(mu), partition2.dx_dep(mu)});
        c_coeff1[mu].resize({partition1.dx_dep(mu), grid.r, grid.r});
        d_coeff1[mu].resize({partition1.dx_dep(mu), grid.r, grid.r});
        c_coeff2[mu].resize({partition2.dx_dep(mu), grid.r, grid.r});
        d_coeff2[mu].resize({partition2.dx_dep(mu), grid.r, grid.r});
    }

    // Calculate the integration weights
    double min_prop, max_prop;
    double norm;
    CalculateWeightDep(w_x_dep, min_prop, max_prop, mysystem, grid, partition1, partition2);

    // Container for the initial values
    multi_array<double, 2> xx1, xx2, ss;
    vector<const double *> x1, x2;
    Index n_basisfunctions;

    // Read initial values and store them in `xx1` and `xx2`
    ReadNC("input/input.nc", xx1, xx2, ss, n_basisfunctions);

    // Objects for setting up X1 and X2 for t = 0
    double *it1 = xx1.begin();
    double *it2 = xx2.begin();
    for (Index i = 0; i < n_basisfunctions; i++)
    {
        x1.push_back(it1);
        x2.push_back(it2);
        it1 += grid.dx1;
        it2 += grid.dx2;
    }

    // Low rank structure (for storing X1, X2 and S)
    lr2<double> lr_sol(grid.r, {grid.dx1, grid.dx2});

    // Inner products
    std::function<double(double *, double *)> ip_xx1;
    std::function<double(double *, double *)> ip_xx2;
    blas_ops blas;

    vector<Index> sigma1(mysystem.mu()), sigma2(mysystem.mu());

    // Set up the low-rank structure and the inner products
    ip_xx1 = inner_product_from_const_weight((double) grid.h1_mult, grid.dx1);
    ip_xx2 = inner_product_from_const_weight((double) grid.h2_mult, grid.dx2);
    initialize(lr_sol, x1, x2, ip_xx1, ip_xx2, blas);

    // TODO: this lines are actually superfluous, provided Ensign works correctly
    if (n_basisfunctions == kR)
    {
        ReadNC("input/input.nc", xx1, xx2, ss, n_basisfunctions);
        lr_sol.S = ss;
    }

    // Calculate the shift amount for all reactions (this has to be done only once)
    CalculateShiftAmount(sigma1, sigma2, mysystem, grid);

    // Check if folder in ../output/ exists, otherwise create folder
    std::stringstream fname;
    fname << "output/" << kFilename;
    std::filesystem::create_directory(fname.str());

    // Store initial values
    fname.str("");
    fname << "output/" << kFilename << "/output_t0.nc";
    double t = 0.0;
    WriteNC(fname.str(), lr_sol, mysystem.species_names, grid, &t, &kTau);

    // Diagnostics
    PrintDiagnostics(grid, min_prop, max_prop, kTau, kSecondOrder, kNSubsteps);

    // Number of time steps
    Index kNsteps = ceil(kTstar / kTau);

    auto start_time(std::chrono::high_resolution_clock::now());

    for (Index ts = 0; ts < kNsteps; ts++)
    {
        if (kTstar - t < kTau)
            kTau = kTstar - t;

        if constexpr (kSecondOrder)
        {
            IntegrateSecondOrder(lr_sol, w_x_dep, c_coeff1, d_coeff1, c_coeff2, d_coeff2, e_coeff, f_coeff, sigma1, sigma2, mysystem, grid, partition1, partition2, ip_xx1, ip_xx2, blas, kTau, kNSubsteps, norm);
        }
        else
        {
            IntegrateFirstOrder(lr_sol, w_x_dep, c_coeff1, d_coeff1, c_coeff2, d_coeff2, e_coeff, f_coeff, sigma1, sigma2, mysystem, grid, partition1, partition2, ip_xx1, ip_xx2, blas, kTau, norm);
        }

        t += kTau;

        // Print progress bar
        PrintProgressBar(ts, kNsteps, start_time, norm);

        // Write snapshot
        if ((ts + 1) % kSnapshot == 0 || (ts + 1) == kNsteps)
        {
            fname.str("");
            fname << "output/" << kFilename << "/output_t" << ts + 1 << ".nc";
            WriteNC(fname.str(), lr_sol, mysystem.species_names, grid, &t, &kTau);
        }
        // get_time::stop("main");
    }
        cout << endl << endl; 
        cout << "TIMER RESULTS" << endl;
        cout << "-------------" << endl;
        cout << get_time::sorted_output();
        // get_time::reset();
        // cout << endl;

    return 0;
}