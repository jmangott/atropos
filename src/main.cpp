#include <chrono>
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
#include "io_functions.hpp"
#include "kl_step_functions.hpp"
#include "parameters.hpp"
#include "partition_class.hpp"
#include "print_functions.hpp"
#include "s_step_functions.hpp"
#include "timer_class.hpp"
#include "weight_functions.hpp"

#include "models/reactions_ts.hpp"
// #include "models/reactions_lp.hpp"
// #include "models/reactions_tgfb6.hpp"

using std::cout;
using std::endl;
using std::string;
using std::vector;


int main()
{
    get_time::start("main");
    /////////////////////////////////////////////
    /////////////////// SETUP ///////////////////
    /////////////////////////////////////////////

    grid_info grid(kM1, kM2, kR, kN1, kN2, kK1, kK2, kLiml1, kLiml2);

    // Perform partition of the network
    partition_info<1> partition1(grid, mysystem);
    partition_info<2> partition2(grid, mysystem);

    // Integration weights
    vector<multi_array<double, 2>> w_x_dep(mysystem.mu());

    // Calculate the integration weights
    CalculateWeightDep(w_x_dep, mysystem, grid, partition1, partition2);

    // Declare LR-specific objects
    // Temporary objects for multiplication
    multi_array<double, 2> tmp_x1({grid.dx1, grid.r});
    multi_array<double, 2> tmp_x2({grid.dx2, grid.r});

    // Container for the initial values
    multi_array<double, 2> xx1, xx2, ss;
    vector<const double *> x1, x2;
    double *it1, *it2;

    // Objects for setting up X1 and X2 for t = 0
    if (kFullMatrixInitialCondition)
    {
        xx1.resize({grid.dx1, grid.r});
        xx2.resize({grid.dx2, grid.r});
        it1 = xx1.begin();
        it2 = xx2.begin();
        ss.resize({grid.r, grid.r});
        ReadInMultiArray(ss, "../input/s_input.csv");
        for (int i = 0; i < kR; i++)
        {
            x1.push_back(it1);
            x2.push_back(it2);
            it1 += grid.dx1;
            it2 += grid.dx2;
        }
    }
    else
    {
        xx1.resize({grid.dx1, kNBasisFunctions});
        xx2.resize({grid.dx2, kNBasisFunctions});
        it1 = xx1.begin();
        it2 = xx2.begin();
        x1.push_back(it1);
        x2.push_back(it2);
    }

    ReadInMultiArray(xx1, "../input/x1_input.csv");
    ReadInMultiArray(xx2, "../input/x2_input.csv");

    // Low rank structure (for storing X1, X2 and S)
    lr2<double> lr_sol(grid.r, {grid.dx1, grid.dx2});

    // Inner products
    std::function<double(double *, double *)> ip_xx1;
    std::function<double(double *, double *)> ip_xx2;
    blas_ops blas;
    gram_schmidt gs(&blas);

    vector<Index> sigma1(mysystem.mu()), sigma2(mysystem.mu());

    // Set up the low-rank structure and the inner products
    ip_xx1 = inner_product_from_const_weight(grid.h1_mult, grid.dx1);
    ip_xx2 = inner_product_from_const_weight(grid.h2_mult, grid.dx2);
    initialize(lr_sol, x1, x2, ip_xx1, ip_xx2, blas);

    // TODO: this line is actually superfluous, provided Ensign works correctly
    if (kFullMatrixInitialCondition) lr_sol.S = ss;

    // Calculate the shift amount for all reactions (this has to be done only once)
    CalculateShiftAmount(sigma1, sigma2, mysystem, grid);

    // Write output files for initial values
    WriteOutMultiArray(lr_sol.X, "../output/toggle_switch_new_250_400/x1_output_t0");
    WriteOutMultiArray(lr_sol.S, "../output/toggle_switch_new_250_400/s_output_t0");
    WriteOutMultiArray(lr_sol.V, "../output/toggle_switch_new_250_400/x2_output_t0");

    auto start_time(std::chrono::high_resolution_clock::now());

    double t = 0.0;
    // int t_int;
    for (Index ts = 0; ts < kNsteps; ts++) 
    {
        if (kTstar - t < kTau)
            kTau = kTstar - t;

        /////////////////////////////////////////////
        ////////////////// K-STEP ///////////////////
        /////////////////////////////////////////////

        tmp_x1 = lr_sol.X;
        blas.matmul(tmp_x1, lr_sol.S, lr_sol.X); // lr_sol.X contains now K
        get_time::start("kstep");
        PerformKLStep<1>(sigma1, sigma2, lr_sol, blas, mysystem, grid, partition1, partition2, w_x_dep, kTau);
        get_time::stop("kstep");
        
        // Perform the QR decomposition K = X * S
        gs(lr_sol.X, lr_sol.S, ip_xx1);

        // Calculate sum of X2 along axis 1
        multi_array<double, 1> x2_bar({grid.r});
        for (Index i = 0; i < grid.r; i++)
            x2_bar(i) = 0.0;
        for (Index i = 0; i < grid.dx2; i++)
        {
            for (Index j = 0; j < grid.r; j++)
            {
                x2_bar(j) += lr_sol.V(i, j);
            }
        }

        /////////////////////////////////////////////
        ////////////////// S-STEP ///////////////////
        /////////////////////////////////////////////
        get_time::start("sstep");
        PerformSStep(sigma1, sigma2, lr_sol, blas, mysystem, grid, partition1, partition2, w_x_dep, kTau);
        get_time::stop("sstep");

        /////////////////////////////////////////////
        ////////////////// L-STEP ///////////////////
        /////////////////////////////////////////////

        tmp_x2 = lr_sol.V;
        blas.matmul_transb(tmp_x2, lr_sol.S, lr_sol.V); // lr_sol.V contains now L
        get_time::start("lstep");
        PerformKLStep<2>(sigma1, sigma2, lr_sol, blas, mysystem, grid, partition2, partition1, w_x_dep, kTau);
        get_time::stop("lstep");
        // Perform the QR decomposition L = V * S^T
        gs(lr_sol.V, lr_sol.S, ip_xx2);
        transpose_inplace(lr_sol.S);

        // Calculate sum of X1 along axis 1
        multi_array<double, 1> x1_bar({grid.r});
        for (Index i = 0; i < grid.r; i++)
            x1_bar(i) = 0.0;
        for (Index i = 0; i < grid.dx1; i++)
        {
            for (Index j = 0; j < grid.r; j++)
            {
                x1_bar(j) += lr_sol.X(i, j);
            }
        }

        // Calculate normalization constant
        double norm = 0.0;
        for (Index i = 0; i < grid.r; i++)
        {
            for (Index j = 0; j < grid.r; j++)
            {
                norm += x1_bar(i) * lr_sol.S(i, j) * x2_bar(j);
            }
        }
        
        // Renormalize S
        lr_sol.S /= norm;

        t += kTau;

        // Print progress bar
        PrintProgressBar(ts, kNsteps, start_time);

        std::stringstream fname_x1_output;
        std::stringstream fname_s_output;
        std::stringstream fname_x2_output;

        if ((ts + 1) % kSnapshot == 0)
        {
            // Write snapshot
            // t_int = (int) (ts + 1) * kTau;
            fname_x1_output << "../output/toggle_switch_new_250_400/x1_output_t" << ts + 1;
            fname_s_output << "../output/toggle_switch_new_250_400/s_output_t" << ts + 1;
            fname_x2_output << "../output/toggle_switch_new_250_400/x2_output_t" << ts + 1;
            WriteOutMultiArray(lr_sol.X, fname_x1_output.str());
            WriteOutMultiArray(lr_sol.S, fname_s_output.str());
            WriteOutMultiArray(lr_sol.V, fname_x2_output.str());
        }
    }

    get_time::stop("main");
    cout << endl << endl << "Timer results: " << endl << get_time::sorted_output();

    return 0;
}