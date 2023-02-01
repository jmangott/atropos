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
// #include "reactions_ts.hpp"
#include "reactions_lp.hpp"
#include "s_step_functions.hpp"
#include "weight_functions.hpp"

using std::cout;
using std::endl;
using std::string;
using std::vector;


int main()
{
    /////////////////////////////////////////////
    /////////////////// SETUP ///////////////////
    /////////////////////////////////////////////

    grid_info grid(kM1, kM2, kR, kN1, kN2, kK1, kK2);

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

    // Objects for setting up X1 and X2 for t = 0
    multi_array<double, 2> xx1({grid.dx1, grid.r});
    multi_array<double, 2> xx2({grid.dx2, grid.r});
    vector<const double *> x1, x2;

    // Low rank structure (for storing X1, X2 and S)
    lr2<double> lr_sol(grid.r, {grid.dx1, grid.dx2});

    // Inner products
    std::function<double(double *, double *)> ip_xx1;
    std::function<double(double *, double *)> ip_xx2;
    blas_ops blas;
    gram_schmidt gs(&blas);

    vector<Index> sigma1, sigma2;

    ReadInMultiArray(xx1, "../input/x1_input.csv");
    ReadInMultiArray(xx2, "../input/x2_input.csv");

    // Point to the beginning of every column of X1 and X2
    double *it1 = xx1.begin();
    double *it2 = xx2.begin();
    for (int i = 0; i < kR; i++)
    {
        x1.push_back(it1);
        x2.push_back(it2);
        it1 += grid.dx1;
        it2 += grid.dx2;
    }

    // Set up the low-rank structure and the inner products
    double h_xx_mult1 = 1, h_xx_mult2 = 1;
    for (Index i = 0; i < grid.m1; i++)
        h_xx_mult1 *= grid.h1(i);
    for (Index i = 0; i < grid.m2; i++)
        h_xx_mult2 *= grid.h2(i);

    ip_xx1 = inner_product_from_const_weight(h_xx_mult1, grid.dx1);
    ip_xx2 = inner_product_from_const_weight(h_xx_mult2, grid.dx2);
    initialize(lr_sol, x1, x2, ip_xx1, ip_xx2, blas);

    // For testing
    // WriteOutMultiArray(lr_sol.S, "../output/s_output_init");
    // TODO: these lines are actually superfluous, provided Ensign works correctly
    multi_array<double, 2> ss({grid.r, grid.r});
    ReadInMultiArray(ss, "../input/s_input.csv");
    lr_sol.S = ss;

    // Calculate the shift amount for all reactions (this has to be done only once)
    CalculateShiftAmount(sigma1, sigma2, mysystem, grid);

    // Write output files for initial values
    WriteOutMultiArray(lr_sol.X, "../output/profiling/x1_output_t0");
    WriteOutMultiArray(lr_sol.S, "../output/profiling/s_output_t0");
    WriteOutMultiArray(lr_sol.V, "../output/profiling/x2_output_t0");

    auto start_time(std::chrono::high_resolution_clock::now());

    double t = 0.0;
    int t_int;
    for (Index ts = 0; ts < kNsteps; ts++) 
    {
        if (kTstar - t < kTau)
            kTau = kTstar - t;

        /////////////////////////////////////////////
        ////////////////// K-STEP ///////////////////
        /////////////////////////////////////////////

        tmp_x1 = lr_sol.X;
        blas.matmul(tmp_x1, lr_sol.S, lr_sol.X); // lr_sol.X contains now K
        PerformKLStep<1>(sigma1, sigma2, lr_sol, blas, mysystem, grid, partition1, partition2, w_x_dep, kTau);
        
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

        // auto start_time_s(std::chrono::high_resolution_clock::now());
        PerformSStep(sigma1, sigma2, lr_sol, blas, mysystem, grid, partition1, partition2, w_x_dep, kTau);
        // auto end_time_s(std::chrono::high_resolution_clock::now());
        // auto duration_s_incr = end_time_s - start_time_s;
        // double duration_s = duration_s_incr.count();
        // cout << "S-step: " << duration_s << endl;

        /////////////////////////////////////////////
        ////////////////// L-STEP ///////////////////
        /////////////////////////////////////////////

        tmp_x2 = lr_sol.V;
        blas.matmul_transb(tmp_x2, lr_sol.S, lr_sol.V); // lr_sol.V contains now L
        PerformKLStep<2>(sigma1, sigma2, lr_sol, blas, mysystem, grid, partition2, partition1, w_x_dep, kTau);
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
        // auto end_time(std::chrono::high_resolution_clock::now());
        // auto duration_main_incr = end_time - start_time;
        // double duration_main = duration_main_incr.count();
        // cout << "main: " << duration_main << endl;
        PrintProgressBar(ts, kNsteps, start_time);

        std::stringstream fname_x1_output;
        std::stringstream fname_s_output;
        std::stringstream fname_x2_output;

        if ((ts + 1) % kSnapshot == 0)
        {
            // Write snapshot
            t_int = (int) (ts + 1) * kTau;
            fname_x1_output << "../output/profiling/x1_output_t" << t_int;
            fname_s_output << "../output/profiling/s_output_t" << t_int;
            fname_x2_output << "../output/profiling/x2_output_t" << t_int;
            WriteOutMultiArray(lr_sol.X, fname_x1_output.str());
            WriteOutMultiArray(lr_sol.S, fname_s_output.str());
            WriteOutMultiArray(lr_sol.V, fname_x2_output.str());
        }
    }

    return 0;
}