#include "integrators.hpp"

using std::vector;

double CalculateNorm(lr2<double> &lr_sol, grid_info &grid)
{
    double norm = 0.0;
    multi_array<double, 1> x1_bar({grid.r});
    multi_array<double, 1> x2_bar({grid.r});

    for (Index i = 0; i < grid.r; i++)
    {
        x1_bar(i) = 0.0;
        x2_bar(i) = 0.0;
    }

    for (Index j = 0; j < grid.r; j++)
    {
        for (Index i = 0; i < grid.dx1; i++)
        {
            x1_bar(j) += lr_sol.X(i, j);
        }
        for (Index i = 0; i < grid.dx2; i++)
        {
            x2_bar(j) += lr_sol.V(i, j);
        }
    }

    for (Index i = 0; i < grid.r; i++)
    {
        for (Index j = 0; j < grid.r; j++)
        {
            norm += x1_bar(i) * lr_sol.S(i, j) * x2_bar(j);
        }
    }

    return norm;
}


void IntegrateFirstOrder(lr2<double> &lr_sol, const vector<multi_array<double, 2>> &w_x_dep, vector<multi_array<double, 3>> &c_coeff1, vector<multi_array<double, 3>> &d_coeff1, vector<multi_array<double, 3>> &c_coeff2, vector<multi_array<double, 3>> &d_coeff2, multi_array<double, 5> &e_coeff, multi_array<double, 5> &f_coeff, const vector<Index> sigma1, const vector<Index> sigma2, mysys &mysystem, grid_info &grid, partition_info<1> &partition1, partition_info<2> &partition2, std::function<double(double *, double *)> ip_xx1, std::function<double(double *, double *)> ip_xx2, blas_ops &blas, double tau)
{
    gram_schmidt gs(&blas);
    double norm;

    // Temporary objects for multiplication and integration
    multi_array<double, 2> tmp_x1({grid.dx1, grid.r});
    multi_array<double, 2> tmp_x2({grid.dx2, grid.r});
    multi_array<double, 2> tmp_s({grid.r, grid.r});

    /////////////////////////////////////////////
    ////////////////// K-STEP ///////////////////
    /////////////////////////////////////////////

    get_time::start("kstep");
    get_time::start("kstep_coeff");
    CalculateCoefficientsKL<1>(c_coeff1, d_coeff1, sigma2, lr_sol, blas, mysystem, grid, partition1, partition2, w_x_dep);
    get_time::stop("kstep_coeff");
    tmp_x1 = lr_sol.X;
    blas.matmul(tmp_x1, lr_sol.S, lr_sol.X); // lr_sol.X contains now K
    PerformKLStep<1>(tmp_x1, lr_sol.X, c_coeff1, d_coeff1, sigma1, blas, mysystem, grid, partition1, w_x_dep, tau);
    lr_sol.X += tmp_x1;
    get_time::stop("kstep");

    // Perform the QR decomposition K = X1 * S
    gs(lr_sol.X, lr_sol.S, ip_xx1);

    /////////////////////////////////////////////
    ////////////////// S-STEP ///////////////////
    /////////////////////////////////////////////

    get_time::start("sstep");
    get_time::start("sstep_coeff");
    CalculateCoefficientsS(e_coeff, f_coeff, sigma1, sigma2, lr_sol, blas, mysystem, grid, partition1, partition2, w_x_dep);
    get_time::stop("sstep_coeff");
    PerformSStep(tmp_s, lr_sol.S, e_coeff, f_coeff, sigma1, sigma2, blas, mysystem, grid, partition1, partition2, w_x_dep, tau);
    lr_sol.S -= tmp_s;
    get_time::stop("sstep");

    /////////////////////////////////////////////
    ////////////////// L-STEP ///////////////////
    /////////////////////////////////////////////

    get_time::start("lstep");
    get_time::start("lstep_coeff");
    CalculateCoefficientsKL<2>(c_coeff2, d_coeff2, sigma1, lr_sol, blas, mysystem, grid, partition2, partition1, w_x_dep);
    get_time::stop("lstep_coeff");
    tmp_x2 = lr_sol.V;
    blas.matmul_transb(tmp_x2, lr_sol.S, lr_sol.V); // lr_sol.V contains now L
    PerformKLStep<2>(tmp_x2, lr_sol.V, c_coeff2, d_coeff2, sigma2, blas, mysystem, grid, partition2, w_x_dep, tau);
    lr_sol.V += tmp_x2;
    get_time::stop("lstep");

    // Perform the QR decomposition L = X2 * S^T
    gs(lr_sol.V, lr_sol.S, ip_xx2);
    transpose_inplace(lr_sol.S);

    /////////////////////////////////////////////
    /////////////// NORMALIZATION ///////////////
    /////////////////////////////////////////////

    // Renormalize S
    norm = CalculateNorm(lr_sol, grid);
    lr_sol.S /= norm;
    cout << "norm(P): " << norm << endl;
}


///////////////////////////////////////////////////////
/////////////////// 100 Euler steps ///////////////////
///////////////////////////////////////////////////////

void IntegrateSecondOrder(lr2<double> &lr_sol, const vector<multi_array<double, 2>> &w_x_dep, vector<multi_array<double, 3>> &c_coeff1, vector<multi_array<double, 3>> &d_coeff1, vector<multi_array<double, 3>> &c_coeff2, vector<multi_array<double, 3>> &d_coeff2, multi_array<double, 5> &e_coeff, multi_array<double, 5> &f_coeff, const vector<Index> sigma1, const vector<Index> sigma2, mysys &mysystem, grid_info &grid, partition_info<1> &partition1, partition_info<2> &partition2, std::function<double(double *, double *)> ip_xx1, std::function<double(double *, double *)> ip_xx2, blas_ops &blas, double tau)
{
    gram_schmidt gs(&blas);
    double norm;

    double tau_sub = 0.01;
    Index n_substeps = ceil(1.0 / tau_sub);

    // Temporary objects for multiplication and integration
    multi_array<double, 2> tmp_x1_0({grid.dx1, grid.r});
    multi_array<double, 2> tmp_x2_0({grid.dx2, grid.r});
    multi_array<double, 2> tmp_s_0({grid.r, grid.r});

    /////////////////////////////////////////////
    ///////////////// 1/2 K-STEP ////////////////
    /////////////////////////////////////////////

    get_time::start("kstep");
    get_time::start("kstep_coeff");
    CalculateCoefficientsKL<1>(c_coeff1, d_coeff1, sigma2, lr_sol, blas, mysystem, grid, partition1, partition2, w_x_dep);
    get_time::stop("kstep_coeff");
    tmp_x1_0 = lr_sol.X;
    blas.matmul(tmp_x1_0, lr_sol.S, lr_sol.X); // lr_sol.X contains now K
    for (Index i = 0; i < n_substeps; i++)
    {
        PerformKLStep<1>(tmp_x1_0, lr_sol.X, c_coeff1, d_coeff1, sigma1, blas, mysystem, grid, partition1, w_x_dep, 0.5 * tau * tau_sub);
        lr_sol.X += tmp_x1_0;
    }
    get_time::stop("kstep");

    // Perform the QR decomposition K = X1 * S
    gs(lr_sol.X, lr_sol.S, ip_xx1);

    /////////////////////////////////////////////
    //////////////// 1/2 S-STEP /////////////////
    /////////////////////////////////////////////

    get_time::start("sstep");
    get_time::start("sstep_coeff");
    CalculateCoefficientsS(e_coeff, f_coeff, sigma1, sigma2, lr_sol, blas, mysystem, grid, partition1, partition2, w_x_dep);
    get_time::stop("sstep_coeff");
    for (Index i = 0; i < n_substeps; i++)
    {
        PerformSStep(tmp_s_0, lr_sol.S, e_coeff, f_coeff, sigma1, sigma2, blas, mysystem, grid, partition1, partition2, w_x_dep, 0.5 * tau * tau_sub);
        lr_sol.S -= tmp_s_0;
    }
    get_time::stop("sstep");

    /////////////////////////////////////////////
    ////////////////// L-STEP ///////////////////
    /////////////////////////////////////////////

    get_time::start("lstep");
    get_time::start("lstep_coeff");
    CalculateCoefficientsKL<2>(c_coeff2, d_coeff2, sigma1, lr_sol, blas, mysystem, grid, partition2, partition1, w_x_dep);
    get_time::stop("lstep_coeff");
    tmp_x2_0 = lr_sol.V;
    blas.matmul_transb(tmp_x2_0, lr_sol.S, lr_sol.V); // lr_sol.V contains now L
    for (Index i = 0; i < n_substeps; i++)
    {
        PerformKLStep<2>(tmp_x2_0, lr_sol.V, c_coeff2, d_coeff2, sigma2, blas, mysystem, grid, partition2, w_x_dep, tau * tau_sub);
        lr_sol.V += tmp_x2_0;
    }
    get_time::stop("lstep");

    // Perform the QR decomposition L = X2 * S^T
    gs(lr_sol.V, lr_sol.S, ip_xx2);
    transpose_inplace(lr_sol.S);

    /////////////////////////////////////////////
    //////////////// 1/2 S-STEP /////////////////
    /////////////////////////////////////////////

    get_time::start("sstep");
    get_time::start("sstep_coeff");
    CalculateCoefficientsS(e_coeff, f_coeff, sigma1, sigma2, lr_sol, blas, mysystem, grid, partition1, partition2, w_x_dep);
    get_time::stop("sstep_coeff");
    for (Index i = 0; i < n_substeps; i++)
    {
        PerformSStep(tmp_s_0, lr_sol.S, e_coeff, f_coeff, sigma1, sigma2, blas, mysystem, grid, partition1, partition2, w_x_dep, 0.5 * tau * tau_sub);
        lr_sol.S -= tmp_s_0;
    }
    get_time::stop("sstep");

    /////////////////////////////////////////////
    ///////////////// 1/2 K-STEP ////////////////
    /////////////////////////////////////////////

    get_time::start("kstep");
    get_time::start("kstep_coeff");
    CalculateCoefficientsKL<1>(c_coeff1, d_coeff1, sigma2, lr_sol, blas, mysystem, grid, partition1, partition2, w_x_dep);
    get_time::stop("kstep_coeff");
    tmp_x1_0 = lr_sol.X;
    blas.matmul(tmp_x1_0, lr_sol.S, lr_sol.X); // lr_sol.X contains now K
    for (Index i = 0; i < n_substeps; i++)
    {
        PerformKLStep<1>(tmp_x1_0, lr_sol.X, c_coeff1, d_coeff1, sigma1, blas, mysystem, grid, partition1, w_x_dep, 0.5 * tau * tau_sub);
        lr_sol.X += tmp_x1_0;
    }
    get_time::stop("kstep");

    // Perform the QR decomposition K = X1 * S
    gs(lr_sol.X, lr_sol.S, ip_xx1);

    /////////////////////////////////////////////
    /////////////// NORMALIZATION ///////////////
    /////////////////////////////////////////////

    // Renormalize S
    norm = CalculateNorm(lr_sol, grid);
    lr_sol.S /= norm;
    cout << "norm(P): " << norm << endl;
}


// ///////////////////////////////////////////////////////
// //////////////////// Modified Euler ///////////////////
// ///////////////////////////////////////////////////////

// void IntegrateSecondOrder(lr2<double> &lr_sol, const vector<multi_array<double, 2>> &w_x_dep, vector<multi_array<double, 3>> &c_coeff1, vector<multi_array<double, 3>> &d_coeff1, vector<multi_array<double, 3>> &c_coeff2, vector<multi_array<double, 3>> &d_coeff2, multi_array<double, 5> &e_coeff, multi_array<double, 5> &f_coeff, const vector<Index> sigma1, const vector<Index> sigma2, mysys &mysystem, grid_info &grid, partition_info<1> &partition1, partition_info<2> &partition2, std::function<double(double *, double *)> ip_xx1, std::function<double(double *, double *)> ip_xx2, blas_ops &blas, double tau)
// {
//     gram_schmidt gs(&blas);
//     double norm;

//     // Temporary objects for multiplication and integration
//     multi_array<double, 2> tmp_x1_0({grid.dx1, grid.r});
//     multi_array<double, 2> tmp_x1_1({grid.dx1, grid.r});
//     multi_array<double, 2> tmp_x2_0({grid.dx2, grid.r});
//     multi_array<double, 2> tmp_x2_1({grid.dx2, grid.r});
//     multi_array<double, 2> tmp_s_0({grid.r, grid.r});
//     multi_array<double, 2> tmp_s_1({grid.r, grid.r});

//     /////////////////////////////////////////////
//     ///////////////// 1/2 K-STEP ////////////////
//     /////////////////////////////////////////////

//     get_time::start("kstep");
//     get_time::start("kstep_coeff");
//     CalculateCoefficientsKL<1>(c_coeff1, d_coeff1, sigma2, lr_sol, blas, mysystem, grid, partition1, partition2, w_x_dep);
//     get_time::stop("kstep_coeff");
//     tmp_x1_1 = lr_sol.X;
//     blas.matmul(tmp_x1_1, lr_sol.S, lr_sol.X); // lr_sol.X contains now K
//     tmp_x1_1 = lr_sol.X;
//     PerformKLStep<1>(tmp_x1_0, lr_sol.X, c_coeff1, d_coeff1, sigma1, blas, mysystem, grid, partition1, w_x_dep, 0.5 * 0.5 * tau);
//     tmp_x1_1 += tmp_x1_0;
//     PerformKLStep<1>(tmp_x1_0, tmp_x1_1, c_coeff1, d_coeff1, sigma1, blas, mysystem, grid, partition1, w_x_dep, 0.5 * tau);
//     lr_sol.X += tmp_x1_0;
//     get_time::stop("kstep");

//     // Perform the QR decomposition K = X1 * S
//     gs(lr_sol.X, lr_sol.S, ip_xx1);

//     /////////////////////////////////////////////
//     //////////////// 1/2 S-STEP /////////////////
//     /////////////////////////////////////////////

//     get_time::start("sstep");
//     get_time::start("sstep_coeff");
//     CalculateCoefficientsS(e_coeff, f_coeff, sigma1, sigma2, lr_sol, blas, mysystem, grid, partition1, partition2, w_x_dep);
//     get_time::stop("sstep_coeff");
//     tmp_s_1 = lr_sol.S;
//     PerformSStep(tmp_s_0, lr_sol.S, e_coeff, f_coeff, sigma1, sigma2, blas, mysystem, grid, partition1, partition2, w_x_dep, 0.5 * 0.5 * tau);
//     tmp_s_1 -= tmp_s_0;
//     PerformSStep(tmp_s_0, tmp_s_1, e_coeff, f_coeff, sigma1, sigma2, blas, mysystem, grid, partition1, partition2, w_x_dep, 0.5 * tau);
//     lr_sol.S -= tmp_s_0;
//     get_time::stop("sstep");

//     /////////////////////////////////////////////
//     ////////////////// L-STEP ///////////////////
//     /////////////////////////////////////////////

//     get_time::start("lstep");
//     get_time::start("lstep_coeff");
//     CalculateCoefficientsKL<2>(c_coeff2, d_coeff2, sigma1, lr_sol, blas, mysystem, grid, partition2, partition1, w_x_dep);
//     get_time::stop("lstep_coeff");
//     tmp_x2_1 = lr_sol.V;
//     blas.matmul_transb(tmp_x2_1, lr_sol.S, lr_sol.V); // lr_sol.V contains now L
//     tmp_x2_1 = lr_sol.V;
//     PerformKLStep<2>(tmp_x2_0, lr_sol.V, c_coeff2, d_coeff2, sigma2, blas, mysystem, grid, partition2, w_x_dep, 0.5 * tau);
//     tmp_x2_1 += tmp_x2_0;
//     PerformKLStep<2>(tmp_x2_0, tmp_x2_1, c_coeff2, d_coeff2, sigma2, blas, mysystem, grid, partition2, w_x_dep, tau);
//     lr_sol.V += tmp_x2_0;
//     get_time::stop("lstep");

//     // Perform the QR decomposition L = X2 * S^T
//     gs(lr_sol.V, lr_sol.S, ip_xx2);
//     transpose_inplace(lr_sol.S);

//     /////////////////////////////////////////////
//     //////////////// 1/2 S-STEP /////////////////
//     /////////////////////////////////////////////

//     get_time::start("sstep");
//     get_time::start("sstep_coeff");
//     CalculateCoefficientsS(e_coeff, f_coeff, sigma1, sigma2, lr_sol, blas, mysystem, grid, partition1, partition2, w_x_dep);
//     get_time::stop("sstep_coeff");
//     tmp_s_1 = lr_sol.S;
//     PerformSStep(tmp_s_0, lr_sol.S, e_coeff, f_coeff, sigma1, sigma2, blas, mysystem, grid, partition1, partition2, w_x_dep, 0.5 * 0.5 * tau);
//     tmp_s_1 -= tmp_s_0;
//     PerformSStep(tmp_s_0, tmp_s_1, e_coeff, f_coeff, sigma1, sigma2, blas, mysystem, grid, partition1, partition2, w_x_dep, 0.5 * tau);
//     lr_sol.S -= tmp_s_0;
//     get_time::stop("sstep");

//     /////////////////////////////////////////////
//     ///////////////// 1/2 K-STEP ////////////////
//     /////////////////////////////////////////////

//     get_time::start("kstep");
//     get_time::start("kstep_coeff");
//     CalculateCoefficientsKL<1>(c_coeff1, d_coeff1, sigma2, lr_sol, blas, mysystem, grid, partition1, partition2, w_x_dep);
//     get_time::stop("kstep_coeff");
//     tmp_x1_1 = lr_sol.X;
//     blas.matmul(tmp_x1_1, lr_sol.S, lr_sol.X); // lr_sol.X contains now K
//     tmp_x1_1 = lr_sol.X;
//     PerformKLStep<1>(tmp_x1_0, lr_sol.X, c_coeff1, d_coeff1, sigma1, blas, mysystem, grid, partition1, w_x_dep, 0.5 * 0.5 * tau);
//     tmp_x1_1 += tmp_x1_0;
//     PerformKLStep<1>(tmp_x1_0, tmp_x1_1, c_coeff1, d_coeff1, sigma1, blas, mysystem, grid, partition1, w_x_dep, 0.5 * tau);
//     lr_sol.X += tmp_x1_0;
//     get_time::stop("kstep");

//     // Perform the QR decomposition K = X1 * S
//     gs(lr_sol.X, lr_sol.S, ip_xx1);

//     /////////////////////////////////////////////
//     /////////////// NORMALIZATION ///////////////
//     /////////////////////////////////////////////

//     // Renormalize S
//     norm = CalculateNorm(lr_sol, grid);
//     lr_sol.S /= norm;
//     cout << "norm(P): " << norm << endl;
// }


// ///////////////////////////////////////////////////////
// //////////////////////// Heun /////////////////////////
// ///////////////////////////////////////////////////////

// void IntegrateSecondOrder(lr2<double> &lr_sol, const vector<multi_array<double, 2>> &w_x_dep, vector<multi_array<double, 3>> &c_coeff1, vector<multi_array<double, 3>> &d_coeff1, vector<multi_array<double, 3>> &c_coeff2, vector<multi_array<double, 3>> &d_coeff2, multi_array<double, 5> &e_coeff, multi_array<double, 5> &f_coeff, const vector<Index> sigma1, const vector<Index> sigma2, mysys &mysystem, grid_info &grid, partition_info<1> &partition1, partition_info<2> &partition2, std::function<double(double *, double *)> ip_xx1, std::function<double(double *, double *)> ip_xx2, blas_ops &blas, double tau)
// {
//     gram_schmidt gs(&blas);
//     double norm;

//     // Temporary objects for multiplication and integration
//     multi_array<double, 2> tmp_x1({grid.dx1, grid.r});
//     multi_array<double, 2> x1_k1({grid.dx1, grid.r});
//     multi_array<double, 2> x1_k2({grid.dx1, grid.r});

//     multi_array<double, 2> tmp_x2({grid.dx2, grid.r});
//     multi_array<double, 2> x2_k1({grid.dx2, grid.r});
//     multi_array<double, 2> x2_k2({grid.dx2, grid.r});

//     multi_array<double, 2> tmp_s({grid.r, grid.r});
//     multi_array<double, 2> s_k1({grid.r, grid.r});
//     multi_array<double, 2> s_k2({grid.r, grid.r});

//     /////////////////////////////////////////////
//     ///////////////// 1/2 K-STEP ////////////////
//     /////////////////////////////////////////////

//     get_time::start("kstep");
//     get_time::start("kstep_coeff");
//     CalculateCoefficientsKL<1>(c_coeff1, d_coeff1, sigma2, lr_sol, blas, mysystem, grid, partition1, partition2, w_x_dep);
//     get_time::stop("kstep_coeff");
//     tmp_x1 = lr_sol.X;
//     blas.matmul(tmp_x1, lr_sol.S, lr_sol.X); // lr_sol.X contains now K

//     PerformKLStep<1>(x1_k1, lr_sol.X, c_coeff1, d_coeff1, sigma1, blas, mysystem, grid, partition1, w_x_dep, 0.5 * tau);
//     tmp_x1 = lr_sol.X + x1_k1;
//     PerformKLStep<1>(x1_k2, tmp_x1, c_coeff1, d_coeff1, sigma1, blas, mysystem, grid, partition1, w_x_dep, 0.5 * tau);
//     lr_sol.X += x1_k1 * 0.5 + x1_k2 * 0.5;

//     get_time::stop("kstep");

//     // Perform the QR decomposition K = X1 * S
//     gs(lr_sol.X, lr_sol.S, ip_xx1);

//     /////////////////////////////////////////////
//     //////////////// 1/2 S-STEP /////////////////
//     /////////////////////////////////////////////

//     get_time::start("sstep");
//     get_time::start("sstep_coeff");
//     CalculateCoefficientsS(e_coeff, f_coeff, sigma1, sigma2, lr_sol, blas, mysystem, grid, partition1, partition2, w_x_dep);
//     get_time::stop("sstep_coeff");

//     PerformSStep(s_k1, lr_sol.S, e_coeff, f_coeff, sigma1, sigma2, blas, mysystem, grid, partition1, partition2, w_x_dep, -0.5 * tau);
//     tmp_s = lr_sol.S + s_k1;
//     PerformSStep(s_k2, tmp_s, e_coeff, f_coeff, sigma1, sigma2, blas, mysystem, grid, partition1, partition2, w_x_dep, -0.5 * tau);
//     lr_sol.S += s_k1 * 0.5 + s_k2 * 0.5;

//     get_time::stop("sstep");

//     /////////////////////////////////////////////
//     ////////////////// L-STEP ///////////////////
//     /////////////////////////////////////////////

//     get_time::start("lstep");
//     get_time::start("lstep_coeff");
//     CalculateCoefficientsKL<2>(c_coeff2, d_coeff2, sigma1, lr_sol, blas, mysystem, grid, partition2, partition1, w_x_dep);
//     get_time::stop("lstep_coeff");
//     tmp_x2 = lr_sol.V;
//     blas.matmul_transb(tmp_x2, lr_sol.S, lr_sol.V); // lr_sol.V contains now L

//     PerformKLStep<2>(x2_k1, lr_sol.V, c_coeff2, d_coeff2, sigma2, blas, mysystem, grid, partition2, w_x_dep, tau);
//     tmp_x2 = lr_sol.V + x2_k1;
//     PerformKLStep<2>(x2_k2, tmp_x2, c_coeff2, d_coeff2, sigma2, blas, mysystem, grid, partition2, w_x_dep, tau);
//     lr_sol.V += x2_k1 * 0.5 + x2_k2 * 0.5;

//     get_time::stop("lstep");

//     // Perform the QR decomposition L = X2 * S^T
//     gs(lr_sol.V, lr_sol.S, ip_xx2);
//     transpose_inplace(lr_sol.S);

//     /////////////////////////////////////////////
//     //////////////// 1/2 S-STEP /////////////////
//     /////////////////////////////////////////////

//     get_time::start("sstep");
//     get_time::start("sstep_coeff");
//     CalculateCoefficientsS(e_coeff, f_coeff, sigma1, sigma2, lr_sol, blas, mysystem, grid, partition1, partition2, w_x_dep);
//     get_time::stop("sstep_coeff");

//     PerformSStep(s_k1, lr_sol.S, e_coeff, f_coeff, sigma1, sigma2, blas, mysystem, grid, partition1, partition2, w_x_dep, -0.5 * tau);
//     tmp_s = lr_sol.S + s_k1;
//     PerformSStep(s_k2, tmp_s, e_coeff, f_coeff, sigma1, sigma2, blas, mysystem, grid, partition1, partition2, w_x_dep, -0.5 * tau);
//     lr_sol.S += s_k1 * 0.5 + s_k2 * 0.5;

//     get_time::stop("sstep");

//     /////////////////////////////////////////////
//     ///////////////// 1/2 K-STEP ////////////////
//     /////////////////////////////////////////////

//     get_time::start("kstep");
//     get_time::start("kstep_coeff");
//     CalculateCoefficientsKL<1>(c_coeff1, d_coeff1, sigma2, lr_sol, blas, mysystem, grid, partition1, partition2, w_x_dep);
//     get_time::stop("kstep_coeff");
//     tmp_x1 = lr_sol.X;
//     blas.matmul(tmp_x1, lr_sol.S, lr_sol.X); // lr_sol.X contains now K

//     PerformKLStep<1>(x1_k1, lr_sol.X, c_coeff1, d_coeff1, sigma1, blas, mysystem, grid, partition1, w_x_dep, 0.5 * tau);
//     tmp_x1 = lr_sol.X + x1_k1;
//     PerformKLStep<1>(x1_k2, tmp_x1, c_coeff1, d_coeff1, sigma1, blas, mysystem, grid, partition1, w_x_dep, 0.5 * tau);
//     lr_sol.X += x1_k1 * 0.5 + x1_k2 * 0.5;

//     get_time::stop("kstep");

//     // Perform the QR decomposition K = X1 * S
//     gs(lr_sol.X, lr_sol.S, ip_xx1);

//     /////////////////////////////////////////////
//     /////////////// NORMALIZATION ///////////////
//     /////////////////////////////////////////////

//     // Renormalize S
//     norm = CalculateNorm(lr_sol, grid);
//     lr_sol.S /= norm;
//     cout << "norm(P): " << norm << endl;
// }


// ///////////////////////////////////////////////////////
// ///////////////////////// RK4 /////////////////////////
// ///////////////////////////////////////////////////////

// void IntegrateSecondOrder(lr2<double> &lr_sol, const vector<multi_array<double, 2>> &w_x_dep, vector<multi_array<double, 3>> &c_coeff1, vector<multi_array<double, 3>> &d_coeff1, vector<multi_array<double, 3>> &c_coeff2, vector<multi_array<double, 3>> &d_coeff2, multi_array<double, 5> &e_coeff, multi_array<double, 5> &f_coeff, const vector<Index> sigma1, const vector<Index> sigma2, mysys &mysystem, grid_info &grid, partition_info<1> &partition1, partition_info<2> &partition2, std::function<double(double *, double *)> ip_xx1, std::function<double(double *, double *)> ip_xx2, blas_ops &blas, double tau)
// {
//     gram_schmidt gs(&blas);
//     double norm;
//     double six_inv = 1.0 / 6.0;

//     // Temporary objects for multiplication and integration
//     multi_array<double, 2> tmp_x1({grid.dx1, grid.r});
//     multi_array<double, 2> x1_k1({grid.dx1, grid.r});
//     multi_array<double, 2> x1_k2({grid.dx1, grid.r});
//     multi_array<double, 2> x1_k3({grid.dx1, grid.r});
//     multi_array<double, 2> x1_k4({grid.dx1, grid.r});

//     multi_array<double, 2> tmp_x2({grid.dx2, grid.r});
//     multi_array<double, 2> x2_k1({grid.dx2, grid.r});
//     multi_array<double, 2> x2_k2({grid.dx2, grid.r});
//     multi_array<double, 2> x2_k3({grid.dx2, grid.r});
//     multi_array<double, 2> x2_k4({grid.dx2, grid.r});

//     multi_array<double, 2> tmp_s({grid.r, grid.r});
//     multi_array<double, 2> s_k1({grid.r, grid.r});
//     multi_array<double, 2> s_k2({grid.r, grid.r});
//     multi_array<double, 2> s_k3({grid.r, grid.r});
//     multi_array<double, 2> s_k4({grid.r, grid.r});

//     /////////////////////////////////////////////
//     ///////////////// 1/2 K-STEP ////////////////
//     /////////////////////////////////////////////

//     get_time::start("kstep");
//     get_time::start("kstep_coeff");
//     CalculateCoefficientsKL<1>(c_coeff1, d_coeff1, sigma2, lr_sol, blas, mysystem, grid, partition1, partition2, w_x_dep);
//     get_time::stop("kstep_coeff");
//     tmp_x1 = lr_sol.X;
//     blas.matmul(tmp_x1, lr_sol.S, lr_sol.X); // lr_sol.X contains now K

//     PerformKLStep<1>(x1_k1, lr_sol.X, c_coeff1, d_coeff1, sigma1, blas, mysystem, grid, partition1, w_x_dep, 0.5 * tau);
//     tmp_x1 = lr_sol.X + x1_k1 * 0.5;
//     PerformKLStep<1>(x1_k2, tmp_x1, c_coeff1, d_coeff1, sigma1, blas, mysystem, grid, partition1, w_x_dep, 0.5 * tau);
//     tmp_x1 = lr_sol.X + x1_k2 * 0.5;
//     PerformKLStep<1>(x1_k3, tmp_x1, c_coeff1, d_coeff1, sigma1, blas, mysystem, grid, partition1, w_x_dep, 0.5 * tau);
//     tmp_x1 = lr_sol.X + x1_k3;
//     PerformKLStep<1>(x1_k4, tmp_x1, c_coeff1, d_coeff1, sigma1, blas, mysystem, grid, partition1, w_x_dep, 0.5 * tau);
//     tmp_x1 = x1_k1 + x1_k2 * 2.0 + x1_k3 * 2.0 + x1_k4;
//     lr_sol.X += tmp_x1 * six_inv;

//     get_time::stop("kstep");

//     // Perform the QR decomposition K = X1 * S
//     gs(lr_sol.X, lr_sol.S, ip_xx1);

//     /////////////////////////////////////////////
//     //////////////// 1/2 S-STEP /////////////////
//     /////////////////////////////////////////////

//     get_time::start("sstep");
//     get_time::start("sstep_coeff");
//     CalculateCoefficientsS(e_coeff, f_coeff, sigma1, sigma2, lr_sol, blas, mysystem, grid, partition1, partition2, w_x_dep);
//     get_time::stop("sstep_coeff");

//     PerformSStep(s_k1, lr_sol.S, e_coeff, f_coeff, sigma1, sigma2, blas, mysystem, grid, partition1, partition2, w_x_dep, -0.5 * tau);
//     tmp_s = lr_sol.S + s_k1 * 0.5;
//     PerformSStep(s_k2, tmp_s, e_coeff, f_coeff, sigma1, sigma2, blas, mysystem, grid, partition1, partition2, w_x_dep, -0.5 * tau);
//     tmp_s = lr_sol.S + s_k2 * 0.5;
//     PerformSStep(s_k3, tmp_s, e_coeff, f_coeff, sigma1, sigma2, blas, mysystem, grid, partition1, partition2, w_x_dep, -0.5 * tau);
//     tmp_s = lr_sol.S + s_k3;
//     PerformSStep(s_k4, tmp_s, e_coeff, f_coeff, sigma1, sigma2, blas, mysystem, grid, partition1, partition2, w_x_dep, -0.5 * tau);
//     tmp_s = s_k1 + s_k2 * 2.0 + s_k3 * 2.0 + s_k4;
//     lr_sol.S += tmp_s * six_inv;

//     get_time::stop("sstep");

//     /////////////////////////////////////////////
//     ////////////////// L-STEP ///////////////////
//     /////////////////////////////////////////////

//     get_time::start("lstep");
//     get_time::start("lstep_coeff");
//     CalculateCoefficientsKL<2>(c_coeff2, d_coeff2, sigma1, lr_sol, blas, mysystem, grid, partition2, partition1, w_x_dep);
//     get_time::stop("lstep_coeff");
//     tmp_x2 = lr_sol.V;
//     blas.matmul_transb(tmp_x2, lr_sol.S, lr_sol.V); // lr_sol.V contains now L

//     PerformKLStep<2>(x2_k1, lr_sol.V, c_coeff2, d_coeff2, sigma2, blas, mysystem, grid, partition2, w_x_dep, tau);
//     tmp_x2 = lr_sol.V + x2_k1 * 0.5;
//     PerformKLStep<2>(x2_k2, tmp_x2, c_coeff2, d_coeff2, sigma2, blas, mysystem, grid, partition2, w_x_dep, tau);
//     tmp_x2 = lr_sol.V + x2_k2 * 0.5;
//     PerformKLStep<2>(x2_k3, tmp_x2, c_coeff2, d_coeff2, sigma2, blas, mysystem, grid, partition2, w_x_dep, tau);
//     tmp_x2 = lr_sol.V + x2_k3;
//     PerformKLStep<2>(x2_k4, tmp_x2, c_coeff2, d_coeff2, sigma2, blas, mysystem, grid, partition2, w_x_dep, tau);
//     tmp_x2 = x2_k1 + x2_k2 * 2.0 + x2_k3 * 2.0 + x2_k4;
//     lr_sol.V += tmp_x2 * six_inv;

//     get_time::stop("lstep");

//     // Perform the QR decomposition L = X2 * S^T
//     gs(lr_sol.V, lr_sol.S, ip_xx2);
//     transpose_inplace(lr_sol.S);

//     /////////////////////////////////////////////
//     //////////////// 1/2 S-STEP /////////////////
//     /////////////////////////////////////////////

//     get_time::start("sstep");
//     get_time::start("sstep_coeff");
//     CalculateCoefficientsS(e_coeff, f_coeff, sigma1, sigma2, lr_sol, blas, mysystem, grid, partition1, partition2, w_x_dep);
//     get_time::stop("sstep_coeff");

//     PerformSStep(s_k1, lr_sol.S, e_coeff, f_coeff, sigma1, sigma2, blas, mysystem, grid, partition1, partition2, w_x_dep, -0.5 * tau);
//     tmp_s = lr_sol.S + s_k1 * 0.5;
//     PerformSStep(s_k2, tmp_s, e_coeff, f_coeff, sigma1, sigma2, blas, mysystem, grid, partition1, partition2, w_x_dep, -0.5 * tau);
//     tmp_s = lr_sol.S + s_k2 * 0.5;
//     PerformSStep(s_k3, tmp_s, e_coeff, f_coeff, sigma1, sigma2, blas, mysystem, grid, partition1, partition2, w_x_dep, -0.5 * tau);
//     tmp_s = lr_sol.S + s_k3;
//     PerformSStep(s_k4, tmp_s, e_coeff, f_coeff, sigma1, sigma2, blas, mysystem, grid, partition1, partition2, w_x_dep, -0.5 * tau);
//     tmp_s = s_k1 + s_k2 * 2.0 + s_k3 * 2.0 + s_k4;
//     lr_sol.S += tmp_s * six_inv;

//     get_time::stop("sstep");

//     /////////////////////////////////////////////
//     ///////////////// 1/2 K-STEP ////////////////
//     /////////////////////////////////////////////

//     get_time::start("kstep");
//     get_time::start("kstep_coeff");
//     CalculateCoefficientsKL<1>(c_coeff1, d_coeff1, sigma2, lr_sol, blas, mysystem, grid, partition1, partition2, w_x_dep);
//     get_time::stop("kstep_coeff");
//     tmp_x1 = lr_sol.X;
//     blas.matmul(tmp_x1, lr_sol.S, lr_sol.X); // lr_sol.X contains now K

//     PerformKLStep<1>(x1_k1, lr_sol.X, c_coeff1, d_coeff1, sigma1, blas, mysystem, grid, partition1, w_x_dep, 0.5 * tau);
//     tmp_x1 = lr_sol.X + x1_k1 * 0.5;
//     PerformKLStep<1>(x1_k2, tmp_x1, c_coeff1, d_coeff1, sigma1, blas, mysystem, grid, partition1, w_x_dep, 0.5 * tau);
//     tmp_x1 = lr_sol.X + x1_k2 * 0.5;
//     PerformKLStep<1>(x1_k3, tmp_x1, c_coeff1, d_coeff1, sigma1, blas, mysystem, grid, partition1, w_x_dep, 0.5 * tau);
//     tmp_x1 = lr_sol.X + x1_k3;
//     PerformKLStep<1>(x1_k4, tmp_x1, c_coeff1, d_coeff1, sigma1, blas, mysystem, grid, partition1, w_x_dep, 0.5 * tau);
//     tmp_x1 = x1_k1 + x1_k2 * 2.0 + x1_k3 * 2.0 + x1_k4;
//     lr_sol.X += tmp_x1 * six_inv;

//     get_time::stop("kstep");

//     // Perform the QR decomposition K = X1 * S
//     gs(lr_sol.X, lr_sol.S, ip_xx1);

//     /////////////////////////////////////////////
//     /////////////// NORMALIZATION ///////////////
//     /////////////////////////////////////////////

//     // Renormalize S
//     norm = CalculateNorm(lr_sol, grid);
//     lr_sol.S /= norm;
//     cout << "norm(P): " << norm << endl;
// }