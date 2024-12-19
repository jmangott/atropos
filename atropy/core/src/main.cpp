#include <cstdlib>
#include <filesystem>

#include <cxxopts.hpp>
#include <generic/timer.hpp>

#include "bug_integrator.hpp"
#include "print_functions.hpp"
#include "ps_integrator.hpp"
#include "tree_class.hpp"

#define STRING(x) #x
#define XSTRING(x) STRING(x)

int main(int argc, char** argv)
{
    std::string source_root;

#ifdef SOURCE_ROOT
    source_root = XSTRING(SOURCE_ROOT);
#else
    std::cout << "Warning: SOURCE_ROOT not set, default directories of input and "
                 "output are set to current working directory"
              << std::endl;
    source_root = ".";
#endif

    cxxopts::Options options(
        "hierarchical-cme",
        "Tree tensor network integrator for the chemical master equation");

    options.add_options()(
        "i,input", "Name of the input .nc file",
        cxxopts::value<std::string>()->default_value(source_root + "/input/input.nc"))(
        "o,output", "Name of the output folder", cxxopts::value<std::string>())(
        "s,snapshot", "Number of steps between two snapshots",
        cxxopts::value<int>())("t,tau", "Time step size", cxxopts::value<double>())(
        "f,tfinal", "Final integration time",
        cxxopts::value<double>())("n,substeps", "Number of integration substeps",
                                  cxxopts::value<unsigned int>()->default_value("1"))(
        "m,method",
        "Integration method (`e` (explicit Euler), `r` (explicit RK4), `i` (implicit "
        "Euler), `c` (Crank-Nicolson))",
        cxxopts::value<char>()->default_value("i"))("h,help", "Print usage");

    auto result = options.parse(argc, argv);
    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        std::exit(EXIT_SUCCESS);
    }

    std::string input = result["input"].as<std::string>();
    std::string output = result["output"].as<std::string>();
    int snapshot = result["snapshot"].as<int>();
    double tau = result["tau"].as<double>();
    double tfinal = result["tfinal"].as<double>();
    unsigned int substeps = result["substeps"].as<unsigned int>();
    char method = result["method"].as<char>();

    std::map<std::string, integration_method*> integration_methods;
    switch (method) {
    case 'e':
        integration_methods["K"] = new explicit_euler(substeps);
        integration_methods["Q"] = new explicit_euler(substeps);
        integration_methods["S"] = new explicit_euler(substeps);
        break;
    case 'r':
        integration_methods["K"] = new rk4(substeps);
        integration_methods["Q"] = new rk4(substeps);
        integration_methods["S"] = new rk4(substeps);
        break;
    case 'i':
        integration_methods["K"] = new implicit_euler(substeps);
        integration_methods["Q"] = new explicit_euler(substeps);
        integration_methods["S"] =
            new explicit_euler(substeps); // S is integrated backwards in time
        break;
    case 'c':
        integration_methods["K"] = new crank_nicolson(substeps);
        integration_methods["Q"] = new crank_nicolson(substeps);
        integration_methods["S"] = new crank_nicolson(substeps);
        break;
    default:
        std::cout
            << "Error: Command line option `m` must be either `e`, `r`, `i` or `c`!"
            << std::endl;
        std::exit(EXIT_FAILURE);
    }

    Ensign::Matrix::blas_ops blas;
    ps_integrator integrator(blas, integration_methods);
    cme_lr_tree tree;

    double t = 0.0;
    double dm = 0.0;
    double dm_max = 0.0;

    const Index kNsteps = std::ceil(tfinal / tau);

    tree.Read(input);
    std::cout << tree;
    tree.Orthogonalize(blas);
    double norm = tree.Normalize();
    tree.InitializeAB_bar(blas);
    std::cout << "Norm: " << norm << std::endl;

    // Check if folder in SOURCE_ROOT/output/ exists, otherwise create folder
    std::string output_fname;
    output_fname = source_root + "/output/" + output;
    std::filesystem::create_directory(output_fname);

    // Store initial values
    std::string fname;
    fname = output_fname + "/output_t0.nc";
    tree.Write(fname, t, tau, dm);

    auto t_start(std::chrono::high_resolution_clock::now());
    Ensign::gt::start("main");
    for (Index ts = 0; ts < kNsteps; ++ts) {
        if (tfinal - t < tau)
            tau = tfinal - t;

        integrator(tree.root, tau);
        norm = tree.Normalize();

        dm = norm - 1.0;
        if (std::abs(dm) > std::abs(dm_max))
            dm_max = dm;
        t += tau;

        PrintProgressBar(ts, kNsteps, t_start, norm);

        // Write snapshot
        if ((ts + 1) % snapshot == 0 || (ts + 1) == kNsteps) {
            fname = output_fname + "/output_t" + std::to_string(ts + 1) + ".nc";
            tree.Write(fname, t, tau, dm);
        }
    }
    Ensign::gt::stop("main");

    auto t_stop(std::chrono::high_resolution_clock::now());
    auto t_elapsed = t_stop - t_start;

    std::cout << "\n\n";
    std::cout << "TIMER RESULTS\n";
    std::cout << "-------------\n";
    std::cout << Ensign::gt::sorted_output();

    std::ofstream diagnostics_file;
    diagnostics dgn{integrator, t_elapsed, tau, dm_max};
    diagnostics_file.open(output_fname + "/diagnostics.txt", std::fstream::out);
    diagnostics_file << dgn;
    diagnostics_file.close();
    std::cout << dgn;

    return EXIT_SUCCESS;
}
