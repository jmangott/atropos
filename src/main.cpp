#include <cstdlib>
#include <filesystem>

#include <cxxopts.hpp>

#include "integrators.hpp"
#include "print_functions.hpp"
#include "timer_class.hpp"
#include "tree_class.hpp"

int main(int argc, char** argv)
{
    cxxopts::Options options("hierarchical-cme", "Tree tensor network integrator for the chemical master equation");

    options.add_options()
        ("i,input", "Name of the input .nc file", cxxopts::value<std::string>()->default_value("input/input.nc"))
        ("o,output", "Name of the output folder", cxxopts::value<std::string>())
        ("s,snapshot", "Number of steps between two snapshots", cxxopts::value<int>())
        ("t,tau", "Time step size", cxxopts::value<double>())
        ("f,tfinal", "Final integration time", cxxopts::value<double>())
        ("n,substeps", "Number of integration substeps", cxxopts::value<unsigned int>()->default_value("1"))
        ("m,method", "Integration method (`e` (explicit Euler), `i` (implicit Euler), `c` (Crank-Nicolson))", cxxopts::value<char>()->default_value("i"))
        ("h,help", "Print usage")
        ;

    auto result = options.parse(argc, argv);
    if (result.count("help"))
    {
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
    switch (method)
    {
    case 'e':
        integration_methods["K"] = new explicit_euler(substeps);
        integration_methods["Q"] = new explicit_euler(substeps);
        integration_methods["S"] = new explicit_euler(substeps);
        break;
    case 'i':
        integration_methods["K"] = new implicit_euler(substeps);
        integration_methods["Q"] = new implicit_euler(substeps);
        integration_methods["S"] = new explicit_euler(substeps);  // S is integrated backwards in time
        break;
    case 'c':
        integration_methods["K"] = new crank_nicolson(substeps);
        integration_methods["Q"] = new crank_nicolson(substeps);
        integration_methods["S"] = new crank_nicolson(substeps);
        break;
    default:
        std::cout << "Error: Command line option `m` must be either `e`, `i` or `c`!" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    blas_ops blas;
    ttn_integrator integrator(blas, integration_methods);
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

    // Check if folder in ../output/ exists, otherwise create folder
    std::string fname;
    fname = "output/" + output;
    std::filesystem::create_directory(fname);

    // Store initial values
    fname = "output/" + output + "/output_t0.nc";
    tree.Write(fname, t, tau, dm);

    auto t_start(std::chrono::high_resolution_clock::now());
    get_time::start("main");
    for (Index ts = 0; ts < kNsteps; ++ts)
    {
        if (tfinal - t < tau)
            tau = tfinal - t;

        integrator(tree.root, tau);
        norm = tree.Normalize();

        dm = norm - 1.0;
        if (std::abs(dm) > std::abs(dm_max)) dm_max = dm;
        t += tau;

        PrintProgressBar(ts, kNsteps, t_start, norm);

        // Write snapshot
        if ((ts + 1) % snapshot == 0 || (ts + 1) == kNsteps)
        {
            fname = "output/" + output + "/output_t" + std::to_string(ts + 1) + ".nc";
            tree.Write(fname, t, tau, dm);
        }
    }
    get_time::stop("main");

    auto t_stop(std::chrono::high_resolution_clock::now());
    auto t_elapsed = t_stop - t_start;

    std::cout << "\n\n";
    std::cout << "TIMER RESULTS\n";
    std::cout << "-------------\n";
    std::cout << get_time::sorted_output();

    std::ofstream diagnostics_file;
    diagnostics dgn{integrator, t_elapsed, tau, dm_max};
    diagnostics_file.open("output/" + output + "/diagnostics.txt", std::fstream::out);
    diagnostics_file << dgn;
    diagnostics_file.close();
    std::cout << dgn;

    return EXIT_SUCCESS;
}