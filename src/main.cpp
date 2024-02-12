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
        ("o,output", "Name of the output folder", cxxopts::value<std::string>())
        ("s,snapshot", "Number of steps between two snapshots", cxxopts::value<int>())
        ("t,tau", "Time step size", cxxopts::value<double>())
        ("f,tfinal", "Final integration time", cxxopts::value<double>())
        ("h,help", "Print usage")
        ;

    auto result = options.parse(argc, argv);
    if (result.count("help"))
    {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    std::string output = result["output"].as<std::string>();
    int snapshot = result["snapshot"].as<int>();
    double tau = result["tau"].as<double>();
    double tfinal = result["tfinal"].as<double>();

    get_time::start("main");

    std::map<std::string, integration_method*> integrations_methods;
    integrations_methods["K"] = new implicit_euler{};
    integrations_methods["S"] = new explicit_euler{};
    integrations_methods["Q"] = new implicit_euler{};

    blas_ops blas;
    TTNIntegrator integrator(blas, integrations_methods);
    cme_lr_tree tree;

    double t = 0.0;
    double dm = 0.0;
    double dm_max = 0.0;

    const Index kNsteps = ceil(tfinal / tau); // Number of time steps

    tree.Read("input/input.nc");
    std::cout << tree;
    tree.Orthogonalize(blas);
    double norm = tree.Normalize();
    std::cout << "Norm: " << norm << std::endl;

    // Check if folder in ../output/ exists, otherwise create folder
    std::string fname;
    fname = "output/" + output;
    std::filesystem::create_directory(fname);

    // Store initial values
    fname = "output/" + output + "/output_t0.nc";
    tree.Write(fname, t, tau, dm);

    auto t_start(std::chrono::high_resolution_clock::now());

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

    auto t_stop(std::chrono::high_resolution_clock::now());
    auto t_elapsed = t_stop - t_start;

    get_time::stop("main");
    std::cout << "\n\n";
    std::cout << "TIMER RESULTS\n";
    std::cout << "-------------\n";
    std::cout << get_time::sorted_output();

    PrintDiagnostics(integrator.integration_methods, t_elapsed, tau, dm_max);

    return EXIT_SUCCESS;
}