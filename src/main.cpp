#include <filesystem>

#include <cxxopts.hpp>

#include "print_functions.hpp"
#include "subflows.hpp"
#include "timer_class.hpp"
#include "tree_class.hpp"

int main(int argc, char* argv[])
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

    double t = 0.0;
    blas_ops blas;
    implicit_euler method;
    cme_lr_tree tree;
    const Index kNsteps = ceil(tfinal / tau); // Number of time steps

    tree.Read("input/input.nc");
    cout << tree;
    tree.Orthogonalize(blas);
    double norm = tree.Normalize();
    cout << "Norm: " << norm << endl;

    // Check if folder in ../output/ exists, otherwise create folder
    std::string fname;
    fname = "output/" + output;
    std::filesystem::create_directory(fname);

    // Store initial values
    fname = "output/" + output + "/output_t0.nc";
    tree.Write(fname, t, tau);

    auto start_time(std::chrono::high_resolution_clock::now());

    for (Index ts = 0; ts < kNsteps; ++ts)
    {
        if (tfinal - t < tau)
            tau = tfinal - t;

        TTNIntegrator(tree.root, blas, tau, method);
        norm = tree.Normalize();

        t += tau;

        // Print progress bar
        PrintProgressBar(ts, kNsteps, start_time, norm);

        // Write snapshot
        if ((ts + 1) % snapshot == 0 || (ts + 1) == kNsteps)
        {
            fname = "output/" + output + "/output_t" + std::to_string(ts + 1) + ".nc";
            tree.Write(fname, t, tau);
        }
    }

    get_time::stop("main");
    cout << endl << endl;
    cout << "TIMER RESULTS" << endl;
    cout << "-------------" << endl;
    cout << get_time::sorted_output();

    return 0;
}