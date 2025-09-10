#include "print_functions.hpp"

// TODO: memory requirement
std::ostream& operator<<(std::ostream& os, const diagnostics& dgn)
{
    const auto [hrs, mins, secs, ms] = ChronoBurst(dgn.t_elapsed);

    os << "DIAGNOSTICS\n"
        << "-----------\n"
        << "Time elapsed: "
        << hrs.count() << "h "
        << mins.count() << "mins "
        << secs.count() << "s "
        << ms.count() << "ms\n"
        << "Integration method (K): " << dgn.integrator.integration_methods.at("K")->get_name() << "\n"
        << "Integration method (S): " << dgn.integrator.integration_methods.at("S")->get_name() << "\n"
        << "Integration method (Q): " << dgn.integrator.integration_methods.at("Q")->get_name() << "\n"
        << "Time step size: " << dgn.tau << "\n"
        << "max(norm - 1.0): " << dgn.dm_max << "\n";
#ifdef __OPENMP__
    os << "[OpenMP activated]: OMP_NUM_THREADS=" << omp_get_max_threads() << "\n";
#else
    os << "[OpenMP not activated]\n";
#endif
    os << "-----------\n" << endl;

    return os;
}