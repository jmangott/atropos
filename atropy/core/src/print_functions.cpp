#include "print_functions.hpp"

// TODO: output width is not always the same
// TODO: calculate `time_left` with a moving mean (of the last n steps)
void PrintProgressBar(const Index ts, const Index kNsteps,
                      const std::chrono::system_clock::time_point t_start,
                      const double norm)
{
    int bar_width = 30;
    double progress = (ts + 1.0) / kNsteps;
    int pos = bar_width * progress;
    string time_unit;
    string progress_bar(bar_width, ' ');

    auto t_stop(std::chrono::high_resolution_clock::now());
    auto duration(std::chrono::duration_cast<std::chrono::seconds>(t_stop - t_start));
    // time_per_step = duration.count() / (ts + 1.0);
    auto time_per_step = duration / (ts + 1.0);
    auto time_left = time_per_step * (kNsteps - 1.0 - ts);
    double time_per_step_count;

    if (time_per_step.count() < 0.01) {
        time_unit = "ms";
        time_per_step_count = time_per_step.count() * 1000.0;
    }
    else if (time_per_step.count() < 60.0) {
        time_unit = "s";
        time_per_step_count = time_per_step.count();
    }
    else {
        time_unit = "min";
        time_per_step_count = time_per_step.count() / 60.0;
    }

    const auto [hrs, mins, secs, ms] = ChronoBurst(time_left);

    std::fill(std::begin(progress_bar), std::begin(progress_bar) + pos, '#');
    std::fill(std::begin(progress_bar) + pos, std::end(progress_bar), '-');

    printf("[%*s], step: %ti/%ti, time per step: %.2f%*s, time left: "
           "%2.2lli:%2.2lli:%2.2lli, progress: %4.2f%%, |norm(P)-1|: %3.2e\r",
           bar_width, progress_bar.c_str(), ts + 1, kNsteps, time_per_step_count,
           (int)time_unit.size(), time_unit.c_str(), (long long int)hrs.count(),
           (long long int)mins.count(), (long long int)secs.count(), progress * 100,
           std::abs(norm - 1.0));
    fflush(stdout);
}

// TODO: memory requirement
std::ostream& operator<<(std::ostream& os, const diagnostics& dgn)
{
    const auto [hrs, mins, secs, ms] = ChronoBurst(dgn.t_elapsed);

    os << "DIAGNOSTICS\n"
       << "-----------\n"
       << "Time elapsed: " << hrs.count() << "h " << mins.count() << "mins "
       << secs.count() << "s " << ms.count() << "ms\n"
       << "Integration method (K): "
       << dgn.integrator.integration_methods.at("K")->get_name() << "\n"
       << "Integration method (S): "
       << dgn.integrator.integration_methods.at("S")->get_name() << "\n"
       << "Integration method (Q): "
       << dgn.integrator.integration_methods.at("Q")->get_name() << "\n"
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
