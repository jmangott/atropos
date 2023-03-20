#include "print_functions.hpp"

using std::cout;
using std::endl;

// TODO: output width is not always the same
// TODO: calculate `time_left` with a moving mean (of the last n steps)
void PrintProgressBar(Index ts, Index kNsteps, std::chrono::_V2::system_clock::time_point start_time, double norm)
{
    int bar_width = 30;
    double progress = (ts + 1.0) / kNsteps;
    int pos = bar_width * progress;
    string time_unit;
    string progress_bar(bar_width, ' ');

    auto end_time(std::chrono::high_resolution_clock::now());
    auto duration(std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time));
    // time_per_step = duration.count() / (ts + 1.0);
    auto time_per_step = duration / (ts + 1.0);
    auto time_left = time_per_step * (kNsteps - 1 - ts);
    double time_per_step_count;

    if (time_per_step.count() < 0.01)
    {
        time_unit = "ms";
        time_per_step_count = std::chrono::duration_cast<std::chrono::milliseconds>(time_per_step).count();
    }
    else if (time_per_step.count() < 60)
    {
        time_unit = "s";
        time_per_step_count = time_per_step.count();
    }
    else
    {
        time_unit = "min";
        time_per_step_count = std::chrono::duration_cast<std::chrono::minutes>(time_per_step).count();
    }

    // hours = (int) (time_left / 3600);
    // minutes = (int) (std::fmod(time_left, 3600.0) / 60);
    // seconds = (int) std::fmod(time_left, 60.0);

    auto hours = std::chrono::duration_cast<std::chrono::hours>(time_left);
    auto minutes = std::chrono::duration_cast<std::chrono::minutes>(time_left - hours);
    auto seconds = std::chrono::duration_cast<std::chrono::seconds>(time_left - hours - minutes);

    std::fill(progress_bar.begin(), progress_bar.begin() + pos, '#');
    std::fill(progress_bar.begin() + pos, progress_bar.end(), '-');

    printf("[%*s], step: %ti/%ti, time per step: %.2f%*s, time left: %2.2lli:%2.2lli:%2.2lli, progress: %4.2f%%, |norm(P)-1|: %3.2e\r", bar_width, progress_bar.c_str(), ts + 1, kNsteps, time_per_step_count, (int)time_unit.size(), time_unit.c_str(), hours.count(), minutes.count(), seconds.count(), progress * 100, std::abs(norm - 1.0));
    fflush(stdout);
}


// TODO: memory requirement
void PrintDiagnostics(grid_info grid, double min_prop, double max_prop, double tau, Index n_substeps)
{
        cout << "DIAGNOSTICS" << endl;
        cout << "-----------" << endl;
        cout << "Memory requirement: "
             << 8.0 * grid.dx1 * grid.r / 1.0e9
             << " GB (X1), "
             << 8.0 * grid.dx2 * grid.r / 1.0e9
             << " GB (X2)" << endl;
        cout << "Min, max propensity: " << min_prop << ", " << max_prop << endl;
        cout << "Effective time step size: " << tau / n_substeps << endl;
        cout << "-----------" << endl;
        cout << endl;
}