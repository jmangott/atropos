#include "timer_class.hpp"

namespace get_time {
std::map<std::string, ntime> timers;

bool is_master()
{
#ifdef _OPENMP
    if (omp_get_thread_num() != 0)
        return false;
#endif
    return true;
}

void reset()
{
    for (auto& el : timers)
        el.second.reset();
}

void print()
{
    for (auto el : timers)
        std::cout << "get_time " << el.first << ": " << el.second.total() << " s"
                  << std::endl;
}

std::string sorted_output()
{
    typedef std::pair<std::string, ntime> pair_nt;
    auto comp = [](pair_nt a1, pair_nt a2) {
        return a1.second.total() > a2.second.total();
    };
    std::set<pair_nt, decltype(comp)> sorted(begin(timers), end(timers), comp);

    std::stringstream ss;
    ss.precision(4);
    ss.setf(std::ios_base::scientific);
    for (auto el : sorted) {
        ntime& t = el.second;
        ss << std::setw(40) << el.first << std::setw(15) << t.total() << std::setw(15)
           << t.count() << std::setw(15) << t.average() << std::setw(15)
           << t.deviation() / t.average() << std::endl;
    }
    return ss.str();
}

void start(std::string name)
{
    if (is_master())
        timers[name].start();
}

void stop(std::string name)
{
    if (is_master())
        timers[name].stop();
}

double total(std::string name) { return timers[name].total(); }

double average(std::string name) { return timers[name].average(); }

double deviation(std::string name) { return timers[name].deviation(); }
} // namespace get_time
