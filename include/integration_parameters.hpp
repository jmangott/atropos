#include <cstddef>

constexpr double kTstar = 10.0;             // final time
double Tau = 0.01;                          // time step size
constexpr std::ptrdiff_t kSnapshot = 100;   // number of time steps between snapshots
constexpr char kFilename[] = "ts";          // name of the output folder