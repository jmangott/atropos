#ifndef REACTIONS_HPP
#define REACTIONS_HPP

#include "../reaction_class.hpp"

const std::vector<std::string> kNN = {"cyclin", "pcyclin_pcdc2", "cdc2", "pcdc2", "pcyclin_cdc2"};
constexpr double kNA_V = 6.02214e3;

constexpr double kA0 = 0.015 * kNA_V;
constexpr double kA1 = 200.0 / kNA_V;
constexpr double kA2 = 0.018;
constexpr double kA3 = 360.0 / kNA_V / kNA_V;
constexpr double kA4 = 1.0;
constexpr double kA5 = 100;
constexpr double kA6 = 10;

mysys mysystem(kNN);

myreact myreact0(
    {1, 0, 0, 0, 0}, {}, [](std::vector<double> y)
    { return kA0; },
    mysystem);
myreact myreact1(
    {-1, 1, 0, -1, 0}, {0, 3}, [](std::vector<double> y)
    { return kA1 * y[0] * y[3]; },
    mysystem);
myreact myreact2(
    {0, -1, 0, 0, 1}, {1}, [](std::vector<double> y)
    { return kA2 * y[1]; },
    mysystem);
myreact myreact3(
    {0, -1, 0, 0, 1}, {1, 4}, [](std::vector<double> y)
    { return kA3 * y[1] * y[4] * (y[4] - 1.0) / 2.0; },
    mysystem);
myreact myreact4(
    {0, 0, 1, 0, -1}, {4}, [](std::vector<double> y)
    { return kA4 * y[4]; },
    mysystem);
myreact myreact5(
    {0, 0, -1, 1, 0}, {2}, [](std::vector<double> y)
    { return kA5 * y[2]; },
    mysystem);
myreact myreact6(
    {0, 0, 1, -1, 0}, {3}, [](std::vector<double> y)
    { return kA6 * y[3]; },
    mysystem);

#endif