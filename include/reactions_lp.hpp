#ifndef REACTIONS_HPP
#define REACTIONS_HPP

#include "reaction_class.hpp"

const std::vector<std::string> kNN = {"S1", "S2", "S3", "S4", "S5"};
constexpr double kA0 = 0.5;
constexpr double kA1 = 1.0;
constexpr double kA2 = 0.15;
constexpr double kA3 = 0.3;
constexpr double kA4 = 0.3;
constexpr double kB0 = 0.12;
constexpr double kB1 = 0.6;
constexpr double kB2 = 1.0;
constexpr double kB3 = 1.0;
constexpr double kB4 = 1.0;
constexpr double kC0 = 0.0025;
constexpr double kC1 = 0.0007;
constexpr double kC2 = 0.0231;
constexpr double kC3 = 0.01;
constexpr double kC4 = 0.01;

mysys mysystem(kNN);

myreact myreact0(
    {1, 0, 0, 0, 0}, {1}, [](std::vector<double> y)
    { return kA0 * kB0 / (kB0 + y[1]); },
    mysystem);

myreact myreact1(
    {0, 1, 0, 0, 0}, {0, 4}, [](std::vector<double> y)
    { return (kA1 + y[4]) * kB1 / (kB1 + y[0]); },
    mysystem);

myreact myreact2(
    {0, 0, 1, 0, 0}, {1}, [](std::vector<double> y)
    { return kA2 * kB2 * y[1] / (kB2 * y[1] + 1.0); },
    mysystem);

myreact myreact3(
    {0, 0, 0, 1, 0}, {2}, [](std::vector<double> y)
    { return kA3 * kB3 * y[2] / (kB3 * y[2] + 1.0); },
    mysystem);

myreact myreact4(
    {0, 0, 0, 0, 1}, {2}, [](std::vector<double> y)
    { return kA4 * kB4 * y[2] / (kB4 * y[2] + 1.0); },
    mysystem);

myreact myreact5(
    {-1, 0, 0, 0, 0}, {0}, [](std::vector<double> y)
    { return kC0 * y[0]; },
    mysystem);

myreact myreact6(
    {0, -1, 0, 0, 0}, {1}, [](std::vector<double> y)
    { return kC1 * y[1]; },
    mysystem);

myreact myreact7(
    {0, 0, -1, 0, 0}, {2}, [](std::vector<double> y)
    { return kC2 * y[2]; },
    mysystem);

myreact myreact8(
    {0, 0, 0, -1, 0}, {3}, [](std::vector<double> y)
    { return kC3 * y[3]; },
    mysystem);

myreact myreact9(
    {0, 0, 0, 0, -1}, {4}, [](std::vector<double> y)
    { return kC4 * y[4]; },
    mysystem);

#endif