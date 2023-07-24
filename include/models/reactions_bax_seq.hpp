#ifndef REACTIONS_HPP
#define REACTIONS_HPP

#include "../reaction_class.hpp"

const std::vector<std::string> kNN = {"Bax6", "Bax4mSmac", "Bax5mSmac", "Bax6mSmac", "mSmac", "cSmac", "Bax1", "Bax2", "Bax3", "Bax4", "Bax5"};

constexpr double kA0 = 0.0002;
constexpr double kA1 = 0.001;
constexpr double kA2 = 0.0002;
constexpr double kA3 = 0.001;
constexpr double kA4 = 0.0002;
constexpr double kA5 = 0.001;
constexpr double kA6 = 0.0002;
constexpr double kA7 = 0.001;
constexpr double kA8 = 0.0002;
constexpr double kA9 = 0.001;

constexpr double kA10 = 3.0e-5;
constexpr double kA11 = 0.001;
constexpr double kA12 = 10.0;
constexpr double kA13 = 3.0e-5;
constexpr double kA14 = 0.001;
constexpr double kA15 = 10.0;
constexpr double kA16 = 3.0e-5;
constexpr double kA17 = 0.001;
constexpr double kA18 = 10.0;

mysys mysystem(kNN);

myreact myreact0(
    {0, 0, 0, 0, 0, 0, -2, 1, 0, 0, 0}, {6}, [](std::vector<double> y)
    { return kA0 * y[6] * (y[6] - 1.0) / 2.0; },
    mysystem);

myreact myreact1(
    {0, 0, 0, 0, 0, 0, 2, -1, 0, 0, 0}, {7}, [](std::vector<double> y)
    { return kA1 * y[7]; },
    mysystem);

myreact myreact2(
    {0, 0, 0, 0, 0, 0, -1, -1, 1, 0, 0}, {6, 7}, [](std::vector<double> y)
    { return kA2 * y[6] * y[7]; },
    mysystem);

myreact myreact3(
    {0, 0, 0, 0, 0, 0, 1, 1, -1, 0, 0}, {8}, [](std::vector<double> y)
    { return kA3 * y[8]; },
    mysystem);

myreact myreact4(
    {0, 0, 0, 0, 0, 0, -1, 0, -1, 1, 0}, {6, 8}, [](std::vector<double> y)
    { return kA4 * y[6] * y[8]; },
    mysystem);

myreact myreact5(
    {0, 0, 0, 0, 0, 0, 1, 0, 1, -1, 0}, {9}, [](std::vector<double> y)
    { return kA5 * y[9]; },
    mysystem);

myreact myreact6(
    {0, 0, 0, 0, 0, 0, -1, 0, 0, -1, 1}, {6, 9}, [](std::vector<double> y)
    { return kA6 * y[6] * y[9]; },
    mysystem);

myreact myreact7(
    {0, 0, 0, 0, 0, 0, 1, 0, 0, 1, -1}, {10}, [](std::vector<double> y)
    { return kA7 * y[10]; },
    mysystem);

myreact myreact8(
    {1, 0, 0, 0, 0, 0, -1, 0, 0, 0, -1}, {6, 10}, [](std::vector<double> y)
    { return kA8 * y[6] * y[10]; },
    mysystem);

myreact myreact9(
    {-1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1}, {0}, [](std::vector<double> y)
    { return kA9 * y[0]; },
    mysystem);

myreact myreact10(
    {0, 1, 0, 0, -1, 0, 0, 0, 0, -1, 0}, {9, 4}, [](std::vector<double> y)
    { return kA10 * y[9] * y[4]; },
    mysystem);

myreact myreact11(
    {0, -1, 0, 0, 1, 0, 0, 0, 0, 1, 0}, {1}, [](std::vector<double> y)
    { return kA11 * y[1]; },
    mysystem);

myreact myreact12(
    {0, -1, 0, 0, 0, 1, 0, 0, 0, 1, 0}, {1}, [](std::vector<double> y)
    { return kA12 * y[1]; },
    mysystem);

myreact myreact13(
    {0, 0, 1, 0, -1, 0, 0, 0, 0, 0, -1}, {10, 4}, [](std::vector<double> y)
    { return kA13 * y[10] * y[4]; },
    mysystem);

myreact myreact14(
    {0, 0, -1, 0, 1, 0, 0, 0, 0, 0, 1}, {2}, [](std::vector<double> y)
    { return kA14 * y[2]; },
    mysystem);

myreact myreact15(
    {0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 1}, {2}, [](std::vector<double> y)
    { return kA15 * y[2]; },
    mysystem);

myreact myreact16(
    {-1, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0}, {0, 4}, [](std::vector<double> y)
    { return kA16 * y[0] * y[4]; },
    mysystem);

myreact myreact17(
    {1, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0}, {3}, [](std::vector<double> y)
    { return kA17 * y[3]; },
    mysystem);

myreact myreact18(
    {1, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0}, {3}, [](std::vector<double> y)
    { return kA18 * y[3]; },
    mysystem);

#endif