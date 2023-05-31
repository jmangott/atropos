#ifndef REACTIONS_HPP
#define REACTIONS_HPP

#include "../reaction_class.hpp"

const std::vector<std::string> kNN = {"Bax5", "Bax6", "Bax4mSmac", "Bax5mSmac", "Bax6mSmac", "mSmac", "cSmac", "Bax1", "Bax2", "Bax3", "Bax4"};

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
    {0, 0, 0, 0, 0, 0, 0, -2, 1, 0, 0}, {7}, [](std::vector<double> y)
    { return kA0 * y[7] * (y[7] - 1.0) / 2.0; },
    mysystem);

myreact myreact1(
    {0, 0, 0, 0, 0, 0, 0, 2, -1, 0, 0}, {8}, [](std::vector<double> y)
    { return kA1 * y[8]; },
    mysystem);

myreact myreact2(
    {0, 0, 0, 0, 0, 0, 0, -1, -1, 1, 0}, {7, 8}, [](std::vector<double> y)
    { return kA2 * y[7] * y[8]; },
    mysystem);

myreact myreact3(
    {0, 0, 0, 0, 0, 0, 0, 1, 1, -1, 0}, {9}, [](std::vector<double> y)
    { return kA3 * y[9]; },
    mysystem);

myreact myreact4(
    {0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 1}, {7, 9}, [](std::vector<double> y)
    { return kA4 * y[7] * y[9]; },
    mysystem);

myreact myreact5(
    {0, 0, 0, 0, 0, 0, 0, 1, 0, 1, -1}, {10}, [](std::vector<double> y)
    { return kA5 * y[10]; },
    mysystem);

myreact myreact6(
    {1, 0, 0, 0, 0, 0, 0, -1, 0, 0, -1}, {7, 10}, [](std::vector<double> y)
    { return kA6 * y[7] * y[10]; },
    mysystem);

myreact myreact7(
    {-1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1}, {0}, [](std::vector<double> y)
    { return kA7 * y[0]; },
    mysystem);

myreact myreact8(
    {-1, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0}, {7, 0}, [](std::vector<double> y)
    { return kA8 * y[7] * y[0]; },
    mysystem);

myreact myreact9(
    {1, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0}, {1}, [](std::vector<double> y)
    { return kA9 * y[1]; },
    mysystem);

myreact myreact10(
    {0, 0, 1, 0, 0, -1, 0, 0, 0, 0, -1}, {10, 5}, [](std::vector<double> y)
    { return kA10 * y[10] * y[5]; },
    mysystem);

myreact myreact11(
    {0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 1}, {2}, [](std::vector<double> y)
    { return kA11 * y[2]; },
    mysystem);

myreact myreact12(
    {0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 1}, {2}, [](std::vector<double> y)
    { return kA12 * y[2]; },
    mysystem);

myreact myreact13(
    {-1, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0}, {0, 5}, [](std::vector<double> y)
    { return kA13 * y[0] * y[5]; },
    mysystem);

myreact myreact14(
    {1, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0}, {3}, [](std::vector<double> y)
    { return kA14 * y[3]; },
    mysystem);

myreact myreact15(
    {1, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0}, {3}, [](std::vector<double> y)
    { return kA15 * y[3]; },
    mysystem);

myreact myreact16(
    {0, -1, 0, 0, 1, -1, 0, 0, 0, 0, 0}, {1, 5}, [](std::vector<double> y)
    { return kA16 * y[1] * y[5]; },
    mysystem);

myreact myreact17(
    {0, 1, 0, 0, -1, 1, 0, 0, 0, 0, 0}, {4}, [](std::vector<double> y)
    { return kA17 * y[4]; },
    mysystem);

myreact myreact18(
    {0, 1, 0, 0, -1, 0, 1, 0, 0, 0, 0}, {4}, [](std::vector<double> y)
    { return kA18 * y[4]; },
    mysystem);

#endif