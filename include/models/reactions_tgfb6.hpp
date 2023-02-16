#ifndef REACTIONS_HPP
#define REACTIONS_HPP

#include "../reaction_class.hpp"

const std::vector<std::string> kNN = {"S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8"};
constexpr double kA0 = 5.525363366159611461e-02;
constexpr double kA1 = 0.1;
constexpr double kA2 = 1.421915991687823499e+01;
constexpr double kA3 = 2.080559410140531806e+00;
constexpr double kA4 = 1.975939657536690888e-02;
constexpr double kA5 = 8.350953931031456223e-01;
constexpr double kA6 = 0.1;
constexpr double kA7 = 5.022274319604386195e-01;
constexpr double kA8 = 8.718972334813199752e-05;
constexpr double kA9 = 2.915856647494313378e-02;
constexpr double kA10 = 1.339559213470352894e+01;
constexpr double kA11 = 1.151177584786399599e+01;
constexpr double kA12 = 2.051775933975554054e-03;
constexpr double kA13 = 2.182132701358179272e+00;
constexpr double kA14 = 7.051867558145050729e-01;
constexpr double kA15 = 1.565395682931424934e-02;
constexpr double kA16 = 1.108779807581306898e+01;
constexpr double kA17 = 3.491475235854992198e-01;

// Check if these constants are correct
constexpr double kX0 = 0.0;
constexpr double kX1 = 1.0;
constexpr double kX2 = 0.0;
constexpr double kX3 = 0.0;
constexpr double kX4 = 0.0;

mysys mysystem(kNN);

myreact myreact0(
    {-1, 1, 0, 0, 0, 0, 0, 0}, {0, 3, 5}, [](std::vector<double> y)
    { return kA14 * kX0 * y[0] * y[5] / ((kA15 + y[0]) * (1.0 + y[3] / kA16)); },
    mysystem);

myreact myreact1(
    {1, -1, 0, 0, 0, 0, 0, 0}, {1}, [](std::vector<double> y)
    { return kA0 * y[1]; },
    mysystem);

myreact myreact2(
    {0, 0, 0, 0, 0, 0, -1, 1}, {6}, [](std::vector<double> y)
    { return kA1 * kX1 * y[6]; },
    mysystem);

myreact myreact3(
    {0, 0, 0, 0, 0, 0, -1, 1}, {5, 6}, [](std::vector<double> y)
    { return kA17 * y[5] * y[6]; },
    mysystem);

myreact myreact4(
    {0, 0, 0, 0, 0, 0, -1, 1}, {6}, [](std::vector<double> y)
    { return kA3 * kX0 * y[6] / (kA4 + y[6] + kA4 * kX2 / kA5); },
    mysystem);

myreact myreact5(
    {0, 0, 0, 0, 0, 0, 1, -1}, {7}, [](std::vector<double> y)
    { return kA2 * y[7]; },
    mysystem);

myreact myreact6(
    {0, 0, 0, 0, -1, 1, 0, 0}, {4}, [](std::vector<double> y)
    { return kA6 * kX1 * y[4]; },
    mysystem);

myreact myreact7(
    {0, 0, 0, 0, 1, -1, 0, 0}, {3, 5}, [](std::vector<double> y)
    { return kA7 * y[3] * y[5]; },
    mysystem);

myreact myreact8(
    {0, 0, 0, 0, -1, 1, 0, 0}, {4}, [](std::vector<double> y)
    { return kA8 * kX0 * y[4] / ((kA9 + y[4]) * (1.0 + kX3 / kA9)); },
    mysystem);

myreact myreact9(
    {0, 0, -1, 1, 0, 0, 0, 0}, {2, 5}, [](std::vector<double> y)
    { return kA10 * y[2] * y[5] / (kA11 + y[2] + kA11 * kX4 / kA12); },
    mysystem);

myreact myreact10(
    {0, 0, 1, -1, 0, 0, 0, 0}, {3}, [](std::vector<double> y)
    { return kA13 * y[3]; },
    mysystem);

// myreact myreact0(
//     {-1, 1, 0, 0, 0, 0, 0, 0}, {0, 3, 5}, [](std::vector<double> y)
//     { return kA14 * kX0 * y[0] * y[5] / ((kA15 + y[0]) * (1.0 + y[3] / kA16)); },
//     mysystem);

// myreact myreact1(
//     {1, -1, 0, 0, 0, 0, 0, 0}, {1}, [](std::vector<double> y)
//     { return kA0 * y[1]; },
//     mysystem);

// myreact myreact2(
//     {0, 0, -1, 1, 0, 0, 0, 0}, {2}, [](std::vector<double> y)
//     { return kA1 * kX1 * y[2]; },
//     mysystem);

// myreact myreact3(
//     {0, 0, -1, 1, 0, 0, 0, 0}, {2, 5}, [](std::vector<double> y)
//     { return kA17 * y[2] * y[5]; },
//     mysystem);

// myreact myreact4(
//     {0, 0, -1, 1, 0, 0, 0, 0}, {2}, [](std::vector<double> y)
//     { return kA3 * kX0 * y[2] / (kA4 + y[2] + kA4 * kX2 / kA5); },
//     mysystem);

// myreact myreact5(
//     {0, 0, 1, -1, 0, 0, 0, 0}, {3}, [](std::vector<double> y)
//     { return kA2 * y[3]; },
//     mysystem);

// myreact myreact6(
//     {0, 0, 0, 0, -1, 1, 0, 0}, {4}, [](std::vector<double> y)
//     { return kA6 * kX1 * y[4]; },
//     mysystem);

// myreact myreact7(
//     {0, 0, 0, 0, 1, -1, 0, 0}, {5, 7}, [](std::vector<double> y)
//     { return kA7 * y[5] * y[7]; },
//     mysystem);

// myreact myreact8(
//     {0, 0, 0, 0, -1, 1, 0, 0}, {4}, [](std::vector<double> y)
//     { return kA8 * kX0 * y[4] / ((kA9 + y[4]) * (1.0 + kX3 / kA9)); },
//     mysystem);

// myreact myreact9(
//     {0, 0, 0, 0, 0, 0, -1, 1}, {5, 6}, [](std::vector<double> y)
//     { return kA10 * y[5] * y[6] / (kA11 + y[6] + kA11 * kX4 / kA12); },
//     mysystem);

// myreact myreact10(
//     {0, 0, 0, 0, 0, 0, 1, -1}, {7}, [](std::vector<double> y)
//     { return kA13 * y[7]; },
//     mysystem);

#endif