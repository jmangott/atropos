#ifndef REACTIONS_HPP
#define REACTIONS_HPP

#include "../reaction_class.hpp"

const std::vector<std::string> kNN = {"S0", "S1", "S2", "S3", "S4", "S5"};
constexpr double kp1 = 0.4;
constexpr double kp2 = 100.0;
constexpr double kp3 = 100.0;
constexpr double km1 = 2.0;
constexpr double km2 = 1.0;
constexpr double km3 = 50.0;


mysys mysystem(kNN);

myreact myreact0(
    {-1, -1, 1, 0, 0, 0}, {0, 1}, [](std::vector<double> y)
    { return kp1 * y[0] * y[1]; },
    mysystem);

myreact myreact1(
    {1, 1, -1, 0, 0, 0}, {2}, [](std::vector<double> y)
    { return kp2 * y[2]; },
    mysystem);

myreact myreact2(
    {0, 0, 0, -1, -1, 1}, {3, 4}, [](std::vector<double> y)
    { return km1 * y[3] * y[4]; },
    mysystem);

myreact myreact3(
    {0, 0, 0, 1, 1, -1}, {5}, [](std::vector<double> y)
    { return km2 * y[5]; },
    mysystem);

myreact myreact4(
    {0, 1, -1, 1, 0, 0}, {2}, [](std::vector<double> y)
    { return kp3 * y[2]; },
    mysystem);

myreact myreact5(
    {1, 0, 0, 0, 1, -1}, {5}, [](std::vector<double> y)
    { return km3 * y[5]; },
    mysystem);

#endif