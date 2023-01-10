#ifndef REACTIONS_HPP
#define REACTIONS_HPP

#include "reaction_class.hpp"

const std::vector<std::string> kNN = {"S1", "S2"};
constexpr double kB = 0.4;
constexpr double kC = 0.05;

mysys mysystem(kNN);

myreact myreact0(
    {-1, 0}, {0}, [](std::vector<double> y)
    { return kC * y[0]; },
    mysystem);
myreact myreact1(
    {0, -1}, {1}, [](std::vector<double> y)
    { return kC * y[1]; },
    mysystem);
myreact myreact2(
    {1, 0}, {1}, [](std::vector<double> y)
    { return kB / (kB + y[1]); },
    mysystem);
myreact myreact3(
    {0, 1}, {0}, [](std::vector<double> y)
    { return kB / (kB + y[0]); },
    mysystem);

#endif