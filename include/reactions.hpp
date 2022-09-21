#ifndef REACTIONS_HPP
#define REACTIONS_HPP

#include "reaction_class.hpp"

std::vector<double> xx;
std::vector<std::string> nn;

mysys mysystem(xx, nn);

myreact myreact0(
    {-1, 0}, [](std::vector<double> y)
    { return y[0]; },
    mysystem);
myreact myreact1(
    {0, -1}, [](std::vector<double> y)
    { return y[1]; },
    mysystem);
myreact myreact2(
    {1, 0}, [](std::vector<double> y)
    { return 1.0 / (1.0 + y[1]); },
    mysystem);
myreact myreact3(
    {0, 1}, [](std::vector<double> y)
    { return 1.0 / (1.0 + y[0]); },
    mysystem);

#endif