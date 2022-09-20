#ifndef REACTIONS_HPP
#define REACTIONS_HPP

#include "reaction_class.hpp"

std::vector<double> xx;
std::vector<std::string> nn;

mysys mysys1(xx, nn);

myreact myreact1(
    {-1, 0}, [](std::vector<double> y)
    { return y[0]; },
    mysys1);
myreact myreact2(
    {0, -1}, [](std::vector<double> y)
    { return y[1]; },
    mysys1);
myreact myreact3(
    {1, 0}, [](std::vector<double> y)
    { return 1.0 / (1.0 + y[1]); },
    mysys1);
myreact myreact4(
    {0, 1}, [](std::vector<double> y)
    { return 1.0 / (1.0 + y[0]); },
    mysys1);

#endif