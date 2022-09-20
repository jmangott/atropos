#ifndef REACTION_CLASS_HPP
#define REACTION_CLASS_HPP

#include <algorithm>
#include <string>
#include <vector>

class myreact;

class mysys
{
public:
    std::vector<double> &x;
    std::vector<std::string> species_names;
    mysys(std::vector<double> &_x, std::vector<std::string> _names);
    std::vector<myreact *> reactions;
    size_t mu();
};

class myreact
{
public:
    std::vector<int> nu;
    double propensity();
    myreact(std::vector<int> _nu, double (*_prop_function)(std::vector<double>), mysys &_ref_system);
    virtual ~myreact();

private:
    double (*prop_function)(std::vector<double>);
    mysys &ref_system;
};

#endif