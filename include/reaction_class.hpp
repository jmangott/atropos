#ifndef REACTION_CLASS_HPP
#define REACTION_CLASS_HPP

#include <algorithm>
#include <string>
#include <vector>

// TODO: replace vector `reactions` in mysys by a map and include in myreact a public variable `name`, which acts as key for the `reactions` map

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
    std::vector<size_t> depends_on; 
    double propensity();
    double propensity(std::vector<double> x);
    myreact(std::vector<int> _nu, std::vector<size_t> _depends_on, double (*_prop_function)(std::vector<double>), mysys &_ref_system);
    virtual ~myreact();

private:
    double (*prop_function)(std::vector<double>);
    mysys &ref_system;
};

#endif