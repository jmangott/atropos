#ifndef REACTION_CLASS_HPP
#define REACTION_CLASS_HPP

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

typedef ptrdiff_t Index;

class myreact;

class mysys
{
public:
    std::vector<std::string> species_names;
    mysys(std::vector<std::string> _names);
    std::vector<myreact *> reactions;
    std::vector<myreact *>::size_type mu();
};

class myreact
{
public:
    std::vector<int> nu;
    std::vector<Index> depends_on; 
    double propensity(std::vector<double> x);
    myreact(std::vector<int> _nu, std::vector<Index> _depends_on, double (*_prop_function)(std::vector<double>), mysys &_ref_system);
    virtual ~myreact();

private:
    double (*prop_function)(std::vector<double>);
    mysys &ref_system;
};

#endif