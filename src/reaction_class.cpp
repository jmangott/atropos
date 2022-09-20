#include "reaction_class.hpp"

using std::string;
using std::vector;

mysys::mysys(vector<double> &_x, vector<string> _names) : x(_x), species_names(_names) {}

size_t mysys::mu()
{
    return reactions.size();
}

myreact::myreact(vector<int> _nu, double (*_prop_function)(vector<double>), mysys &_ref_system) : nu(_nu), prop_function(_prop_function), ref_system(_ref_system)
{
    ref_system.reactions.push_back(this);
}

myreact::~myreact()
{
    ref_system.reactions.erase(std::remove(ref_system.reactions.begin(), ref_system.reactions.end(), this), ref_system.reactions.end());
}

double myreact::propensity()
{
    return prop_function(ref_system.x);
}