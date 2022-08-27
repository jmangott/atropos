#include <algorithm>
#include <chrono>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

using namespace std::chrono;
using std::cout;
using std::endl;
using std::function;
using std::string;
using std::vector;

class sys
{
public:
    vector<double> &x;
    vector<string> species_names;
    sys(vector<double> &_x, vector<string> _names);
    // static set<react *> react_instances;

    template <class T>
    double propensity(T &prop_functor)
    {
        return prop_functor(x);
    }
};

sys::sys(vector<double> &_x, vector<string> _names) : x(_x), species_names(_names) {}

class react
{
public:
    vector<int> population_change;
    // double propensity();
    react(vector<int> _nu, double (*_prop)(vector<double>)) : population_change(_nu), prop(_prop)
    {
        // reference_sys.react_instances.insert(this);
    }
    virtual ~react()
    {
        // reference_sys.react_instances.erase(this);
    }
    double operator()(const vector<double> x)
    {
        return prop(x);
    }

private:
    double (*prop)(vector<double>);
};

/////////////////////////
// ALTERNATIVE CLASSES //
/////////////////////////

class react1;

class sys1
{
public:
    vector<double> &x;
    vector<string> species_names;
    sys1(vector<double> &_x, vector<string> _names);
    vector<react1*> reactions;
    size_t mu();
};

sys1::sys1(vector<double> &_x, vector<string> _names) : x(_x), species_names(_names) {}

size_t sys1::mu()
{
    return reactions.size();
}

class react1
{
public:
    vector<int> nu;
    double propensity();
    react1(sys1& _ref_system, vector<int> _nu, double (*_prop_function)(vector<double>));
    virtual ~react1();
private:
    double (*prop_function)(vector<double>);
    sys1& ref_system;
};

react1::react1(sys1& _ref_system, vector<int> _nu, double (*_prop_function)(vector<double>)) : ref_system(_ref_system), nu(_nu), prop_function(_prop_function)
{
    ref_system.reactions.push_back(this);
}

react1::~react1()
{
    ref_system.reactions.erase(std::remove(ref_system.reactions.begin(), ref_system.reactions.end(), this), ref_system.reactions.end());
}

double react1::propensity()
{
    return prop_function(ref_system.x);
}

//////////////////
// PROPENSITIES //
//////////////////

double prop0(vector<double> y)
{
    return y[0] * y[1];
}

int main()
{
    vector<double> xx = {1, 2, 3, 4};
    vector<string> nn = {"a", "b", "c", "d"};
    vector<int> nunu = {0, 0, 0, 0};

    sys mysys(xx, nn);
    react myreact0(nunu, prop0);
    react myreact1(nunu, [](vector<double> y){return y[0] * y[1];});

    sys1 mysys1(xx, nn);
    react1 myreact10(mysys1, nunu, prop0);
    react1 myreact11(mysys1, nunu, [](vector<double> y){ return y[0] * y[1]; });

    auto start1 = high_resolution_clock::now();
    for (size_t i = 0; i < 1000; i++)
    {
        xx[0] = i;
        mysys.propensity(myreact0);
    }
    auto stop1 = high_resolution_clock::now();
    auto duration1 = duration_cast<microseconds>(stop1 - start1);

    cout << "Time taken by loop 1:" << duration1.count() << " microseconds" << endl;

    auto start2 = high_resolution_clock::now();
    for (size_t i = 0; i < 1000; i++)
    {
        xx[0] = i;
        mysys.propensity(myreact1);
    }
    auto stop2 = high_resolution_clock::now();
    auto duration2 = duration_cast<microseconds>(stop2 - start2);

    cout << "Time taken by loop 2:" << duration2.count() << " microseconds" << endl;

    auto start3 = high_resolution_clock::now();
    for (size_t i = 0; i < 1000; i++)
    {
        xx[0] = i;
        mysys1.reactions[1]->propensity();
    }
    auto stop3 = high_resolution_clock::now();
    auto duration3 = duration_cast<microseconds>(stop3 - start3);

    cout << "Time taken by loop 3:" << duration3.count() << " microseconds" << endl;

    // for (auto it : mysys1.reactions)
    // {
    //     cout << it->propensity() << endl;
    // }

    return 0;
}