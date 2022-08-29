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

typedef long ind;

//////////////////////////
// LOW_RANK_CME OBJECTS //
//////////////////////////

#define RULE_SET(__name, __d, ...)            \
    struct __name                             \
    {                                         \
        template <ind i>                      \
        static double rule(vector<double> x); \
                                              \
        template <ind i>                      \
        vector<ind> static depends_on();      \
                                              \
        static constexpr size_t d = (__d);    \
                                              \
        static vector<string> names()         \
        {                                     \
            return {__VA_ARGS__};             \
        }                                     \
    }

#define RULE(__name, __i, __rule, ...)              \
    template <>                                     \
    double __name::rule<(__i)>(vector<__name::d> x) \
    {                                               \
        return (__rule);                            \
    }                                               \
    template <>                                     \
    vector<ind> __name::depends_on<(__i)>()         \
    {                                               \
        return {__VA_ARGS__};                       \
    }

struct system_ep
{
    template <ind i>
    static double rule(vector<double> y);

    template <ind i>
    vector<ind> static depends_on();

    static vector<std::string> names()
    {
        return {"a", "b", "c", "d"};
    }
};

template <>
double system_ep::rule<0>(vector<double> y)
{
    return y[0] * y[1];
}

////////////////
// MY CLASSES //
////////////////

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


////////////////////////////////////////////
// COMPARISON OF THE DIFFERENT APPROACHES //
////////////////////////////////////////////

int main()
{
    int runs = 20;
    int max_it = 100000;


    vector<double> xx = {1, 2, 3, 4};
    vector<string> nn = {"a", "b", "c", "d"};
    vector<int> nunu = {0, 0, 0, 0};

    sys mysys(xx, nn);
    react myreact0(nunu, prop0);
    // react myreact1(nunu, [](vector<double> y){return y[0] * y[1];});

    sys1 mysys1(xx, nn);
    react1 myreact10(mysys1, nunu, prop0);
    // react1 myreact11(mysys1, nunu, [](vector<double> y){ return y[0] * y[1]; });

    system_ep mysys2;

    int64_t d1 = 0;
    int64_t d2 = 0;
    int64_t d3 = 0;

    for (size_t ii = 0; ii < runs; ii++)
    {
        auto start1 = high_resolution_clock::now();
        for (size_t i = 0; i < max_it; i++)
        {
            xx[0] = i;
            mysys.propensity(myreact0);
        }
        auto stop1 = high_resolution_clock::now();
        auto duration1 = duration_cast<microseconds>(stop1 - start1);
        d1 += duration1.count();

        auto start2 = high_resolution_clock::now();
        for (size_t i = 0; i < max_it; i++)
        {
            xx[0] = i;
            mysys1.reactions[0]->propensity();
        }
        auto stop2 = high_resolution_clock::now();
        auto duration2 = duration_cast<microseconds>(stop2 - start2);
        d2 += duration2.count();

        auto start3 = high_resolution_clock::now();
        for (size_t i = 0; i < max_it; i++)
        {
            xx[0] = i;
            mysys2.rule<0>(xx);
        }
        auto stop3 = high_resolution_clock::now();
        auto duration3 = duration_cast<microseconds>(stop3 - start3);
        d3 += duration3.count();

        // for (auto it : mysys1.reactions)
        // {
        //     cout << it->propensity() << endl;
        // }

    }

    cout << "Time taken by functor class: " << d1 << " microseconds." << endl;
    cout << "Time taken by my class: " << d2 << " microseconds." << endl;
    cout << "Time taken by system_ep: " << d3 << " microseconds." << endl;

    return 0;
}