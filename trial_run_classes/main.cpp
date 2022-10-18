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

class myreact;

class mysys
{
public:
    vector<string> species_names;
    mysys(vector<string> _names);
    vector<myreact*> reactions;
    size_t mu();
};

mysys::mysys(vector<string> _names) : species_names(_names) {}

size_t mysys::mu()
{
    return reactions.size();
}

class myreact
{
public:
    vector<int> nu;
    double propensity(vector<double> xx);
    // myreact(mysys& _ref_system, vector<int> _nu, std::function<double (vector<double>)> _prop_function);
    myreact(mysys& _ref_system, vector<int> _nu, double (*_prop_function)(vector<double>));
    virtual ~myreact();
private:
    std::function<double (vector<double>)> prop_function;
    mysys& ref_system;
};

// myreact::myreact(mysys& _ref_system, vector<int> _nu, std::function<double (vector<double>)> _prop_function) : ref_system(_ref_system), nu(_nu), prop_function(_prop_function)
// {
//     ref_system.reactions.push_back(this);
// }

myreact::myreact(mysys& _ref_system, vector<int> _nu, double (*_prop_function)(vector<double>)) : ref_system(_ref_system), nu(_nu), prop_function(_prop_function)
{
    ref_system.reactions.push_back(this);
}

myreact::~myreact()
{
    ref_system.reactions.erase(std::remove(ref_system.reactions.begin(), ref_system.reactions.end(), this), ref_system.reactions.end());
}

double myreact::propensity(vector<double> xx)
{
    return prop_function(xx);
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

    mysys mysys1(nn);
    myreact myreact10(mysys1, nunu, prop0);
    // react1 myreact11(mysys1, nunu, [](vector<double> y){ return y[0] * y[1]; });

    system_ep mysys2;

    int64_t d1 = 0;
    int64_t d2 = 0;

    for (size_t ii = 0; ii < runs; ii++)
    {
        auto start1 = high_resolution_clock::now();
        for (size_t i = 0; i < max_it; i++)
        {
            xx[0] = i;
            mysys1.reactions[0]->propensity(xx);
        }
        auto stop1 = high_resolution_clock::now();
        auto duration1 = duration_cast<microseconds>(stop1 - start1);
        d1 += duration1.count();

        auto start2 = high_resolution_clock::now();
        for (size_t i = 0; i < max_it; i++)
        {
            xx[0] = i;
            mysys2.rule<0>(xx);
        }
        auto stop2 = high_resolution_clock::now();
        auto duration2 = duration_cast<microseconds>(stop2 - start2);
        d2 += duration2.count();

        // for (auto it : mysys1.reactions)
        // {
        //     cout << it->propensity() << endl;
        // }

    }

    cout << "Time taken by my class: " << d1 << " microseconds." << endl;
    cout << "Time taken by system_ep: " << d2 << " microseconds." << endl;

    return 0;
}