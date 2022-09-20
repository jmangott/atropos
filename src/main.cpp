#include <iostream>
#include <string>
#include <vector>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>
#include <lr/coefficients.hpp>
#include <lr/lr.hpp>

#include "reactions.hpp"

using std::cout;
using std::endl;
using std::string;
using std::vector;

int main()
{
    vector<double> xx = {1, 2, 3, 4};
    vector<string> nn = {"a", "b", "c", "d"};
    vector<int> nunu = {0, 0, 0, 0};
    mysys mysys1(xx, nn);
    myreact myreact1(nunu, [](vector<double> y) { return y[0] * y[1]; }, mysys1);
    multi_array<double, 1> my_mat({2});
    cout << "Hello!" << endl;

    return 0;
}