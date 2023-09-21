#ifndef COEFF_CLASS_HPP
#define COEFF_CLASS_HPP

#include<memory>
#include <vector>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>

#include "grid_class.hpp"


struct coeff {};

struct root_coeff : coeff
{
    std::vector<multi_array<double, 2>> propensity;
};

struct internal_coeff : coeff
{
    multi_array<double, 6> g_coeff, h_coeff;
};

struct external_coeff : coeff
{
    std::vector<multi_array<double, 3>> c_coeff, d_coeff;
    multi_array<double, 4> e_coeff, f_coeff;
};


#endif