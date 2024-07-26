#include "matrix.hpp"

template <>
void Matrix::Matricize<0, 3, double>(const multi_array<double, 3> &input, multi_array<double, 2> &output)
{
    const auto n = input.shape();
    for (Index k = 0; k < n[2]; ++k)
    {
        for (Index j = 0; j < n[1]; ++j)
        {
            for (Index i = 0; i < n[0]; ++i)
            {
                output(j + n[1] * k, i) = input(i, j, k);
            }
        }
    }
}

template <>
void Matrix::Matricize<1, 3, double>(const multi_array<double, 3> &input, multi_array<double, 2> &output)
{
    const auto n = input.shape();
    for (Index k = 0; k < n[2]; ++k)
    {
        for (Index j = 0; j < n[1]; ++j)
        {
            for (Index i = 0; i < n[0]; ++i)
            {
                output(k + n[2] * i, j) = input(i, j, k);
            }
        }
    }
}

template <>
void Matrix::Matricize<2, 3, double>(const multi_array<double, 3> &input, multi_array<double, 2> &output)
{
    const auto n = input.shape();
    for (Index k = 0; k < n[2]; ++k)
    {
        for (Index j = 0; j < n[1]; ++j)
        {
            for (Index i = 0; i < n[0]; ++i)
            {
                output(i + n[0] * j, k) = input(i, j, k);
            }
        }
    }
}

template <>
void Matrix::Tensorize<0, 3, double>(const multi_array<double, 2> &input, multi_array<double, 3> &output)
{
    const auto n = output.shape();
    for (Index k = 0; k < n[2]; ++k)
    {
        for (Index j = 0; j < n[1]; ++j)
        {
            for (Index i = 0; i < n[0]; ++i)
            {
                output(i, j, k) = input(j + n[1] * k, i);
            }
        }
    }
}

template <>
void Matrix::Tensorize<1, 3, double>(const multi_array<double, 2> &input, multi_array<double, 3> &output)
{
    const auto n = output.shape();
    for (Index k = 0; k < n[2]; ++k)
    {
        for (Index j = 0; j < n[1]; ++j)
        {
            for (Index i = 0; i < n[0]; ++i)
            {
                output(i, j, k) = input(k + n[2] * i, j);
            }
        }
    }
}

template <>
void Matrix::Tensorize<2, 3, double>(const multi_array<double, 2> &input, multi_array<double, 3> &output)
{
    const auto n = output.shape();
    for (Index k = 0; k < n[2]; ++k)
    {
        for (Index j = 0; j < n[1]; ++j)
        {
            for (Index i = 0; i < n[0]; ++i)
            {
                output(i, j, k) = input(i + n[0] * j, k);
            }
        }
    }
}

// Shift operator
template <>
void Matrix::ShiftRows<1>(multi_array<double, 2> &output_array, const multi_array<double, 2> &input_array, const grid_parms grid, const Index mu);

// Inverse shift operator
template <>
void Matrix::ShiftRows<-1>(multi_array<double, 2> &output_array, const multi_array<double, 2> &input_array, const grid_parms grid, const Index mu);