# DLR approximation for the kinetic CME

## Requirements
- CMake (3.22.1 or later)
- C++17 compatible C++ compiler
- Fortran compiler (if OpenBLAS is used)
- HDF5 (1.10.x)
- netCDF-4

Optionally:
- OpenMP
- Intel MKL

Check via `nc-config --has-hdf5`, whether HDF5 was used in the netCDF-4 build.

## Installation
Clone the repository via
```shell
git clone https://git.uibk.ac.at/c7021158/kinetic-cme
```
and build the program by executing
```shell
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

The generated executable `kinetic-cme` can be found in `bin`.

### Intel MKL
If you prefer to use Intel MKL as the BLAS and LAPACK backend instead of OpenBLAS set 
```shell
export MKLROOT=/path/to/intel/mkl
cmake -DCMAKE_BUILD_TYPE=Release -DMKL_ENABLED=ON ..
make
```
and make sure to add the MKL libraries to your `LD_LIBRARY_PATH`, i.e.
```shell
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/intel/mkl/lib/intel64_lin/
```
before running the executable.

### OpenMP
OpenMP can be activated via
```shell
cmake -DCMAKE_BUILD_TYPE=Release -DOPENMP=ON ..
```
MacOS: Note that XCode compilers do not support OpenMP. For using OpenMP on macOS, a manual installation (e.g. of `gcc11`) is required and the `CXX`, `CC` and `FC` environment variables have to be set accordingly.

## Preparing input data
A python script called `initial_condition.py` is provided in the `input` folder to generate input data.
This is still work in progress.

