# DLR approximation for the kinetic CME

## Objective
Solves the chemical master equation with the projector-splitting based dynamical low-rank approximation.

## Requirements
- CMake (3.22.1 or later)
- C++17 compatible C++ compiler
- Fortran compiler (if OpenBLAS is used)
- HDF5 (1.10.x)
- netCDF4
- Python (>3.8)

Check via `nc-config --has-hdf5`, whether HDF5 was used in the netCDF-4 build.

Optionally:
- OpenMP
- Intel MKL
- StochKit2 (for generating reference solutions for the example problems)

## Installation
Clone the repository via
```shell
git clone https://git.uibk.ac.at/c7021158/kinetic-cme
```
and build the program by executing
```shell
cd kinetic-cme
cmake -B build -DCMAKE_BUILD_TYPE=Release ..
cmake --build build
```

The generated executable `kinetic-cme` can be found in `bin`.
To enable compiler options for debugging, use `-DCMAKE_BUILD_TYPE=Debug` instead.

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
A python script called `initial_condition.py` is provided in `scripts/input` to generate input data.

## Python environment
To use the Python scripts included in `scripts`, a Python environment with the packages specified in `pip-requirements.txt` needs to be configured and enabled.
For Python venv:
```shell
$ python -m venv path/to/my_venv
$ source path/to/my_venv/bin/activate
$ pip install -r pip-requirements.txt
$ pip install -e .
```
For anaconda:
```shell
$ conda create -n my_venv python --file pip-requirements.txt
$ conda activate my_venv
$ pip install -e .
```
In order to generate reference solutions for the example problems, the PySB package (and StochKit2) has to be installed additionally. PySB is installed in an anaconda environment via
```shell
$ conda install -c alubbock pysb
```

## Reference solutions for the example problems

