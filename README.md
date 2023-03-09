# DLR approximation for the kinetic CME

## Requirements
- git (2.32.1 or later)
- CMake (3.22.1 or later)

For performing snapshots, installation of
- HDF5 (1.10.x)
- netCDF-4
is required. Check via `nc-config --has-hdf5`, whether HDF5 was used in the netCDF-4 build.

macOS: Note that XCode compilers do not support OpenMP. For using OpenMP on macOS, a manual installation (e.g. of `gcc11`) is required and the `CXX`, `CC` and `FC` environment variables have to be set accordingly.

## Installation
Clone the repository via
```shell
git clone https://git.uibk.ac.at/c7021158/kinetic-cme
```
and build the program by executing
```shell
mkdir build
cd build
cmake ..
make
```