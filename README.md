# DLR approximation for the kinetic CME

## Requirements
- git (2.32.1 or later)
- CMake (3.22.1 or later)

macOS: Note that XCode compilers do not support OpenMP. To overcome this, a manual installation (e.g. of `gcc`) is required and the `CXX`, `CC` and `FC` environment variables have to be set accordingly.

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