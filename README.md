# DLR approximation for the kinetic CME

## Requirements
- git (2.32.1 or later)
- CMake (3.22.1 or later)

## Installation
The repository makes use of Ensign (cf. https://github.com/leinkemmer/Ensign) as a git submodule. The easiest way to obtain all the relevant sources is to execute
```shell
    git clone --recursive https://git.uibk.ac.at/c7021158/kinetic-cme
```

To build the program execute
```shell
    mkdir build
    cd build
    cmake ..
    make
```