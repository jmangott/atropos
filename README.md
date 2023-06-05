# DLR approximation for the kinetic CME

- [DLR approximation for the kinetic CME](#dlr-approximation-for-the-kinetic-cme)
  - [Objective](#objective)
  - [Requirements](#requirements)
  - [Installation](#installation)
    - [Intel MKL](#intel-mkl)
    - [OpenMP](#openmp)
    - [Python environment](#python-environment)
  - [Input](#input)
    - [Integrators](#integrators)
      - [First-order method](#first-order-method)
      - [Second-order method](#second-order-method)
    - [Writing a model file](#writing-a-model-file)
    - [Preparing input data](#preparing-input-data)
  - [Output](#output)
  - [Example problems](#example-problems)
  - [References](#references)

## Objective
`kinetic-cme` solves the chemical master equation (CME),
```math
\partial_t{P}(t,\,x) = \sum_{\mu = 1}^{M}\left(a_\mu(x-\nu_\mu)\,P(t,\,x-\nu_\mu) - a_\mu(x)\,P(t,\,x)\right)
```
with the projector-splitting based dynamical low-rank (DLR) approximation.[^fn1]

<!-- The CME describes the time evolution of the probability distribution `$P(t,\,x)$` in a chemical reaction network 

It makes use of the low-rank framework `Ensign`.[^fn2] -->


## Requirements
- CMake (3.22.1 or later)
- C++17 compatible C++ compiler
- Fortran compiler (if OpenBLAS is used)
- HDF5 (1.10.x)
- netCDF4
- Python (>3.8)

Check via `nc-config --has-hdf5`, whether HDF5 was used in the netCDF4 build.

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
Unit tests are provided in the `tests` folder and can be run with 
```shell
ctest
```

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
Make sure that the `OMP_NUM_THREADS` environment variable is in accordance with your hardware specification and run the unit tests in the `build` folder via 
```shell
ctest
```
to ensure that OpenMP and `kinetic-cme` work correctly.

MacOS: Note that XCode compilers do not support OpenMP. For using OpenMP on macOS, a manual installation (e.g. of `gcc11`) is required and the `CXX`, `CC` and `FC` environment variables have to be set accordingly.

### Python environment
To use the Python programs included in `scripts`, a Python environment with the packages specified in `pip-requirements.txt` needs to be configured and enabled.
For Python venv:
```shell
python -m venv path/to/my_venv
source path/to/my_venv/bin/activate
pip install -r pip-requirements.txt
pip install -e .
```
For anaconda:
```shell
conda create -n my_venv python --file pip-requirements.txt
conda activate my_venv
pip install -e .
```
In order to generate reference solutions for the example problems, the PySB package (and StochKit2) has to be installed additionally. PySB is installed in an anaconda environment via
```shell
conda install -c alubbock pysb
```

## Input
`kinetic-cme` reads the parameters from `parameters.hpp` (located in `include/`) and the initial condition for the low-rank factors and the coupling coefficients from `input.nc` (located in `input/`).

Following model parameters can be set in `parameters.hpp`:

|Parameter name|Type|Description|Condition|
| --- | --- | --- | --- |
|`kR`|`Index`|rank|>0|
|`kM1`|`Index`|# species in partition 1|>0|
|`kM2`|`Index`|# species in partition 2|>0|
|`kN1`|`array<Index>`|# grid points for species in partition 1|size `kM1`|
|`kN2`|`array<Index>`|# grid points for species in partition 2|size `kM2`|
|`kBinsize1`|`array<Index>`|binsize for species in partition 1|size `kM1`|
|`kBinsize2`|`array<Index>`|binsize for species in partition 2|size `kM2`|
|`kLiml1`|`array<Index>`|left population number limit for species in partition 1|size `kM1`|
|`kLiml2`|`array<Index>`|left population number limit for species in partition 2|size `kM2`|
|`kTstar`|`double`|final integration time||
|`kTau`|`double`|time step size|>0|
|`kSnapshot`|`Index`|# timesteps between snapshots||
|`kSecondOrder`|`bool`|flag for activiating the second-order method||
|`kNSubsteps`|`bool`|# explicit Euler substeps for second-order method||
|`kFilename`|`array<char>`|name of the output folder||

`Index` is a type name alias for `ptrdiff_t`.

Moreover, the model file (see [below](#writing-a-model-file)) for the reactions has to be included in `parameters.hpp`.

The structure of the `input.nc` file is as follows:

|Variable name|Dimension|Description|
| --- | --- | --- |
|`xx1`|`(n_basisfunctions, dx1)`|Low-rank factors $`X_i^1`$|
|`ss`|`(r, r)`|Coupling coefficients $`S_{ij}`$|
|`xx2`|`(n_basisfunctions, dx2)`|Low-rank factors $`X_j^2`$|

`n_basisfunctions` denotes the rank of the initial condition.
In general, `n_basisfunctions` = `r`, but the initial condition can be smaller than `r` (the fixed rank choosen for the entire simulation), for example when choosing Kronecker-$`\delta`$-like initial conditions (see [below](#preparing-input-data)).
In that case, `Ensign` will construct the remaining `r`-`n_basisfunctions` low-rank factors via Gram-Schmidt orthogonalization.

### Integrators
The user can choose between a first-order and a second-order method by setting `kSecondOrder` accordingly.

#### First-order method
When setting `kSecondOrder = false`, Lie-Trotter splitting and an explicit Euler method will be used for integration. The parameter `kNSubsteps` will not be effective in that case.

#### Second-order method
When setting `kSecondOrder = true`, Strang splitting and an explicit Euler method with substeps will be used for integration. `kNSubsteps` determines the number of substeps.

### Writing a model file
The model file contains all reactions $`R_\mu`$ ($`\mu=1,\dots,M`$) of the specific problem and has to be stored as a `.hpp` file in `include/models` in order to work with the input scripts. Reactions are passed to the instance `mysystem` of the  `mysys` class by constructing an instance of the `myreact` class with
```c++
myreact(vector<int> nu, vector<Index> depends_on, double (*prop_function)(vector<double>), mysys mysystem)
```
with stoichiometric vector `nu` (corresponding to $`\nu_\mu`$), propensity function `prop_function` (corresponding to $`a_\mu(x)`$). `depends_on` is a vector which references all species on which the propensity function $`a_\mu(x)`$ explicitly depends.

The `mysystem` class keeps track of all associated reactions, an instance can be constructed via
```c++
mysys(vector<string> species_names)
```
where `species_names` constains the names of all species in the system.
It is recommended to use the exisiting model files for the example problems as a template for your own model.

### Preparing input data
A template Python script called `set_input_template.py` is provided in `scripts/input` in order to facilitate the generation of the parameters (`parameters.hpp`) and the inital condition (`input.nc`). Code marked with `TODO` has to be modified according to the specific needs.
There are two different types of 

## Output
`kinetic-cme` automatically creates a folder in `output/` named `kFilename`, which is set in the `parameter.hpp` file.
The low-rank factors and coupling coefficients as well as the chosen model parameters are stored in theis folder as `output_t<ts>.nc` (`<ts>` denotes the time step) in intervals according to the parameter `kSnapshot`.

## Example problems
Input generation scripts for the example problems (toggle switch, lambda phage and BAX pore aggregation) are provided in `scripts/input/example_setups`, model files can be found in `include/models`.

## References
- [^fn1]: Lubich. C., Oseledets, I.: "A projector-splitting integrator for dynamical low-rank approximation", BIT Numerical Mathematics **54** (2014)
- [^fn2]: Cassini, F., Einkemmer, L.: "Efficient 6D Vlasov simulation using the dynamical low-rank framework Ensign", Computer Physics Communications **280** (2022)