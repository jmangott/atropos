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

$`P(t,\,x)\,\mathrm{d}t`$ is the probability of finding a population number of $`x = (x_1, \dots, x_N)`$ molecules of species $`S_1, \dots, S_N`$ in the time interval $`[t, t + \mathrm{d}t]`$.
The CME describes the time evolution of this probability distribution $`P(t,\,x)`$ in a chemical reaction network with $`N`$ different species $`S_1, \dots, S_N`$, which can react via $`M`$ reaction channels $`R_1, \dots, R_M`$. For a given reaction $`\mu`$, the stoichiometric vector $`\nu_\mu`$ denotes the population change by that reaction and the propensity functions $`a_\mu(x)`$ and $`a_\mu(x)`$ can be interpreted as transition probabilities $`T(x+\nu_\mu|x)`$ and $`T(x|x-\nu_\mu)`$.

In our DLR approach, the reaction network has to be separated into two parts, such that $`x=(x_{(1)},\,x_{(2)})`$. The probability distribution is then approximated by
```math
P \approx \sum_{i,j=1}^r X_i^1(t,\,x_{(1)})\,S_{ij}(t)\,X_i^2(t,\,x_{(2)})
```
with rank $`r`$, low-rank factors $`X_i^1(t,\,x_{(1)})`$ and $`X_i^2(t,\,x_{(2)})`$ and coupling coefficients $`S_{ij}(t)`$. The rank is usually a small number.

`kinetic-cme` makes use of the low-rank framework `Ensign`.[^fn2]

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
ctest --test-dir build
```

### Intel MKL
If you prefer to use Intel MKL as the BLAS and LAPACK backend instead of OpenBLAS set 
```shell
export MKLROOT=/path/to/intel/mkl
cmake -B build -DCMAKE_BUILD_TYPE=Release -DMKL_ENABLED=ON ..
cmake --build build
```
and make sure to add the MKL libraries to your `LD_LIBRARY_PATH`, i.e.
```shell
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/intel/mkl/lib/intel64_lin/
```
before running the executable.

### OpenMP
OpenMP can be activated via
```shell
cmake -B build -DCMAKE_BUILD_TYPE=Release -DOPENMP=ON ..
```
Make sure that the `OMP_NUM_THREADS` environment variable is in accordance with your hardware specification and run the unit tests folder via 
```shell
ctest --test-dir build
```
to ensure that OpenMP and `kinetic-cme` work correctly.

MacOS: Note that XCode compilers do not support OpenMP. For using OpenMP on macOS, a manual installation (e.g. of `gcc11`) is required and the `CXX`, `CC` and `FC` environment variables have to be set accordingly.

### Python environment
To use the Python notebooks and programs included in `scripts`, a Python environment with the packages specified in `pip-requirements.txt` needs to be configured and enabled.
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
All scripts and notebooks have to be executed from the project root. When using a IDE, make sure to adjust the settings accordingly. In Microsoft Visual Studio Code one has to set "Notebook File Root" to `{workspaceFolder}` to run the notebooks.

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

<!-- TODO: ### Binning -->

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

Note that when writing the model file and implementing the propensity functions one has to keep in mind the following partition convention of `kinetic-cme`, i.e. species with indices from `0` to `kM1`-1 belong to partition 1, species with indices from `kM1` to `kM1`+`kM2`-1 to partition 2.

It is recommended to use the exisiting model files for the example problems as a template for your own model.

`kinetic-cme` has to be recompiled if the model file is changed.

### Preparing input data
A template Python script called `set_input_template.py` is provided in `scripts/input` in order to facilitate the generation of the parameters (`parameters.hpp`) and the inital condition (`input.nc`). Code marked with `TODO` has to be modified according to the specific needs.
There are currently two different methods for generating the initial condition implemented:
1. `SetInputKD`: $`P(t=0,\,x) = \delta_{x,x_0}`$, where $`x_0`$ has to be specified. In this implementation the low-rank structure is exploited.
2. `SetInputGeneral`: $`P(t=0,\,x) = P_0(x)`$, where $`P_0(x)`$ has to be specified. In this implementation no low-rank structure is exploited, therefore it should be used only for small problems.

## Output
`kinetic-cme` automatically creates a folder in `output/` named `kFilename`, which is set in the `parameter.hpp` file.
The low-rank factors and coupling coefficients as well as the chosen model parameters are stored in this folder as `output_t<ts>.nc` (`<ts>` denotes the time step) in intervals according to the parameter `kSnapshot`.

The structure of the `output_t<ts>.nc` files is as follows:

|Variable name|Dimension|Description|
| --- | --- | --- |
|`X`|`(n_basisfunctions, dx1)`|Low-rank factors $`X_i^1`$|
|`S`|`(r, r)`|Coupling coefficients $`S_{ij}`$|
|`V`|`(n_basisfunctions, dx2)`|Low-rank factors $`X_j^2`$|
|`names`|`(n_basisfunctions, dx2)`|Species names|
|`n1`|`(n_basisfunctions, dx2)`|Population numbers in partition 1|
|`n2`|`(n_basisfunctions, dx2)`|Population numbers in partition 2|
|`binsize`|`(d)`|Concatenation of `binsize1` and `binsize2`|
|`liml`|`(d)`|Concatenation of `liml1` and `liml2`|
|`t`|`0`|Time point|
|`dt`|`0`|Time step size `kTau`|

## Example problems
Input generation scripts for the example problems (toggle switch, lambda phage and BAX pore aggregation) are provided in `scripts/input/example_setups`, model files can be found in `include/models`.

Interactive Python notebooks for comparison of the DLR results with reference solutions for the example problems can be found in `scripts/output`.

Before executing the notebooks, one has to generate the output files with `kinetic-cme` and the solvers for the reference solutions.

For toggle switch:
```shell
python3 scripts/reference_solutions/ode_ts.py
python3 scripts/input/example_setups/set_ts.py --tstar 500 --tau 0.01 --snapshot 100 --fname ts
cmake --build build
./bin/kinetic-cme
```
For lambda phage:
```shell
python3 scripts/reference_solutions/pysb_stochkit.py
python3 scripts/input/example_setups/set_lp.py --tstar 10 --tau 0.01 --so --substeps 10 --snapshot 100 --fname lp
cmake --build build
./bin/kinetic-cme
```

## References
[^fn1]: Lubich. C., Oseledets, I.: "A projector-splitting integrator for dynamical low-rank approximation", BIT Numerical Mathematics **54** (2014)

[^fn2]: Cassini, F., Einkemmer, L.: "Efficient 6D Vlasov simulation using the dynamical low-rank framework Ensign", Computer Physics Communications **280** (2022)