# DLR approximation for the kinetic CME

- [DLR approximation for the kinetic CME](#dlr-approximation-for-the-kinetic-cme)
  - [Objective](#objective)
  - [Requirements](#requirements)
  - [Installation](#installation)
    - [Intel MKL](#intel-mkl)
    - [OpenMP](#openmp)
    - [Python environment](#python-environment)
  - [Run the program](#run-the-program)
  - [Input](#input)
    - [Integrators](#integrators)
    - [Preparing input data](#preparing-input-data)
      - [Writing a model file with the `ReactionSystem` class](#writing-a-model-file-with-the-reactionsystem-class)
  - [Output](#output)
  - [Example problems](#example-problems)
  - [References](#references)

## Objective
`hierarchical-cme` solves the chemical master equation (CME),
```math
\partial_t{P}(t,x) = \sum_{\mu = 1}^{M}\left(a_\mu(x-\nu_\mu)P(t,x-\nu_\mu) - a_\mu(x)P(t,x)\right)
```
with the projector-splitting integrator for Tree Tensor networks.[^fn1]

$`P(t,x)\,\mathrm{d}t`$ is the probability of finding a population number of $`x = (x_1, \dots, x_N)`$ molecules of species $`S_1, \dots, S_N`$ in the time interval $`[t,\,t + \mathrm{d}t]`$.
The CME describes the time evolution of this probability distribution $`P(t,x)`$ in a chemical reaction network with $`N`$ different species $`S_1, \dots, S_N`$, which can react via $`M`$ reaction channels $`R_1, \dots, R_M`$. For a given reaction $`\mu`$, the stoichiometric vector $`\nu_\mu`$ denotes the population change by that reaction and the propensity functions $`a_\mu(x)`$ and $`a_\mu(x-\nu_\mu)`$ are proportional to the transition probabilities $`T(x+\nu_\mu|x)`$ and $`T(x|x-\nu_\mu)`$.

<!-- In our approach, the reaction network has to be separated recursively into two parts, such that $`x=(x_{(1)},x_{(2)})`$. The probability distribution is then approximated by
```math
P(t,x) \approx \sum_{i,j=1}^r X_i^1(t,x_{(1)})\,S_{ij}(t)\,X_i^2(t,x_{(2)})
```
with rank $`r`$, low-rank factors $`X_i^1(t,x_{(1)})`$ and $`X_i^2(t,x_{(2)})`$ and coupling coefficients $`S_{ij}(t)`$. The rank is usually a small number. -->

`hierarchical-cme` makes use of the low-rank framework `Ensign`.[^fn2]

## Requirements
- CMake (3.22.1 or later)
- C++20 compatible C++ compiler
- Eigen 3.4 (if the implicit Euler integrator is used)
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
git clone git@github.com:jmangott/CME-integrator.git
```

<!-- ```shell
# git clone https://git.uibk.ac.at/c7021158/kinetic-cme.git
``` -->

and build the program by executing
```shell
cd CME-integrator
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

<!-- ```shell
cd kinetic-cme
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
``` -->

The generated executable `hierarchical-cme` can be found in `bin`.

To enable compiler options for debugging, use `-DCMAKE_BUILD_TYPE=Debug` instead.
Unit tests for C++ files are provided in the `tests` folder. They can be run with 
```shell
ctest --test-dir build
```

### Intel MKL
If you prefer to use Intel MKL as the BLAS and LAPACK backend instead of OpenBLAS set 
```shell
export MKLROOT=/path/to/intel/mkl
cmake -B build -DCMAKE_BUILD_TYPE=Release -DMKL_ENABLED=ON
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
cmake -B build -DCMAKE_BUILD_TYPE=Release -DOPENMP=ON
```
Make sure that the `OMP_NUM_THREADS` environment variable is in accordance with your hardware specification and run the unit tests via 
```shell
ctest --test-dir build
```
to ensure that OpenMP and `hierarchical-cme` work correctly.

MacOS: Note that XCode compilers do not support OpenMP. For using OpenMP on macOS, a manual installation (e.g. of `gcc11`) is required and the `CXX`, `CC` and `FC` environment variables have to be set accordingly.

### Python environment
To use the Python notebooks and programs included in `scripts`, a Python environment with external packages specified in `pip-requirements.txt` needs to be configured and enabled.
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
All scripts and notebooks have to be executed from the project root. When using a IDE, make sure to adjust the settings accordingly. In Microsoft Visual Studio Code one has to set the "Notebook File Root" to `{workspaceFolder}` to run the notebooks.
Unit tests for Python files are located in the `scripts/tests` folder. They can be run in the Python environment via
```shell
pytest scripts/tests
```

## Run the program
`hierarchical-cme` has to be run with
```
  ./bin/hierarchical-cme [OPTION...]
```
and expects the following command line arguments:
- `-o`, `--output`: Name of the output folder, stored in `output/`
- `-s`, `--snapshot`: Number of steps between two snapshots
- `-t`, `--tau`: Time step size
- `-f`, `--tfinal`: Final integration time
- `-h`, `--help`: Print usage

## Input
Input netCDF files have to be stored as `input/input.nc` and can be generated with the input scripts provided in `scripts/input`.

**Caution:** As `Ensign` stores arrays in column-major (Fortran) order, it is assumed that input arrays also follow this convention.
<!-- TODO: Give more detais -->

### Integrators
The user can choose between a first-order explicit and implicit Euler integrator for the K step, where the low-rank factors depending on the population numbers are updated.
<!-- TODO: ### Binning -->

### Preparing input data
The `set_bax.py` script located in `scripts/input` folder generates input data for the BAX pore assembly model. It gives an example on how the initial conditions have to be set up. The `input/input.nc` file is generated with the `set_bax.py` script via
```shell
python3 scripts/input/set_bax.py --partition "(0 1 2)(((3 4 6 7)(5 8))(9 10))" --rank 5 15 15
```
and a short documentation for this script is provided by
```shell
python3 scripts/input/set_bax.py --help
```
<!-- TODO: ### Describe examples in more detail -->

Note that `hierarchical-cme` assumes that the propensity function is factorizable for the species in different partitions. However, the input scripts rely on the `ReactionSystem` class (cf. `scripts/reaction_class.py`), which assumes that the propensity function is factorizable in *all* species. This is a valid assumption for most scenarios. For problems where species in a partition are not factorizable, the propensity function can be adjusted manually after initializing the `Tree` with the method `initialize`.

#### Writing a model file with the `ReactionSystem` class
The model file contains all reactions $`R_\mu`$ ($`\mu=1,\dots,M`$) of the specific problem and has to be imported in the input scripts.

<!-- TODO: More detailed description. -->


## Output
`hierarchical-cme` automatically creates a folder in `output/` with a name set by the `-o`/`--output` parameter.
The low-rank factors and coupling coefficients as well as the chosen model parameters are stored in this folder as `output_t<ts>.nc` (`<ts>` denotes the time step) in intervals according to the `-s`/`--snapshot` parameter.

<!-- TODO: Describe the structure of the .netCDF file -->

## Example problems
Input generation scripts for the example problems (toggle switch, lambda phage and BAX pore aggregation) are provided in `scripts/input/` and model files can be found in `scripts/models`.

Interactive Python notebooks for comparison of the DLR results with reference solutions for the example problems can be found in `scripts/output`.

Before executing the notebooks, one has to generate the output files with `hierarchical-cme` and the corresponding reference solutions. These can in turn be generated with the scripts located in `scripts/reference_solutions`.

## References
[^fn1]: Ceruti, G., Lubich, C., and Walach, H.: Time integration of Tree Tensor networks", SIAM J. Numer. Anal. **59** (2021)
<!-- Lubich, C., Oseledets, I.: "A projector-splitting integrator for dynamical low-rank approximation", BIT Numerical Mathematics **54** (2014) -->

[^fn2]: Cassini, F., and Einkemmer, L.: "Efficient 6D Vlasov simulation using the dynamical low-rank framework Ensign", Comp. Phys. Comm. **280** (2022)