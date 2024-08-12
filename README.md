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
    - [Preparing input data](#preparing-input-data)
  - [Output](#output)
  - [Example problems](#example-problems)
  - [References](#references)

## Objective
`hierarchical-cme` solves the chemical master equation (CME),
```math
\partial_t{P}(t,x) = \sum_{\mu = 1}^{M}\left(a_\mu(x-\nu_\mu)P(t,x-\nu_\mu) - a_\mu(x)P(t,x)\right)
```
according to the algorithm proposed in https://arxiv.org/abs/2407.11792, which is based on the projector-splitting integrator for Tree Tensor networks.[^fn1]

$`P(t,x)\,\mathrm{d}t`$ is the probability of finding a population number of $`x = (x_1, \dots, x_N)`$ molecules of species $`S_1, \dots, S_N`$ in the time interval $`[t,\,t + \mathrm{d}t]`$.
The CME describes the time evolution of this probability distribution $`P(t,x)`$ in a chemical reaction network with $`N`$ different species $`S_1, \dots, S_N`$, which can react via $`M`$ reaction channels $`R_1, \dots, R_M`$. For a given reaction $`\mu`$, the stoichiometric vector $`\nu_\mu`$ denotes the population change by that reaction and the propensity functions $`a_\mu(x)`$ and $`a_\mu(x-\nu_\mu)`$ are proportional to the transition probabilities $`T(x+\nu_\mu|x)`$ and $`T(x|x-\nu_\mu)`$.

`hierarchical-cme` makes use of the low-rank framework `Ensign`.[^fn2]

## Requirements
- CMake (3.22.1 or later)
- C++20 compatible C++ compiler
- Eigen 3.4 (if the implicit Euler or Crank-Nicolson integrators are used)
- Fortran compiler (if OpenBLAS is used)
- HDF5 (1.10.x)
- netCDF4
- Python (>3.8)

Check via `nc-config --has-hdf5`, whether HDF5 was used in the netCDF4 build.

Optionally:
- OpenMP
- Intel MKL

## Installation
Clone the repository via
```shell
git clone https://git.uibk.ac.at/c7021158/kinetic-cme.git
```

<!-- ```shell
# git clone git@github.com:jmangott/CME-integrator.git
``` -->

and build the program by executing
```shell
cd kinetic-cme
cmake -B <build> -DCMAKE_BUILD_TYPE=Release
cmake --build <build>
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
ctest --test-dir <build>
```

### Intel MKL
If you prefer to use Intel MKL as the BLAS and LAPACK backend instead of OpenBLAS set 
```shell
export MKLROOT=/path/to/intel/mkl
cmake -B <build> -DCMAKE_BUILD_TYPE=Release -DMKL_ENABLED=ON
cmake --build <build>
```
and make sure to add the MKL libraries to your `LD_LIBRARY_PATH`, i.e.
```shell
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/intel/mkl/lib/intel64_lin/
```
before running the executable.

### OpenMP
OpenMP can be activated via
```shell
cmake -B <build> -DCMAKE_BUILD_TYPE=Release -DOPENMP=ON
```
Make sure that the `OMP_NUM_THREADS` environment variable is in accordance with your hardware specification and run the unit tests via 
```shell
ctest --test-dir <build>
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
- `-i`, `--input`: Name of the input .nc file (default: `input/input.nc`)
- `-o`, `--output`: Name of the output folder, stored in `output/`
- `-s`, `--snapshot`: Number of steps between two snapshots
- `-t`, `--tau`: Time step size
- `-f`, `--tfinal`: Final integration time
- `-n`, `--substeps`: Number of integration substeps (default: `1`)
- `-m`, `--method`: Integration method (`e` (explicit Euler), `r` (explicit RK4), `i` 
                      (implicit Euler), `c` (Crank-Nicolson)) (default: `i`)
- `-h`, `--help`: Print usage

## Input
Input netCDF files have to be stored as `input/input.nc` (the directory can be changed using the `-i` flag) and can be generated with the input scripts provided in `scripts/input_generation`.

**Caution:** As `Ensign` stores arrays in column-major (Fortran) order, it is assumed that input arrays also follow this convention.
<!-- TODO: Give more detais -->

<!-- TODO: ### Binning -->

### Preparing input data
Let us consider the input script `set_lambda_phage.py` located in the `scripts/input_generation` folder, which generates input data for the lambda phage model. It gives an example on how the initial conditions have to be set up. The `input/input.nc` file is generated via
```shell
python3 scripts/input_generation/set_lamba_phage.py --partition "(0 1)((2 3)(4))" --rank 5 4
```
and a short documentation for this script is provided by
```shell
python3 scripts/input_generation/set_lambda_phage.py --help
```
<!-- TODO: ### Describe examples in more detail -->

Note that `hierarchical-cme` assumes that the propensity function is factorizable for the species in different partitions. However, the input scripts rely on the `ReactionSystem` class (cf. `scripts/reaction_class.py`), which assumes that the propensity function is factorizable in *all* species. This is a valid assumption for most scenarios. For problems where species in a partition are not factorizable, the propensity function can be adjusted manually after initializing the `Tree` with the method `initialize`.

<!-- #### Writing a model file with the `ReactionSystem` class
The model file contains all reactions $`R_\mu`$ ($`\mu=1,\dots,M`$) of the specific problem and has to be imported in the input scripts. -->

<!-- TODO: More detailed description. -->


## Output
`hierarchical-cme` automatically creates a folder in `output/` with a name set by the `-o`/`--output` parameter.
The low-rank factors and coupling coefficients as well as the chosen model parameters are stored in this folder as `output_t<ts>.nc` (`<ts>` denotes the time step) in intervals according to the `-s`/`--snapshot` parameter.

<!-- TODO: Describe the structure of the .netCDF file -->

## Example problems
Input generation scripts for the example problems (lambda phage and reaction cascade) are provided in `scripts/input_generation` and the corresponding model files can be found in `scripts/models`.

All required output files and reference solutions for reproducing the plots in https://arxiv.org/abs/2407.11792 can be computed with the shell scripts provided in `scripts/shell`. Before generating the plots with the interactive Python notebooks provided in `scripts/output/notebooks`, a `plots` folder has to be created in the project root:
```shell
mkdir plots
```
Then, for the lambda phage example one has to run
```shell
sh scripts/shell/run_lambda_phage.sh
```
and for the cascade reaction example
```shell
sh scripts/shell/run_cascade.sh
```

## References
[^fn1]: Ceruti, G., Lubich, C., and Walach, H.: Time integration of Tree Tensor networks", SIAM J. Numer. Anal. **59** (2021)
<!-- Lubich, C., Oseledets, I.: "A projector-splitting integrator for dynamical low-rank approximation", BIT Numerical Mathematics **54** (2014) -->

[^fn2]: Cassini, F., and Einkemmer, L.: "Efficient 6D Vlasov simulation using the dynamical low-rank framework Ensign", Comp. Phys. Comm. **280** (2022)