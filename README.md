[![Licence: GPL v3](https://img.shields.io/badge/Licence-GPLv3-blue.svg)](LICENCE)

# xy-type-models

xy-type-models is an open-source Fortran/Python application that implements the event-chain and Metropolis-Hastings 
Monte Carlo algorithms for the simulation of two-dimensional XY-type models in statistical physics. 
Event-chain and Metropolis-Hastings simulation is available for the XY and harmonic XY (HXY) spin models. 
Metropolis-Hastings simulation is available for the Maggs lattice-field electrolyte model in the grand canonical 
ensemble (for particles). Each model is defined on a two-dimensional square lattice. 

We provide multivalued and elementary versions of the lattice-field electrolyte. In the former, the charge value at 
each lattice site can be any integer multiple of the elementary charge, q; in the latter, the charge values are zero or 
±q. In order to easily compare with the XY and HXY models, we have set q = 2\pi.

For an introduction to the three-dimensional Maggs lattice-field electrolyte model in the canonical ensemble (for 
particles), see [\[Maggs2002\]](https://doi.org/10.1103/PhysRevLett.88.196402). For an introduction to both the 
two-dimensional lattice-field electrolyte model in the grand canonical ensemble (for particles) and its equivalence 
with the two-dimensional Villain model of magnetism ([\[Villain1975\]](
https://doi.org/10.1051/jphys:01975003606058100)), see [\[Faulkner2015\]](https://doi.org/10.1103/PhysRevB.91.155412). 
For an analysis of the similarities between the HXY model and two-dimensional lattice-field electrolyte, see 
[\[Faulkner2017\]](https://doi.org/10.1088/1361-648X/aa523f).

## Installation

To install xy-type-models, clone this repository.

The code that simulates the Markov processes was written in Fortran 90. The code that analyses the resultant samples 
(i.e., that contained in the [`output`](output) directory) was written in Python and is likely to support any Python 
version >= 3.6 (though we need to check this). The sample-analysis code has been tested with CPython.

The sample-analysis code depends on [`numpy`](https://numpy.org). Some of it also depends on [`matplotlib`](
https://matplotlib.org) and [`rpy2`](https://rpy2.github.io). To manage external Python packages, we use [conda](
https://docs.conda.io/projects/conda/en/latest/) environments via the [miniconda distribution](
https://docs.conda.io/en/latest/miniconda.html). However, we found [`rpy2`](https://rpy2.github.io) to be buggy when 
installed via conda. Instead, we `pip install rpy2` from within the project's conda environment (after having `conda 
install`ed [`numpy`](https://numpy.org) and [`matplotlib`](https://matplotlib.org)).

[`markov_chain_diagnostics.py`](output/markov_chain_diagnostics.py) depends on the R packages [`LaplacesDemon`](
https://cran.r-project.org/web/packages/LaplacesDemon/) and [`mcmcse`](
https://cran.r-project.org/web/packages/mcmcse/). To install these R packages: download the binaries [here](
https://cran.r-project.org/web/packages/LaplacesDemon/) and [here](https://cran.r-project.org/web/packages/mcmcse/) 
and then run `R CMD INSTALL <binary location>` in your terminal.

## Implementation

To create the Fortran executables, open your terminal, navigate to your xy-type-models directory and enter `make`. This 
command creates all six executables (`xy_ecmc_algorithm.exe`, `xy_metropolis_algorithm.exe`, `hxy_ecmc_algorithm.exe`, 
`hxy_metropolis_algorithm.exe`, `elementary_electrolyte_algorithm.exe` and `multivalued_electrolyte_algorithm.exe`) by 
running six different makefiles.

Each makefile is located in the youngest child directory corresponding to the relevant 
algorithm, e.g., the [`makefile`](src/xy_models/xy/ecmc/makefile) for `xy_ecmc_algorithm.exe` is contained in [the 
xy-ecmc directory](src/xy_models/xy/ecmc). To create a single Fortran executable, open your terminal, navigate to your 
xy-type-models directory and enter `make xy-ecmc`, `make xy-metropolis`, `make hxy-ecmc`, `make hxy-metropolis`, 
`make elementary-electrolyte` or `make multivalued-electrolyte`. This will make the corresponding executable.

The user interface of xy-type-models consists of an executable and a configuration file. Each executable expects the 
path to the configuration file as the first positional argument. Configuration files should be located in the [
`config_files`](config_files) directory. To run an executable, stay in your xy-type-models directory and enter, e.g., 
`./xy_ecmc_algorithm <configuration file>`. The generated sample data will appear in the [`output`](output) directory 
(at a location given in the configuration file). Sample analysis can then be performed via scripts contained in the [
`output`](output) directory.

## Convergence tests

In order to test the convergence of our code, we run `./xy_ecmc_algorithm config_files/convergence_tests/xy/ecmc.txt` 
followed by `python output/convergence_tests/test_convergence.py config_files/convergence_tests/xy/ecmc.txt` (and 
similarly for the other five executables). This computes the effective sample size of the sample generated by the Markov 
process and produces a plot of the cumulative distribution functions of both a reference sample and the generated 
sample. If the two cumulative distribution functions agree, we can be fairly sure that the code is working.

Unit tests would be preferable, but we didn't have time to create any.

## Code structure

The structure of the Fortran 90 code is very much non-beautiful! This is because Fortran 90 does not lend itself well 
to object-oriented programming. However, Fortran is a very fast language, so xy-type-models is useful for numerical 
experiments that require large system and sample sizes, e.g., the analysis of nonergodic correlations in the 
Berezinskii-Kosterlitz-Thouless phase.

While xy-type-models does not use classes, it is still relatively modular. The main program scripts are [
`ecmc_algorithm.f90`](src/ecmc_algorithm.f90) and [`metropolis_algorithm.f90`](src/metropolis_algorithm.f90). Each main 
program calls many subroutines, each of which is contained in its eponymous file, e.g., [
`pre_simulation_processes.f90`](src/pre_simulation_processes.f90) (for subroutines common to all models and simulation 
methods) and [`xy_ecmc_read_config_file.f90`](src/xy_models/xy/ecmc/xy_ecmc_read_config_file.f90) (for subroutines 
common to, e.g., the xy-ecmc algorithm). There are a couple of exceptions: each subroutine that is only used in one 
other program or subroutine is located beneath that program or subroutine.

In the future, we would like to integrate these XY-type models into [super-aLby](
https://github.com/michaelfaulkner/super-aLby), our object-oriented Python application for super-relativistic Monte 
Carlo simulation. super-aLby uses the [mediator-design pattern](https://en.wikipedia.org/wiki/Mediator_pattern) and is 
therefore highly modular and beautiful. This would require additional super-aLby functionality for both the models and 
Metropolis-Hastings/event-chain simulation. The slow functions would then be rewritten using the Fortran (or perhaps 
equivalent C) code of xy-type-models. We would then benchmark super-aLby against xy-type-models; if the total CPU times 
prove to be similar, this would provide evidence for Fortran/C code contained within a mediator-based, object-oriented 
Python structure being the optimal approach to coding in statistical physics.
