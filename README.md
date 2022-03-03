[![Licence: GPL v3](https://img.shields.io/badge/Licence-GPLv3-blue.svg)](LICENCE)

# xy-type-models

xy-type-models is an open-source, hybrid Fortran-Python application that implements the event-chain and Metropolis 
Monte Carlo algorithms for the simulation of two-dimensional XY-type models in statistical physics. 

Event-chain and Metropolis simulation is available for the XY and harmonic XY (HXY) spin models.  Metropolis simulation 
is available for the Maggs lattice-field electrolyte model in the grand canonical ensemble (for particles). Each model 
is defined on a two-dimensional square lattice. 

We provide uniform- and Gaussian-noise versions of the Metropolis algorithm for the XY and HXY models. The latter are 
labelled appropriately.

We provide multivalued and elementary versions of the lattice-field electrolyte. In the former, the charge value at 
each lattice site can be any integer multiple of the elementary charge, q; in the latter, the charge values are zero or 
±q. In order to easily compare with the XY and HXY models, we have set q = 2\pi.

For an introduction to the three-dimensional Maggs lattice-field electrolyte model in the canonical ensemble (for 
particles), see [\[Maggs2002\]](https://doi.org/10.1103/PhysRevLett.88.196402). For an introduction to both the 
two-dimensional lattice-field electrolyte model in the grand canonical ensemble (for particles) and its equivalence 
with the two-dimensional Villain model of magnetism ([\[Villain1975\]](
https://doi.org/10.1051/jphys:01975003606058100)), see [\[Faulkner2015\]](https://doi.org/10.1103/PhysRevB.91.155412) 
(Section II and Appendices B and C are particularly useful; the paper itself demonstrates the power of the model in 
characterising the nonergodic phase of the Berezinskii-Kosterlitz-Thouless transition). For an analysis of the 
similarities between the HXY model and two-dimensional lattice-field electrolyte, see [\[Faulkner2017\]](
https://doi.org/10.1088/1361-648X/aa523f).


## Installation

To install xy-type-models, clone this repository, navigate to the top xy-type-models directory and run `make`. The 
`make` command creates the eight Fortran executables, and stores them in a new directory called `executables`. The 
Fortran executables simulate the Markov processes. Their corresponding source code is contained in the [`src`](src) 
directory.

The code that analyses the resultant samples (i.e., that contained in the [`sample_analysis`](sample_analysis) 
directory) was written in Python and depends on [`numpy`](https://numpy.org). Some of it also depends on 
[`matplotlib`](https://matplotlib.org), [`rpy2`](https://rpy2.github.io) or [`scipy`](https://www.scipy.org). To manage 
external Python packages, we use [conda](https://docs.conda.io/projects/conda/en/latest/) environments via the 
[miniconda distribution](https://docs.conda.io/en/latest/miniconda.html). However, we found [`rpy2`](
https://rpy2.github.io) to be buggy when installed via conda. Instead, we `pip install rpy2` from within the project's 
conda environment (after having `conda install`ed [`numpy`](https://numpy.org), [`matplotlib`](https://matplotlib.org) 
and [`scipy`](https://www.scipy.org)).

[`markov_chain_diagnostics.py`](sample_analysis/markov_chain_diagnostics.py) depends on the R packages 
[`LaplacesDemon`](https://cran.r-project.org/web/packages/LaplacesDemon/) and [`mcmcse`](
https://cran.r-project.org/web/packages/mcmcse/). To install these R packages: download the binaries [here](
https://cran.r-project.org/web/packages/LaplacesDemon/) and [here](https://cran.r-project.org/web/packages/mcmcse/) 
and then run `R CMD INSTALL <binary location>` in your terminal.

The Fortran code was written in Fortran 90. The Python code is likely to support any Python version >= 3.6 (though we 
need to check this). We tested the Fortran / Python code with GNU Fortran (Homebrew GCC 10.2.0_4) 10.2.0 / CPython.


## Implementation 

The user interface of xy-type-models consists of the [`run.py`](run.py) script and a configuration file. The [`run.py`](
run.py) script expects the path to the configuration file as the first positional argument. Configuration files should 
be located in the [`config_files`](config_files) directory. 

To run xy-type-models, open your terminal, navigate to the top directory and enter `python run.py 
<configuration file>`. The generated sample data will appear in the location corresponding to the value of 
`output_directory` in the configuration file. Sample analysis can then be performed via the Python scripts contained in 
the [`sample_analysis`](sample_analysis) directory. Most sample-analysis scripts expect the path to the configuration 
file as the first positional argument, but there are some exceptions.


## Configuration files

Configuration files are located in the [`config_files`](config_files) directory. Their contents must follow the 
strict algorithm-dependent orders given in the subsections below. All floats must be given in non-exponential form to 
14 significant figures and with the suffix `d0`. All strings must be enclosed between two apostrophes, e.g., 
`'string'`. Typically, the value of `output_directory` starts with `'output/`, i.e., the entire value is of the form 
`'output/rest_of_string'`, but this is not a requirement ([`output`](output) is the xy-type-models output 
directory, but we often replace this with `../bc4-output` – the relative location of the output directory on our 
cluster machine). For the value of `output_directory`, i) refrain from giving long strings, as this can lead to Fortran 
runtime errors, and ii) ensure three tabbed spaces between the value of `output_directory` and the word 
`output_directory` itself, as this avoids Fortran errors when performing multiple runs of the same simulation (the 
number of runs is set by the value of `no_of_jobs`). The minimum / maximum value of `no_of_temperature_increments` is 
0 / 99. The values of `no_of_jobs` and `max_no_of_cpus` must be positive integers. For the value of `max_no_of_cpus`, 
we recommend giving half the number of CPUs available on your personal machine, e.g., for a four-core machine with two 
threads per core, we set `max_no_of_cpus = 4`; if `no_of_jobs = 8`, xy-type-models will perform two sets of four 
parallel runs of the same simulation, and similarly for certain sample-analysis processes.

### hxy-ecmc configuration file (an example)

```
'hxy-ecmc'                                  algorithm_name
'output/convergence_tests/hxy/ecmc'         output_directory
8                                           integer_lattice_length
100                                         no_of_equilibration_sweeps
40000                                       no_of_observations
1.3d0                                       initial_temperature
1.3d0                                       final_temperature
0                                           no_of_temperature_increments
0                                           vacuum_permittivity_sum_cutoff
.false.                                     randomise_initial_field_configuration
.false.                                     use_external_global_moves
.false.                                     calculate_potential_minimising_twists
1                                           no_of_jobs
4                                           max_no_of_cpus
```

### hxy-metropolis configuration file (an example)

```
'hxy-metropolis'                                    algorithm_name
'output/convergence_tests/hxy/metropolis'           output_directory
8                                                   integer_lattice_length
10000                                               no_of_equilibration_sweeps
100000                                              no_of_observations
1.3d0                                               initial_temperature
1.3d0                                               final_temperature
0                                                   no_of_temperature_increments
1.0d0                                               width_of_proposal_interval (initial)
0.44d0                                              target_acceptance_rate_of_field_rotations
0                                                   vacuum_permittivity_sum_cutoff
.false.                                             randomise_initial_field_configuration
.false.                                             use_external_global_moves
.false.                                             calculate_potential_minimising_twists
1                                                   no_of_jobs
4                                                   max_no_of_cpus
```

### hxy-gaussian-noise-metropolis configuration file (an example)

```
'hxy-gaussian-noise-metropolis'                                     algorithm_name
'output/convergence_tests/hxy/gaussian_noise_metropolis'            output_directory
8                                                                   integer_lattice_length
10000                                                               no_of_equilibration_sweeps
100000                                                              no_of_observations
1.3d0                                                               initial_temperature
1.3d0                                                               final_temperature
0                                                                   no_of_temperature_increments
1.0d0                                                               width_of_proposal_interval (initial)
0.44d0                                                              target_acceptance_rate_of_field_rotations
0                                                                   vacuum_permittivity_sum_cutoff
.false.                                                             randomise_initial_field_configuration
.false.                                                             use_external_global_moves
.false.                                                             calculate_potential_minimising_twists
1                                                                   no_of_jobs
4                                                                   max_no_of_cpus
```

### xy-ecmc configuration file (an example)

```
'xy-ecmc'                                   algorithm_name
'output/convergence_tests/xy/ecmc'          output_directory
8                                           integer_lattice_length
100                                         no_of_equilibration_sweeps
40000                                       no_of_observations
0.8d0                                       initial_temperature
0.8d0                                       final_temperature
0                                           no_of_temperature_increments
.false.                                     randomise_initial_field_configuration
.false.                                     use_external_global_moves
1                                           no_of_jobs
4                                           max_no_of_cpus
```

### xy-metropolis configuration file (an example)

```
'xy-metropolis'                                     algorithm_name
'output/convergence_tests/xy/metropolis'            output_directory
8                                                   integer_lattice_length
10000                                               no_of_equilibration_sweeps
100000                                              no_of_observations
0.8d0                                               initial_temperature
0.8d0                                               final_temperature
0                                                   no_of_temperature_increments
1.0d0                                               width_of_proposal_interval (initial)
0.44d0                                              target_acceptance_rate_of_field_rotations
.false.                                             randomise_initial_field_configuration
.false.                                             use_external_global_moves
1                                                   no_of_jobs
4                                                   max_no_of_cpus
```

### xy-gaussian-noise-metropolis configuration file (an example)

```
'xy-gaussian-noise-metropolis'                                  algorithm_name
'output/convergence_tests/xy/gaussian_noise_metropolis'         output_directory
8                                                               integer_lattice_length
10000                                                           no_of_equilibration_sweeps
100000                                                          no_of_observations
0.8d0                                                           initial_temperature
0.8d0                                                           final_temperature
0                                                               no_of_temperature_increments
1.0d0                                                           width_of_proposal_interval (initial)
0.44d0                                                          target_acceptance_rate_of_field_rotations
.false.                                                         randomise_initial_field_configuration
.false.                                                         use_external_global_moves
1                                                               no_of_jobs
4                                                               max_no_of_cpus
```

### elementary-electrolyte configuration file (an example)

```
'elementary-electrolyte'                                    algorithm_name
'output/convergence_tests/electrolyte/elementary'           output_directory
8                                                           integer_lattice_length
10000                	                                    no_of_equilibration_sweeps
100000	                                                    no_of_observations
1.5d0                 	                                    initial_temperature
1.5d0                                                       final_temperature
0                                                           no_of_temperature_increments
1.0d0       	                                            width_of_proposal_interval (initial)
0.44d0                                                      target_acceptance_rate_of_field_rotations
0.66666666666666d0                                          charge_hop_proportion
.false.                                                     use_external_global_moves
1                                                           no_of_jobs
4                                                           max_no_of_cpus
```

### multivalued-electrolyte configuration file (an example)

```
'multivalued-electrolyte'                                   algorithm_name
'output/convergence_tests/electrolyte/multivalued'          output_directory
8                                                           integer_lattice_length
10000                	                                    no_of_equilibration_sweeps
100000	                                                    no_of_observations
1.5d0                 	                                    initial_temperature
1.5d0                                                       final_temperature
0                                                           no_of_temperature_increments
1.0d0       	                                            width_of_proposal_interval (initial)
0.44d0                                                      target_acceptance_rate_of_field_rotations
0.66666666666666d0                                          charge_hop_proportion
.false.                                                     use_external_global_moves
1                                                           no_of_jobs
4                                                           max_no_of_cpus
```


## Convergence tests

In order to test the convergence of our code, we run (for example) 
`python run.py config_files/convergence_tests/xy/ecmc.txt` followed by 
`python sample_analysis/test_convergence.py config_files/convergence_tests/xy/ecmc.txt` (and similarly for the other 
seven algorithms). 

This computes the effective sample size of the sample generated by the Markov process and produces a plot of the 
cumulative distribution functions of both a reference sample and the generated sample. If the two cumulative 
distribution functions agree, we can be fairly sure that the code is working. Unit tests would be preferable, but we 
didn't have time to create any.

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
Metropolis/event-chain simulation. The slow functions would then be rewritten using the Fortran (or perhaps equivalent 
C) code of xy-type-models. We would then benchmark super-aLby against xy-type-models; if the total CPU times prove to 
be similar, this would provide evidence for Fortran/C code contained within a mediator-based, object-oriented Python 
structure being the optimal approach to coding in statistical physics (since Python is relatively easy to read and 
write and contains a lot of functionality).


## Makefiles

In the top directory, the `make` command runs the [`makefile`](makefile) contained there. By running eight different 
makefiles, this creates all eight executables (`xy_ecmc_algorithm.exe`, `xy_metropolis_algorithm.exe`, 
`xy_gaussian_noise_metropolis_algorithm.exe`, `hxy_ecmc_algorithm.exe`, `hxy_metropolis_algorithm.exe`, 
`hxy_gaussian_noise_metropolis_algorithm.exe`, `elementary_electrolyte_algorithm.exe` and 
`multivalued_electrolyte_algorithm.exe`) and stores them in a new directory called `executables`.

Each makefile is located in the youngest child directory corresponding to the relevant algorithm, e.g., the 
[`makefile`](src/xy_models/xy/ecmc/makefile) for `xy_ecmc_algorithm.exe` is contained in [the xy-ecmc directory](
src/xy_models/xy/ecmc). To create a single Fortran executable, open your terminal, navigate to the top xy-type-models 
directory and enter `make xy-ecmc`, `make xy-metropolis`, `make xy-gaussian-noise-metropolis`, `make hxy-ecmc`, 
`make hxy-metropolis`, `make hxy-gaussian-noise-metropolis`, `make elementary-electrolyte` or 
`make multivalued-electrolyte`. This will create the corresponding executable and store it in the `executables` 
directory.
