# xy-type-models

Code Title:
ECMC algorithm for the 2dXY model (with Metropolis MCMC global spin twists to directly propose TSFs)

Author:
Michael Faulkner (with help from Manon Michel and Werner Krauth)

Purpose:
Draws samples from the PDF of the 2dXY model in order to estimate estimators of macro. quantities

MF - 09 Aug 2017

---------------------------------------------------------------------------

0. ENSURE that all of the files listed at the bottom of this README_ECMC_XY are in the correct directories;
1. Execute the make command (BASH COMMAND: make) to produce the executable ECMC_XY.exe in the source directory;
2. Copy ECMC_XY.exe to the directory in which you wish to run the Markov chain (the run directory)
   (BASH COMMAND (e.g.): cp ECMC_XY.exe ../L8run1);
3. Set the lattice length, no. of samples (measurements), etc. in initial.in in the run directory;
4. Execute ECMC_XY.exe (BASH COMMAND: ./ECMC_XY.exe) in the run directory: this will produce subdirectories for each temperature, then run the ECMC and draw samples from PDF at each temperatures: the samples are saved to files in the subdirectories;
5. Run the python script analyse_samples_ECMC_XY.py (BASH COMMAND: python analyse_samples_ECMC_XY.py) in the run directory to produce the estimates of the quantities of interest, which are printed to a file.

---------------------------------------------------------------------------

The initial input should be of the form:

side         	(length of one side of the lattice)
therm_sweeps	(number of thermalization sweeps)
measurements	(number of chains/measurements)
Tmin             (Lowest temperature)
Tmax             (Highest temperature)
Tstep            (Number of temperature steps between Tmin and Tmax)
start			(initialize spins: 0 = ordered, 1 = random)
twist      		(directly propose TSFs via global spin twists: 0 = no; 1 = yes) 
seed      		(seed to initialize the random number generator)

---------------------------------------------------------------------------

Files used in the program:

SOURCE FILES:
make                    Compiles the following files and produces the executable ECMC_XY.exe
ECMC_XY.f90             Main program
measure_ECMC_XY.f90     Draws samples when chain stops
setup_ECMC_XY           Sets PDF from initial.in; also sets all common variables in module variables
random_ECMC_XY.f90  	  RANMAR lagged Fibonacci random number generator
global_twist_XY.f90        Contains the global_twist subroutine, which proposes TSFs via global phase twists

ANALYSIS FILES (IN THE RUN DIRECTORY):
initial.in                  Sets lattice size, chain length, etc.
analyse_samples_ECMC_XY.py  Computes estimates of estimators from the saved Markov chain
MCMC_functions              External functions used in analyse_event_chain.py

---------------------------------------------------------------------------



Code Title:
Metropolis algorithm for the HXY model (with Metropolis MCMC global spin twists to directly propose TSFs)

Author:
Michael Faulkner

Purpose:
Draws samples from the PDF of the HXY model in order to estimate estimators of macro. quantities

MF - 09 Aug 2017

---------------------------------------------------------------------------

0. ENSURE that all of the files listed at the bottom of this README_Metrop_HXY are in the correct directories;
1. Execute the make command (BASH COMMAND: make) to produce the executable Metrop_HXY.exe in the source directory;
2. Copy Metrop_HXY.exe to the directory in which you wish to run the Markov chain (the run directory)
   (BASH COMMAND (e.g.): cp Metrop_HXY.exe ../L8run1);
3. Set the lattice length, no. of samples (measurements), etc. in initial.in in the run directory;
4. Execute Metrop_HXY.exe (BASH COMMAND: ./Metrop_HXY.exe) in the run directory: this will produce subdirectories for each temperature, then run the Markov chain and draw samples from PDF at each temperatures: the samples are saved to files in the subdirectories;
5. Run the python script analyse_samples_Metrop_HXY.py (BASH COMMAND: python analyse_samples_Metrop_HXY.py)
   to produce the estimates of the quantities of interest, which are printed to a file.

---------------------------------------------------------------------------

The initial input should be of the form:

side         	(length of one side of the lattice)
therm_sweeps	(number of thermalization sweeps)
measurements	(number of chains/measurements)
Tmin             (Lowest temperature)
Tmax             (Highest temperature)
Tstep            (Number of temperature steps between Tmin and Tmax)
start			(initialize spins: 0 = ordered, 1 = random)
twist      		(directly propose TSFs via global spin twists: 0 = no; 1 = yes) 
seed      		(seed to initialize the random number generator)

---------------------------------------------------------------------------

Files used in the program:

SOURCE FILES:
make                    Compiles the following files and produces the executable Metrop_HXY.exe
Metrop_HXY.f90             Main program
measure_Metrop_HXY.f90     Draws samples when chain stops
setup_Metrop_HXY.f90           Sets PDF from initial.in; also sets all common variables in module 
random_Metrop_HXY.f90  	  RANMAR lagged Fibonacci random number generator
global_twist_HXY.f90        Contains the global_twist subroutine, which proposes TSFs via global phase twists

ANALYSIS FILES (IN THE RUN DIRECTORY):
initial.in                  Sets lattice size, chain length, etc.
analyse_samples_Metrop_HXY.py Computes estimates of estimators from the saved Markov chain
MCMC_functions              External functions used in analyse_markov_chain.py

---------------------------------------------------------------------------

An optional change that also measures the TSFs by finding the number of global twists that minimizes the potential (USING AS DEFAULT ON 16 FEB 2018):

measure_Metrop_HXY.f90_with_twist_TSF_measure.f90
setup_Metrop_HXY_with_twist_TSF_measure.f90





Code Title:
Event-chain Monte Carlo (ECMC) algorithm for the HXY model (with Metropolis MCMC global spin twists to directly propose TSFs)

Author:
Michael Faulkner

Purpose:
Draws samples from the PDF of the HXY model in order to estimate estimators of macro. quantities

MF - 09 Aug 2017

---------------------------------------------------------------------------

0. ENSURE that all of the files listed at the bottom of this README_ECMC_HXY are in the correct directories;
1. Execute the make command (BASH COMMAND: make) to produce the executable ECMC_HXY.exe in the source directory;
2. Copy ECMC_HXY.exe to the directory in which you wish to run the Markov chain (the run directory)
   (BASH COMMAND (e.g.): cp ECMC_HXY.exe ../L8run1);
3. Set the lattice length, no. of samples (measurements), etc. in initial.in in the run directory;
4. Execute ECMC_HXY.exe (BASH COMMAND: ./ECMC_HXY.exe) in the run directory: this will produce subdirectories for each temperature, then run the Markov chain and draw samples from PDF at each temperatures: the samples are saved to files in the subdirectories;
5. Run the python script analyse_samples_ECMC_HXY.py (BASH COMMAND: python analyse_samples_ECMC_HXY.py)
   to produce the estimates of the quantities of interest, which are printed to a file.

---------------------------------------------------------------------------

The initial input should be of the form:

side         	(length of one side of the lattice)
therm_sweeps	(number of thermalization sweeps)
measurements	(number of chains/measurements)
Tmin             (Lowest temperature)
Tmax             (Highest temperature)
Tstep            (Number of temperature steps between Tmin and Tmax)
start			(initialize spins: 0 = ordered, 1 = random)
twist      		(directly propose TSFs via global spin twists: 0 = no; 1 = yes) 
seed      		(seed to initialize the random number generator)

---------------------------------------------------------------------------

Files used in the program:

SOURCE FILES:
make                    Compiles the following files and produces the executable ECMC_HXY.exe
ECMC_HXY.f90             Main program
measure_ECMC_HXY.f90     Draws samples when chain stops
setup_ECMC_HXY.f90           Sets PDF from initial.in; also sets all common variables in module 
random_ECMC_HXY.f90  	  RANMAR lagged Fibonacci random number generator
global_twist_HXY.f90        Contains the global_twist subroutine, which proposes TSFs via global phase twists

ANALYSIS FILES (IN THE RUN DIRECTORY):
initial.in                  Sets lattice size, chain length, etc.
analyse_samples_ECMC_HXY.py Computes estimates of estimators from the saved Markov chain
MCMC_functions              External functions used in analyse_markov_chain.py

---------------------------------------------------------------------------

An optional change that also measures the TSFs by finding the number of global twists that minimizes the potential (USING AS DEFAULT ON 16 FEB 2018):

measure_ECMC_HXY.f90_with_twist_TSF_measure.f90
setup_ECMC_HXY_with_twist_TSF_measure.f90








Code Title:
Metropolis algorithm for 2D generalized lattice electrostatics (with global TSF proposals)

Author:
Michael Faulkner

Purpose:
Draws samples from the PDF of the HXY model in order to estimate estimators of macro. quantities

MF - 09 Aug 2017

---------------------------------------------------------------------------

0. ENSURE that all of the files listed at the bottom of this README_Metrop_HXY are in the correct directories;
1. Execute the make command (BASH COMMAND: make) to produce the executable Metrop_HXY.exe in the source directory;
2. Copy Metrop_HXY.exe to the directory in which you wish to run the Markov chain (the run directory)
   (BASH COMMAND (e.g.): cp Metrop_HXY.exe ../L8run1);
3. Set the lattice length, no. of samples (measurements), etc. in initial.in in the run directory;
4. Execute Metrop_HXY.exe (BASH COMMAND: ./Metrop_HXY.exe) in the run directory: this will produce subdirectories for each temperature, then run the Markov chain and draw samples from PDF at each temperatures: the samples are saved to files in the subdirectories;
5. Run the python script analyse_samples_Metrop_HXY.py (BASH COMMAND: python analyse_samples_Metrop_HXY.py)
   to produce the estimates of the quantities of interest, which are printed to a file.

---------------------------------------------------------------------------

The initial input should be of the form:

side         	(length of one side of the lattice)
therm_sweeps	(number of thermalization sweeps)
measurements	(number of chains/measurements)
Tmin             (Lowest temperature)
Tmax             (Highest temperature)
Tstep            (Number of temperature steps between Tmin and Tmax)
start			(initialize spins: 0 = ordered, 1 = random)
twist      		(directly propose TSFs via global spin twists: 0 = no; 1 = yes) 
seed      		(seed to initialize the random number generator)

---------------------------------------------------------------------------

Files used in the program:

SOURCE FILES:
make                    Compiles the following files and produces the executable Metrop_HXY.exe
Metrop_HXY.f90             Main program
measure_Metrop_HXY.f90     Draws samples when chain stops
setup_Metrop_HXY.f90           Sets PDF from initial.in; also sets all common variables in module 
random_Metrop_HXY.f90  	  RANMAR lagged Fibonacci random number generator
global_twist_HXY.f90        Contains the global_twist subroutine, which proposes TSFs via global phase twists

ANALYSIS FILES (IN THE RUN DIRECTORY):
initial.in                  Sets lattice size, chain length, etc.
analyse_samples_Metrop_HXY.py Computes estimates of estimators from the saved Markov chain
MCMC_functions              External functions used in analyse_markov_chain.py

---------------------------------------------------------------------------

An optional change that also measures the TSFs by finding the number of global twists that minimizes the potential (USING AS DEFAULT ON 16 FEB 2018):

measure_Metrop_HXY.f90_with_twist_TSF_measure.f90
setup_Metrop_HXY_with_twist_TSF_measure.f90