vpath %.f90 ../ ../../ ../../../ ../../../../ # asks make to look for .f90 files in the two directories above the current one (as well as in the current one)

CFLAGS = gfortran -c -O
OFLAGS = gfortran -O -o

.SUFFIXES:
.SUFFIXES: .o .f90
.f90.o:
	$(CFLAGS) $<

all: create-executable clean move-executable

OBJECTS = random_number_generator.o 3dxy_variables.o 3dxy_read_config_file.o \
			randomise_array_of_sites.o setup_periodic_boundaries.o 3dxy_initialise_field_configuration.o \
			3dxy_create_sample_files.o 3dxy_get_and_print_observation.o 3dxy_gaussian_noise_metropolis_sweep.o \
			3dxy_output_metropolis_acceptance_rates.o 3dxy_reset_metropolis_acceptance_counters.o \
			3dxy_attempt_single_external_global_move.o pre_simulation_processes.o attempt_external_global_moves.o \
			metropolis_algorithm.o

create-executable: $(OBJECTS)
	$(OFLAGS) 3dxy_gaussian_noise_metropolis_algorithm.exe $(OBJECTS) -lm

clean:
	rm *.o *.mod

move-executable:
	mkdir -p ../../executables
	mv *.exe ../../executables
