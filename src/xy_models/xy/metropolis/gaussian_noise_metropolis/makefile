vpath %.f90 ../ ../../ ../../../ ../../../../ # asks make to look for .f90 files in the two directories above the current one (as well as in the current one)

CFLAGS = gfortran -c -O
OFLAGS = gfortran -O -o

.SUFFIXES:
.SUFFIXES: .o .f90
.f90.o:
	$(CFLAGS) $<

all: create-executable clean move-executable

OBJECTS = random_number_generator.o xy_metropolis_variables.o xy_metropolis_read_config_file.o \
			randomise_array_of_sites.o pre_simulation_processes.o xy_metropolis_do_checkpointing.o \
			xy_metropolis_get_checkpoint.o delete_last_checkpoint.o trim_existing_sample_file.o \
			setup_periodic_boundaries.o xy_initialise_field_configuration.o xy_create_sample_files.o \
			xy_models_reset_metropolis_acceptance_counters.o xy_get_and_print_observation.o \
			xy_models_output_metropolis_acceptance_rates.o attempt_external_global_moves.o \
			xy_attempt_single_external_global_move.o get_spin_difference.o xy_get_potential_difference.o \
			xy_anneal_sweep.o xy_gaussian_noise_metropolis_sweep.o metropolis_algorithm.o

create-executable: $(OBJECTS)
	$(OFLAGS) xy_gaussian_noise_metropolis_algorithm.exe $(OBJECTS) -lm

clean:
	rm *.o *.mod

move-executable:
	mkdir -p ../../../../../executables
	mv *.exe ../../../../../executables
