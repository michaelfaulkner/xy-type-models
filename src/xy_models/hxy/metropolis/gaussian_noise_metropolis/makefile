vpath %.f90 ../ ../../ ../../../ ../../../../ # asks make to look for .f90 files in the two directories above the current one (as well as in the current one)

CFLAGS= gfortran -c -O
OFLAGS= gfortran -O -o

.SUFFIXES :
.SUFFIXES : .o .f90
.f90.o :
	$(CFLAGS) $<

all: create-executable clean move-executable

OBJECTS = random_number_generator.o hxy_metropolis_variables.o hxy_metropolis_read_config_file.o \
			randomise_array_of_sites.o pre_simulation_processes.o hxy_metropolis_do_checkpointing.o \
			hxy_metropolis_get_checkpoint.o delete_last_checkpoint.o trim_existing_sample_file.o \
			xy_models_record_running_sums.o xy_models_set_running_sums.o get_monte_carlo_error.o \
			get_xy_models_topological_sector_component.o xy_models_output_summary_statistics.o \
			setup_periodic_boundaries.o hxy_initialise_field_configuration.o hxy_create_sample_files.o \
			xy_models_reset_metropolis_acceptance_counters.o hxy_process_sample.o \
			xy_models_output_metropolis_acceptance_rates.o attempt_external_global_moves.o \
			hxy_attempt_single_external_global_move.o get_spin_difference.o calculate_emergent_field.o \
			hxy_gaussian_noise_metropolis_sweep.o metropolis_algorithm.o

create-executable: $(OBJECTS)
	$(OFLAGS) hxy_gaussian_noise_metropolis_algorithm.exe $(OBJECTS) -lm

clean:
	rm *.o *.mod

move-executable:
	mkdir -p ../../../../../executables
	mv *.exe ../../../../../executables
