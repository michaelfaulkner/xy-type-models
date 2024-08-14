vpath %.f90 ../ ../../ ../../../ # asks make to look for .f90 files in the two directories above the current one (as well as in the current one)

CFLAGS = gfortran -c -O
OFLAGS = gfortran -O -o

.SUFFIXES:
.SUFFIXES: .o .f90
.f90.o:
	$(CFLAGS) $<

all: create-executable clean move-executable

OBJECTS = random_number_generator.o xy_ecmc_variables.o xy_ecmc_read_config_file.o pre_simulation_processes.o \
			xy_ecmc_do_checkpointing.o xy_ecmc_get_checkpoint.o delete_last_checkpoint.o trim_existing_sample_file.o \
			xy_models_record_running_sums.o xy_models_set_running_sums.o get_monte_carlo_error.o \
			xy_models_get_topological_sector_component.o xy_models_output_summary_statistics.o setup_periodic_boundaries.o \
			xy_initialise_field_configuration.o xy_create_sample_files.o xy_process_sample.o \
			attempt_external_global_moves.o xy_attempt_single_external_global_move.o get_spin_difference.o \
			calculate_emergent_field.o xy_get_potential_difference.o xy_anneal_sweep.o reset_event_counters.o \
			single_xy_event_chain.o ecmc_algorithm.o

create-executable: $(OBJECTS)
	$(OFLAGS) xy_ecmc_algorithm.exe $(OBJECTS) -lm

clean:
	rm *.o *.mod

move-executable:
	mkdir -p ../../../../executables
	mv *.exe ../../../../executables