vpath %.f90 ../ ../../ # asks make to look for .f90 files in the two directories above the current one (as well as in the current one)

CFLAGS= gfortran -c -O
OFLAGS= gfortran -O -o

.SUFFIXES :
.SUFFIXES : .o .f90
.f90.o :
	$(CFLAGS) $<

all: create-executable clean move-executable

OBJECTS = random_number_generator.o elementary_electrolyte_variables.o elementary_electrolyte_read_config_file.o \
			randomise_array_of_sites.o electrolyte_record_running_sums.o electrolyte_set_running_sums.o \
			get_monte_carlo_error.o electrolyte_output_summary_statistics.o setup_periodic_boundaries.o \
			electrolyte_create_sample_files.o elementary_electrolyte_initialise_field_configuration.o \
			electrolyte_process_sample.o electrolyte_attempt_single_external_global_move.o \
			elementary_electrolyte_attempt_charge_hop.o electrolyte_metropolis_sweep.o \
			electrolyte_output_metropolis_acceptance_rates.o electrolyte_reset_metropolis_acceptance_counters.o \
			pre_simulation_processes.o elementary_electrolyte_do_checkpointing.o \
			elementary_electrolyte_get_checkpoint.o delete_last_checkpoint.o trim_existing_sample_file.o \
			attempt_external_global_moves.o metropolis_algorithm.o

create-executable : $(OBJECTS)
	$(OFLAGS) elementary_electrolyte_algorithm.exe $(OBJECTS) -lm

clean:
	rm *.o *.mod

move-executable:
	mkdir -p ../../../executables
	mv *.exe ../../../executables
