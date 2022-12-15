vpath %.f90 ../ ../../ # asks make to look for .f90 files in the two directories above the current one (as well as in the current one)

CFLAGS= gfortran -c -O
OFLAGS= gfortran -O -o

.SUFFIXES :
.SUFFIXES : .o .f90
.f90.o :
	$(CFLAGS) $<

all: create-executable clean move-executable

OBJECTS = random_number_generator.o multivalued_electrolyte_variables.o multivalued_electrolyte_read_config_file.o \
			randomise_array_of_sites.o setup_periodic_boundaries.o electrolyte_create_sample_files.o \
			multivalued_electrolyte_initialise_field_configuration.o electrolyte_get_and_print_observation.o \
			electrolyte_attempt_single_external_global_move.o multivalued_electrolyte_attempt_charge_hop.o \
			electrolyte_metropolis_sweep.o electrolyte_output_metropolis_acceptance_rates.o \
			electrolyte_reset_metropolis_acceptance_counters.o pre_simulation_processes.o \
			multivalued_electrolyte_do_checkpointing.o multivalued_electrolyte_get_checkpoint.o \
			delete_last_checkpoint.o trim_existing_sample_file.o attempt_external_global_moves.o metropolis_algorithm.o

create-executable : $(OBJECTS)
	$(OFLAGS) multivalued_electrolyte_algorithm.exe $(OBJECTS) -lm

clean:
	rm *.o *.mod

move-executable:
	mkdir -p ../../../executables
	mv *.exe ../../../executables
