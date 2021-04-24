vpath %.f90 ../ ../../ ../../../ # asks make to look for .f90 files in the two directories above the current one (as well as in the current one)

CFLAGS = gfortran -c -O
OFLAGS = gfortran -O -o

.SUFFIXES:
.SUFFIXES: .o .f90
.f90.o:
	$(CFLAGS) $<

all: create-executable clean move-executable

OBJECTS = random_number_generator.o hxy_ecmc_variables.o hxy_ecmc_read_config_file.o periodic_boundaries.o \
			hxy_initialise_field_configuration.o xy_models_create_sample_file.o calculate_emergent_field.o \
			get_spin_difference.o hxy_get_and_print_observation.o hxy_attempt_external_global_move.o \
			single_hxy_event_chain.o pre_simulation_processes.o ecmc_algorithm.o

create-executable: $(OBJECTS)
	$(OFLAGS) hxy_ecmc_algorithm.exe $(OBJECTS) -lm

clean:
	rm *.o *.mod

move-executable:
	mv *.exe ../../../..