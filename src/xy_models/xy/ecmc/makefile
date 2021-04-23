vpath %.f90 ../ ../../ ../../../ # asks make to look for .f90 files in the two directories above the current one (as well as in the current one)

CFLAGS = gfortran -c -O
OFLAGS = gfortran -O -o

.SUFFIXES:
.SUFFIXES: .o .f90
.f90.o:
	$(CFLAGS) $<

all: create-executable clean move-executable

OBJECTS = random_number_generator.o xy_ecmc_variables.o xy_ecmc_read_config_file.o periodic_boundaries.o \
			xy_models_initialise_field_configuration.o xy_models_create_sample_file.o xy_draw_and_print_observation.o \
			xy_attempt_external_global_move.o single_xy_event_chain.o pre_simulation_processes.o ecmc_algorithm.o

create-executable: $(OBJECTS)
	$(OFLAGS) xy_ecmc_algorithm.exe $(OBJECTS) -lm

clean:
	rm *.o *.mod

move-executable:
	mv *.exe ../../../..