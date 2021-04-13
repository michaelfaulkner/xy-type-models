vpath %.f90 ../ ../../ # asks make to look for .f90 files in the two directories above the current one (as well as in the current one)

CFLAGS= gfortran -c -O
OFLAGS= gfortran -O -o

.SUFFIXES :
.SUFFIXES : .o .f90
.f90.o :
	$(CFLAGS) $<

all: create-executable clean move-executable

OBJECTS = random_number_generator.o setup_maggs_electrolyte_algorithm.o randomise_array_of_sites.o \
			periodic_boundaries.o maggs_electrolyte_initialise_field_configuration.o \
			maggs_electrolyte_create_sample_files.o maggs_electrolyte_sampler.o \
			attempt_external_topological_sector_fluctuation.o  multivalued_charge_maggs_electrolyte_metropolis_sweep.o \
			maggs_electrolyte_output_metropolis_acceptance_rates.o \
			maggs_electrolyte_reset_metropolis_acceptance_counters.o metropolis_algorithm.o

create-executable : $(OBJECTS)
	$(OFLAGS) multivalued_charge_maggs_electrolyte_algorithm.exe $(OBJECTS) -lm

clean:
	rm *.o *.mod

move-executable:
	mv *.exe ../../..
