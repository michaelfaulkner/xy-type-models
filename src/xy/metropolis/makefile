vpath %.f90 ../ ../../

CFLAGS = gfortran -c -O
OFLAGS = gfortran -O -o

.SUFFIXES:
.SUFFIXES: .o .f90
.f90.o:
	$(CFLAGS) $<

all: create-executable clean move-executable

SETUP = setup_xy_metropolis_algorithm.o
MEASURE = get_xy_observation.o
RANDOM = random.o
TWIST = attempt_external_xy_twist.o
output-acceptance-rates = output_xy_metropolis_acceptance_rates.o
initialise-system = initialise_xy_system_and_sample_files.o
main-program = xy_metropolis_algorithm.o

XY_OBJECT_FILES = $(SETUP) $(MEASURE) $(RANDOM) $(TWIST) $(output-acceptance-rates) $(initialise-system) $(main-program)

create-executable: $(XY_OBJECT_FILES)
	$(OFLAGS) xy_metropolis_algorithm.exe $(XY_OBJECT_FILES) -lm

clean:
	rm *.o *.mod

move-executable:
	mv *.exe ../../..