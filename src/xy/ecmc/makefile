vpath %.f90 ../ ../../ # asks make to look for .f90 files in the two directories above the current one (as well as in the current one)

CFLAGS = gfortran -c -O
OFLAGS = gfortran -O -o

.SUFFIXES:
.SUFFIXES: .o .f90
.f90.o:
	$(CFLAGS) $<

all: create-executable clean move-executable

SETUP = setup_xy_ecmc_algorithm.o
MEASURE = get_xy_observation.o
RANDOM = random.o
TWIST = attempt_external_xy_twist.o
initialise-system = initialise_xy_system_and_sample_files.o
main-program = xy_ecmc_algorithm.o

XY_OBJECT_FILES = $(SETUP) $(MEASURE) $(RANDOM) $(TWIST) $(initialise-system) $(main-program)

create-executable: $(XY_OBJECT_FILES)
	$(OFLAGS) xy_ecmc_algorithm.exe $(XY_OBJECT_FILES) -lm

clean:
	rm *.o *.mod

move-executable:
	mv *.exe ../../..