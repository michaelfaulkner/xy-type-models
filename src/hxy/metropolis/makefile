vpath %.f90 ../ ../../ # asks make to look for .f90 files in the two directories above the current one (as well as in the current one)

CFLAGS= gfortran -c -O
OFLAGS= gfortran -O -o

.SUFFIXES :
.SUFFIXES : .o .f90
.f90.o :
	$(CFLAGS) $<

all: create-executable clean move-executable

SETUP = setup_hxy_metropolis_algorithm.o
MEASURE = hxy_sampler.o
RANDOM = random.o
TWIST = attempt_external_hxy_twist.o
initialise-system = initialise_hxy_system_and_sample_files.o
main-program = hxy_metropolis_algorithm.o

OBJECT_FILES = $(SETUP) $(MEASURE) $(RANDOM) $(TWIST) $(initialise-system) $(main-program)

create-executable : $(OBJECT_FILES)
	$(OFLAGS) hxy_metropolis_algorithm.exe $(OBJECT_FILES) -lm

clean:
	rm *.o *.mod

move-executable:
	mv *.exe ../../..
