vpath %.f90 ../ ../../

CFLAGS = gfortran -c -O
OFLAGS = gfortran -O -o

.SUFFIXES:
.SUFFIXES: .o .f90
.f90.o:
	$(CFLAGS) $<

all: create-executable clean move-executable

SETUP = setup_ECMC_HXY.o
MEASURE = measure_ECMC_HXY.o
RANDOM = random.o
TWIST = attempt_external_hxy_twist.o
main-program = hxy_ecmc_algorithm.o

XY_OBJECT_FILES = $(SETUP) $(MEASURE) $(RANDOM) $(TWIST) $(initialise-system) $(main-program)

create-executable: $(XY_OBJECT_FILES)
	$(OFLAGS) hxy_ecmc_algorithm.exe $(XY_OBJECT_FILES) -lm

clean:
	rm *.o *.mod

move-executable:
	mv *.exe ../../..