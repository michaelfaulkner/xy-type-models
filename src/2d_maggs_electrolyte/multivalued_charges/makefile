vpath %.f90 ../ ../../ # asks make to look for .f90 files in the two directories above the current one (as well as in the current one)

CFLAGS= gfortran -c -O
OFLAGS= gfortran -O -o

.SUFFIXES :
.SUFFIXES : .o .f90
.f90.o :
	$(CFLAGS) $<

all: create-executable clean move-executable

SETUP=setup_maggs_electrolyte_algorithm.o
MEASURE=maggs_electrolyte_sampler.o
RAND=random.o
TWIST=global_TSF_GLE.o
main-program = multivalued_charge_maggs_electrolyte_algorithm.o

OBJECT_FILES = $(SETUP) $(MEASURE) $(RAND) $(TWIST) $(main-program)

create-executable : $(OBJECT_FILES)
	$(OFLAGS) multivalued_charge_maggs_electrolyte_algorithm.exe $(OBJECT_FILES) -lm

clean:
	rm *.o *.mod

move-executable:
	mv *.exe ../../..
