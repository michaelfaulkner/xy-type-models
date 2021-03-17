vpath %.f90 ../ ../../ # asks make to look for .f90 files in the two directories above the current one (as well as in the current one)

CFLAGS= gfortran -c -O
OFLAGS= gfortran -O -o

.SUFFIXES :
.SUFFIXES : .o .f90
.f90.o :
	$(CFLAGS) $<

all: create-executable clean move-executable

OBJECTS = random.o setup_maggs_electrolyte_algorithm.o maggs_electrolyte_sampler.o \
			attempt_external_topological_sector_fluctuation.o elementary_charge_maggs_electrolyte_algorithm.o

create-executable : $(OBJECTS)
	$(OFLAGS) elementary_charge_maggs_electrolyte_algorithm.exe $(OBJECTS) -lm

clean:
	rm *.o *.mod

move-executable:
	mv *.exe ../../..
