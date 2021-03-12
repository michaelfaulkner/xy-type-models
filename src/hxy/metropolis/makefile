vpath %.f90 ../ ../../ # asks make to look for .f90 files in the two directories above the current one (as well as in the current one)

CFLAGS= gfortran -c -O
OFLAGS= gfortran -O -o

.SUFFIXES :
.SUFFIXES : .o .f90
.f90.o :
	$(CFLAGS) $<

all: create-executable clean move-executable

OBJECTS = random.o setup_hxy_metropolis_algorithm.o initialise_hxy_system_and_sample_files.o hxy_sampler.o \
			 attempt_external_hxy_twist.o hxy_metropolis_algorithm.o

create-executable: $(OBJECTS)
	$(OFLAGS) hxy_metropolis_algorithm.exe $(OBJECTS) -lm

clean:
	rm *.o *.mod

move-executable:
	mv *.exe ../../..
