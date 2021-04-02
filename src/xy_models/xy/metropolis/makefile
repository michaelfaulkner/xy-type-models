vpath %.f90 ../ ../../ ../../../ # asks make to look for .f90 files in the two directories above the current one (as well as in the current one)

CFLAGS = gfortran -c -O
OFLAGS = gfortran -O -o

.SUFFIXES:
.SUFFIXES: .o .f90
.f90.o:
	$(CFLAGS) $<

all: create-executable clean move-executable

OBJECTS = random.o setup_xy_metropolis_algorithm.o periodic_boundaries.o initialise_xy_system_and_sample_files.o \
			xy_sampler.o attempt_external_xy_twist.o xy_metropolis_algorithm.o

create-executable: $(OBJECTS)
	$(OFLAGS) xy_metropolis_algorithm.exe $(OBJECTS) -lm

clean:
	rm *.o *.mod

move-executable:
	mv *.exe ../../..