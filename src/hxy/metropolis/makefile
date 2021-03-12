vpath %.f90 ../ ../../ # asks make to look for .f90 files in the two directories above the current one (as well as in the current one)

CFLAGS= gfortran -c -O
OFLAGS= gfortran -O -o

.SUFFIXES :
.SUFFIXES : .o .f90
.f90.o :
	$(CFLAGS) $<

all: create-executable clean move-executable

SETUP=setup_metrop_HXY.o
MEASURE=measure_metrop_HXY.o
RAND=random.o
TWIST=global_twist_HXY.o
main-program=metrop_HXY.o

OBJECT_FILES= $(SETUP) $(MEASURE) $(RAND) $(TWIST) $(main-program)

create-executable : $(OBJECT_FILES)
	$(OFLAGS) metrop_HXY.exe $(OBJECT_FILES) -lm

clean:
	rm *.o *.mod

move-executable:
	mv *.exe ../../..
