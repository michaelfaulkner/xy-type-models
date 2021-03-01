CFLAGS= gfortran -c -O
OFLAGS= gfortran -O -o

.SUFFIXES :
.SUFFIXES : .o .f90
.f90.o :
	$(CFLAGS) $<

all : GLE

SETUP=setup_metrop_GLE.o
MEASURE=measure_metrop_GLE.o
RAND=random.o
TWIST=global_TSF_GLE.o
MARKOVCHAIN=metrop_GLE.o

GLE= $(SETUP) $(MEASURE) $(RAND) $(TWIST) $(MARKOVCHAIN)

GLE : $(GLE)
	$(OFLAGS) metrop_GLE.exe $(GLE) -lm
