CFLAGS= gfortran -c -O
OFLAGS= gfortran -O -o

.SUFFIXES :
.SUFFIXES : .o .f90
.f90.o :
	$(CFLAGS) $<

all : XY

SETUP=setup_ECMC_XY.o
MEASURE=measure_ECMC_XY.o
RAND=random.o
TWIST=global_twist_XY.o
EVENTCHAIN=ECMC_XY.o

XY= $(SETUP) $(MEASURE) $(RAND) $(TWIST) $(EVENTCHAIN)

XY : $(XY)
	$(OFLAGS) ECMC_XY.exe $(XY) -lm
