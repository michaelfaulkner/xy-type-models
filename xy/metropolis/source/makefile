CFLAGS= gfortran -c -O
OFLAGS= gfortran -O -o

.SUFFIXES :
.SUFFIXES : .o .f90
.f90.o :
	$(CFLAGS) $<

all : XY

SETUP=setup_metrop_XY.o
MEASURE=measure_metrop_XY.o
RAND=random.o
TWIST=global_twist_XY.o
EVENTCHAIN=metrop_XY.o

XY= $(SETUP) $(MEASURE) $(RAND) $(TWIST) $(EVENTCHAIN)

XY : $(XY)
	$(OFLAGS) metrop_XY.exe $(XY) -lm
