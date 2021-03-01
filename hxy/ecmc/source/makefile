CFLAGS= gfortran -c -O
OFLAGS= gfortran -O -o

.SUFFIXES :
.SUFFIXES : .o .f90
.f90.o :
	$(CFLAGS) $<

all : HXY

# SETUP=setup_ECMC_HXY.o
SETUP=setup_ECMC_HXY_with_twist_TSF_measure.o
# MEASURE=measure_ECMC_HXY.o
MEASURE=measure_ECMC_HXY_with_twist_TSF_measure.o
RAND=random.o
TWIST=global_twist_HXY.o
EVENTCHAIN=ECMC_HXY.o

HXY= $(SETUP) $(MEASURE) $(RAND) $(TWIST) $(EVENTCHAIN)

HXY : $(HXY)
	$(OFLAGS) ECMC_HXY.exe $(HXY) -lm
