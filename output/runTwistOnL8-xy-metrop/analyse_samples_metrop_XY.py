from math import *
import random
import numpy as np
import pylab
import string
import MCMC_functions as MCMC_fns
from matplotlib.backends.backend_pdf import PdfPages as pdf

# SET PRECISION OF MODULAR OPERATIONS
epsilon = 10 ** (-6)

# READ IN HYPERPARAMETERS FROM initial.in
initial = []
f = open('initial.in')
for line in f:
    var = string.split(line)[0]
    if var != 'end':
        initial.append(var)
f.close()
length = int(initial[0])                                                        # LATTICE LENGTH
measurements = int(initial[2])                                                  # No of chains/measurements
Tmin = float(initial[3])                                                        # BASE TEMPERATURE
Tmax = float(initial[4])                                                        # MAX TEMP
Tsteps = int(initial[5])                                                        # NO OF TEMP STEPS
TwistOn = int(initial[8])                                                       # TSF TWISTS MOVE ON OR OFF: 0 = OFF; 1 = ON

volume = float(length ** 2)
maxchainlength = volume * pi
temp = Tmin

if Tsteps==0:
    Tincr = 0.0
else:
    Tincr = (Tmax - Tmin) / Tsteps

# CREATE OUTPUT FILES FOR MAGNETIZATION AND TOPOLOGICAL PROPERTIES [I.E., HELICITY MOD. AND TOPOLOGICAL-SECTOR FLUCTUATIONS (TSFs)]
output1 = open('magn_metrop_XY_' + str('{0:0{width}}'.format(length,width=3)) + 'x'+str('{0:0{width}}'.format(length,width=3)) + '.dat','w')
output2 = open('topology_metrop_XY_' + str('{0:0{width}}'.format(length,width=3)) + 'x'+str('{0:0{width}}'.format(length,width=3)) + '.dat','w')
output3 = open('spec_heat_ECMC_XY_' + str('{0:0{width}}'.format(length,width=3)) + 'x'+str('{0:0{width}}'.format(length,width=3)) + '.dat','w')

if TwistOn == 1:
    output1.write("Temp.".ljust(5) + '\t' + "Magnetization".ljust(20) + '\t' + "Mag. error".ljust(20) + '\t' + "Susceptibility".ljust(20) + '\t' + "Susc. error".ljust(20) + '\t' + "Av. sweeps".ljust(20) + '\t' + "Twist accept rate".ljust(20) + '\n')
    output2.write("Temp.".ljust(5) + '\t' + "Helicity".ljust(20) + '\t' + "Hel. error".ljust(20) + '\t' + "TSF".ljust(20) + '\t' + "TSF error".ljust(20) + '\t' + "Accept rate".ljust(20) + '\t' + "Twist accept rate".ljust(20) + '\n')
    output3.write("Temp.".ljust(5) + '\t' + "Specific heat".ljust(20) + '\t' + "SH error".ljust(20) + '\t' + "Av. sweeps".ljust(20) + '\t' + "Twist accept rate".ljust(20) + '\n')
else:
    output1.write("Temp.".ljust(5) + '\t' + "Magnetization".ljust(20) + '\t' + "Mag. error".ljust(20) + '\t' + "Susceptibility".ljust(20) + '\t' + "Susc. error".ljust(20) + '\t' + "Accept rate".ljust(20) + '\n')
    output2.write("Temp.".ljust(5) + '\t' + "Helicity".ljust(20) + '\t' + "Hel. error".ljust(20) + '\t' + "TSF".ljust(20) + '\t' + "TSF error".ljust(20) + '\t' + "Accept rate".ljust(20) + '\n')
    output3.write("Temp.".ljust(5) + '\t' + "Specific heat".ljust(20) + '\t' + "SH error".ljust(20) + '\t' + "Av. sweeps".ljust(20) + '\n')

# COMPUTE ESTIMATES OF MAGN, SUSC, AND TOPOLOGICAL PROPERTIES FROM EVENT-CHAIN SAMPLES, AS WELL AS THE ACF OF THE NORM OF MAGNETIZATION SQUARE
for i in range(Tsteps + 1):

    print temp
    beta = 1 / temp
    
    acceptanceRates = []
    
    input = open('temp_eq_'+str('{0:.2f}'.format(temp)) + '/acceptance_rates_XY_' + str('{0:0{width}}'.format(length,width=3))
    + 'x'+str('{0:0{width}}'.format(length,width=3)) + '_temp'+str('{0:.2f}'.format(temp)) + '.dat','r')
    for line in input:
        acceptanceRates.append(float(line))
    input.close()

    acceptanceRate = acceptanceRates[0]
    if TwistOn == 1:
        twistAcceptanceRate = acceptanceRates[1]

    # READ IN DATA
    magn = []
    magnsq = []
    input = open('temp_eq_' + str('{0:.2f}'.format(temp)) + '/magn_samples_XY_'+str('{0:0{width}}'.format(length,width=3))
    + 'x'+str('{0:0{width}}'.format(length,width=3)) + '_temp'+str('{0:.2f}'.format(temp)) + '.dat','r')
    for line in input:
        magn.append(float(line))
        magnsq.append(float(line) ** 2)
    input.close()

    # CONVERT TO numpy ARRAY FOR DATA ANALYTICS IN R
    magn = np.asarray(magn)
    magnsq = np.asarray(magnsq)

    # ESTIMATE MEAN AND VARIANCE OF NORM OF MAGNETIZATION (WITH ERRORS)
    MeanAndVar = MCMC_fns.MCMC_MEAN_VAR(magn)
    magn = float(MeanAndVar[0])
    err_magn = float(MeanAndVar[1])
    susc = beta * float(MeanAndVar[2])                                          # HAD Chi ~ sqrt(Var) FOR SOME REASON; NOW CORRECTED
    err_susc = beta * float(MeanAndVar[3])

    # WRITE MAGNETIZATION DATA TO FILE
    if TwistOn == 1:
        output1.write(str(temp).ljust(5) + '\t' + str(magn).ljust(20) + '\t' + str(err_magn).ljust(20)
        + '\t' + str(susc).ljust(20) + '\t' + str(err_susc).ljust(20) + '\t' + str(acceptanceRate).ljust(20) + '\t' + str(twistAcceptanceRate).ljust(20) + '\n')
    else:
        output1.write(str(temp).ljust(5) + '\t' + str(magn).ljust(20) + '\t' + str(err_magn).ljust(20)
        + '\t' + str(susc).ljust(20) + '\t' + str(err_susc).ljust(20) + '\t' + str(acceptanceRate).ljust(20) + '\n')

    # AUTOCORRELATION FUNCTION
    av_magnsq = np.mean(magnsq)
    aco = MCMC_fns.autocorrelator(magnsq - av_magnsq) / np.var(magnsq)
    f = open('temp_eq_'+str('{0:.2f}'.format(temp)) + '/autocorr_magnsq_metrop_XY_'+str('{0:0{width}}'.format(length,width=3))
    + 'x'+str('{0:0{width}}'.format(length,width=3)) + '_temp'+str('{0:.2f}'.format(temp)) + '.dat','w')
    for a in aco:
        f.write(str(a) + '\n')
    f.close()

# COMMENT THE FOLLOWING BACK IN TO PLOT THE AUTOCORRELATION FUNCTION

#    pylab.title('Metropolis XY autocorrelation' + '\n' + 'L = ' + str(length) + '; T = ' + str(temp) +
#    '; $ \\ell =$' + str(maxchainlength) + ';' + '\n' + 'sweeps / measurement  = '
#    + str(acceptanceRate) + '; measurements = ' + str(measurements),  fontsize=16)
#    pylab.xlabel('time in sweeps (events / spin)', fontsize=18)
#    pylab.ylabel('$M^2 -  M^2$ autocorrelation', fontsize=18)
#    XLimit = np.argmax(aco < 0.02)
#    XScale = np.array(range(XLimit)) * AverageSweeps
#    pylab.semilogy(XScale, aco[0:XLimit], marker='.', markersize=15, color='k', linestyle='None')
#    pylab.tight_layout()
#    pylab.ylim((0.02,1.1))
#    pp = pdf('temp_eq_'+str('{0:.2f}'.format(temp))+'/autocorr_magnsq_XY_'+str('{0:0{width}}'.format(length,width=3))+'x'+str('{0:0{width}}'.format(length,width=3))+'_temp'+str('{0:.2f}'.format(temp))+'.pdf')
#    pp.savefig(bbox_inches='tight')
#    pp.close()
#    pylab.clf()

    # CALCULATE TOPOLOGICAL PROPERTIES: HELICITY MODULUS AND TSFs

    # READ IN DATA
    cos_top = [[] for i in range(2)]
    sin_top = [[] for i in range(2)]
    Ebar = [[] for i in range(2)]
    Potential = []
    
    input1 = open('temp_eq_'+str('{0:.2f}'.format(temp))+'/cos_top_x_samples_XY_'+str('{0:0{width}}'.format(length,width=3))+'x'+str('{0:0{width}}'.format(length,width=3))+'_temp'+str('{0:.2f}'.format(temp))+'.dat','r')
    input2 = open('temp_eq_'+str('{0:.2f}'.format(temp))+'/sin_top_x_samples_XY_'+str('{0:0{width}}'.format(length,width=3))+'x'+str('{0:0{width}}'.format(length,width=3))+'_temp'+str('{0:.2f}'.format(temp))+'.dat','r')
    input3 = open('temp_eq_'+str('{0:.2f}'.format(temp))+'/cos_top_y_samples_XY_'+str('{0:0{width}}'.format(length,width=3))+'x'+str('{0:0{width}}'.format(length,width=3))+'_temp'+str('{0:.2f}'.format(temp))+'.dat','r')
    input4 = open('temp_eq_'+str('{0:.2f}'.format(temp))+'/sin_top_y_samples_XY_'+str('{0:0{width}}'.format(length,width=3))+'x'+str('{0:0{width}}'.format(length,width=3))+'_temp'+str('{0:.2f}'.format(temp))+'.dat','r')
    input5 = open('temp_eq_'+str('{0:.2f}'.format(temp))+'/Ebar_x_samples_XY_'+str('{0:0{width}}'.format(length,width=3))+'x'+str('{0:0{width}}'.format(length,width=3))+'_temp'+str('{0:.2f}'.format(temp))+'.dat','r')
    input6 = open('temp_eq_'+str('{0:.2f}'.format(temp))+'/Ebar_y_samples_XY_'+str('{0:0{width}}'.format(length,width=3))+'x'+str('{0:0{width}}'.format(length,width=3))+'_temp'+str('{0:.2f}'.format(temp))+'.dat','r')
    input7 = open('temp_eq_'+str('{0:.2f}'.format(temp))+'/potential_samples_XY_'+str('{0:0{width}}'.format(length,width=3))+'x'+str('{0:0{width}}'.format(length,width=3))+'_temp'+str('{0:.2f}'.format(temp))+'.dat','r')
    
    for line in input1:
        cos_top[0].append(float(line))
    input1.close()
    for line in input2:
        sin_top[0].append(float(line))
    input2.close()
    for line in input3:
        cos_top[1].append(float(line))
    input3.close()
    for line in input4:
        sin_top[1].append(float(line))
    input4.close()
    for line in input5:
        Ebar[0].append(float(line))
    input5.close()
    for line in input6:
        Ebar[1].append(float(line))
    input6.close()
    for line in input7:
        Potential.append(float(line))
    input7.close()

    # TOPOLOGICAL-SECTOR FLUCTUATIONS
    EbarW = [[] for i in range(2)]
    for i in range(measurements):
        xTwist = 0.0
        n = 1
        while True:
            if - cos(n * 2 * pi / length) * cos_top[0][i] - sin(n * 2 * pi / length) * sin_top[0][i] - ( - cos((n - 1) * 2 * pi / length) * cos_top[0][i] - sin((n - 1) * 2 * pi / length) * sin_top[0][i]) < - epsilon:
                xTwist = n * 2 * pi / length
                n += 1
            else:
                break
    
        n = 1
        while True:
            if - cos(n * 2 * pi / length) * cos_top[0][i] + sin(n * 2 * pi / length) * sin_top[0][i] - ( - cos((n - 1) * 2 * pi / length) * cos_top[0][i] + sin((n - 1) * 2 * pi / length) * sin_top[0][i]) < - epsilon:
                xTwist = - n * 2 * pi / length
                n += 1
            else:
                break

        EbarW[1].append(xTwist)

        yTwist = 0.0
        n = 1
        while True:
            if - cos(n * 2 * pi / length) * cos_top[1][i] - sin(n * 2 * pi / length) * sin_top[1][i] - ( - cos((n - 1) * 2 * pi / length) * cos_top[1][i] - sin((n - 1) * 2 * pi / length) * sin_top[1][i]) < - epsilon:
                yTwist = n * 2 * pi / length
                n += 1
            else:
                break

        n = 1
        while True:
            if - cos(n * 2 * pi / length) * cos_top[1][i] + sin(n * 2 * pi / length) * sin_top[1][i] - ( - cos((n - 1) * 2 * pi / length) * cos_top[1][i] + sin((n - 1) * 2 * pi / length) * sin_top[1][i]) < - epsilon:
                yTwist = - n * 2 * pi / length
                n += 1
            else:
                break

        EbarW[0].append(- yTwist) # MINUS SIGN TO MIMIC FOR EMERGENT-FIELD DEF: E_x(i) = + (theta(i) - theta(neg_y(i)));
                                  # E_y(i) = - (theta(i) - theta(neg_x(i))) (WITH MODULAR OPERATION)
    
    # CONVERT TO numpy ARRAYS FOR DATA ANALYTICS IN R
    cos_top = np.asarray(cos_top)
    sin_top = np.asarray(sin_top)
    Ebar = np.asarray(Ebar)
    EbarW = np.asarray(EbarW)
    Potential = np.asarray(Potential)    

    # ESTIMATE MEAN AND VARIANCE OF HELICITY MODULUS (WITH ERRORS)
    MeanAndVar1 = MCMC_fns.MCMC_MEAN_VAR(sin_top[0])
    MeanAndVar2 = MCMC_fns.MCMC_MEAN_VAR(sin_top[1])
    UpsilonSamples = 0.5 * ((cos_top[0] + cos_top[1]) - volume * beta * ((sin_top[0] - float(MeanAndVar1[0])) ** 2 + (sin_top[1] - float(MeanAndVar2[0])) ** 2))
    MeanAndVar = MCMC_fns.MCMC_MEAN_VAR(UpsilonSamples)
    Upsilon = float(MeanAndVar[0])
    err_Upsilon = float(MeanAndVar[1])    

    # ESTIMATE MEAN AND VARIANCE OF TSFs (WITH ERRORS)
    MeanAndVar1 = MCMC_fns.MCMC_MEAN_VAR(EbarW[0])
    MeanAndVar2 = MCMC_fns.MCMC_MEAN_VAR(EbarW[1])
    TSFsamples = (EbarW[0] - float(MeanAndVar1[0])) ** 2 + (EbarW[1] - float(MeanAndVar2[0])) ** 2
    MeanAndVar = MCMC_fns.MCMC_MEAN_VAR(TSFsamples)
    TSF = volume * beta * float(MeanAndVar[0])
    err_TSF = volume * beta * float(MeanAndVar[1])

    # WRITE DATA TO FILE
    if TwistOn == 1:
        output2.write(str(temp).ljust(5) + '\t' + str(Upsilon).ljust(20) + '\t' + str(err_Upsilon).ljust(20) + '\t' + str(TSF).ljust(20) + '\t' + str(err_TSF).ljust(20) + '\t' + str(acceptanceRate).ljust(20) + '\t' + str(twistAcceptanceRate).ljust(20) + '\n')
    else:
        output2.write(str(temp).ljust(5) + '\t' + str(Upsilon).ljust(20) + '\t' + str(err_Upsilon).ljust(20) + '\t' + str(TSF).ljust(20) + '\t' + str(err_TSF).ljust(20) + '\t' + str(acceptanceRate).ljust(20) + '\n')

    # SPECIFIC HEAT
    MeanAndVar = MCMC_fns.MCMC_MEAN_VAR(Potential)
    SpecificHeatSamples = (Potential - float(MeanAndVar[0])) ** 2
    MeanAndVar = MCMC_fns.MCMC_MEAN_VAR(SpecificHeatSamples)
    SpecificHeat = beta * beta * float(MeanAndVar[0]) / volume
    err_SpecificHeat = beta * beta * float(MeanAndVar[1]) / volume

    # WRITE DATA TO FILE
    if TwistOn == 1:
        output3.write(str(temp).ljust(5) + '\t' + str(SpecificHeat).ljust(20) + '\t' + str(err_SpecificHeat).ljust(20) + '\t' + str(acceptanceRate).ljust(20) + '\t' + str(twistAcceptanceRate).ljust(20) + '\n')
    else:
        output3.write(str(temp).ljust(5) + '\t' + str(SpecificHeat).ljust(20) + '\t' + str(err_SpecificHeat).ljust(20) + '\t' + str(acceptanceRate).ljust(20) + '\n')

    temp += Tincr

output1.close()
output2.close()
