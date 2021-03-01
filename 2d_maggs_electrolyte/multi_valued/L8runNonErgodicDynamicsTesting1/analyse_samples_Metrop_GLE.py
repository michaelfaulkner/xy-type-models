from math import *
import numpy as np
import pylab
import matplotlib.pyplot as plt
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
chargeValue = float(initial[1]) 
measurements = int(initial[3])                                                  # No of chains/measurements
Tmin = float(initial[4])                                                        # BASE TEMPERATURE
Tmax = float(initial[5])                                                        # MAX TEMP
Tsteps = int(initial[6])                                                        # NO OF TEMP STEPS
globalTSFon = int(initial[11])                                                  # TSF TWISTS MOVE ON OR OFF: 0 = OFF; 1 = ON
  
volume = float(length ** 2)
temp = Tmin
TSFunit = chargeValue * length
TSFunitOverTwo = TSFunit / 2

if Tsteps==0:
    Tincr = 0.0
else:
    Tincr = (Tmax - Tmin) / Tsteps

# CREATE OUTPUT FILES FOR MAGNETIZATION AND TOPOLOGICAL PROPERTIES [I.E., HELICITY MOD. AND TOPOLOGICAL-SECTOR FLUCTUATIONS (TSFs)]

output1 = open('topology_Metrop_GLE_' + str('{0:0{width}}'.format(length,width=3)) + 'x'+str('{0:0{width}}'.format(length,width=3)) + '.dat','w')
output2 = open('spec_heat_Metrop_GLE_' + str('{0:0{width}}'.format(length,width=3)) + 'x'+str('{0:0{width}}'.format(length,width=3)) + '.dat','w')

if globalTSFon == 1:
    output1.write("Temp.".ljust(5) + '\t' + "Inverse perm.".ljust(20) + '\t' + "IP error".ljust(20) + '\t' + "TSF".ljust(20) + '\t' + "TSF error".ljust(20) + '\t' + "Aux-field accept rate".ljust(20) + '\t' + "Charge accept rate".ljust(20) + '\t' + "Twist accept rate".ljust(20) + '\n')
    output2.write("Temp.".ljust(5) + '\t' + "Specific heat".ljust(20) + '\t' + "SH error".ljust(20) + '\t' + "Aux-field accept rate".ljust(20) + '\t' + "Charge accept rate".ljust(20) + '\t' + "Twist accept rate".ljust(20) + '\n')
else:
    output1.write("Temp.".ljust(5) + '\t' + "Inverse perm.".ljust(20) + '\t' + "IP error".ljust(20) + '\t' + "TSF".ljust(20) + '\t' + "TSF error".ljust(20) + '\t' + "Aux-field accept rate".ljust(20) + '\t' + "Charge accept rate".ljust(20) + '\n')
    output2.write("Temp.".ljust(5) + '\t' + "Specific heat".ljust(20) + '\t' + "SH error".ljust(20) + '\t' + "Aux-field accept rate".ljust(20) + '\t' + "Charge accept rate".ljust(20) + '\n')

# COMPUTE ESTIMATES OF TOPOLOGICAL AND SPECIFIC HEAT PROPERTIES FROM MARKOV-CHAIN SAMPLES

for i in range(Tsteps + 1):

    print(temp)
    beta = 1 / temp
    
    acceptanceRates = []
    
    input = open('temp_eq_'+str('{0:.2f}'.format(temp)) + '/acceptance_rates_GLE_' + str('{0:0{width}}'.format(length,width=3))
    + 'x'+str('{0:0{width}}'.format(length,width=3)) + '_temp'+str('{0:.2f}'.format(temp)) + '.dat','r')
    for line in input:
        acceptanceRates.append(float(line))
    input.close()

    chargeAcceptanceRate = acceptanceRates[0]
    auxFieldAcceptanceRate = acceptanceRates[1]
    if globalTSFon == 1:
        twistAcceptanceRate = acceptanceRates[2]

    # CALCULATE TOPOLOGICAL PROPERTIES: HELICITY MODULUS AND TSFs

    # READ IN DATA
    input1 = open('temp_eq_' + str('{0:.2f}'.format(temp)) + '/Esum_x_samples_GLE_' + str('{0:0{width}}'.format(length,width=3))
    + 'x' + str('{0:0{width}}'.format(length,width=3)) + '_temp'+str('{0:.2f}'.format(temp)) + '.dat','r')
    input2 = open('temp_eq_' + str('{0:.2f}'.format(temp)) + '/Esum_y_samples_GLE_' + str('{0:0{width}}'.format(length,width=3))
    + 'x'+str('{0:0{width}}'.format(length,width=3)) + '_temp'+str('{0:.2f}'.format(temp)) + '.dat','r')
    input3 = open('temp_eq_' + str('{0:.2f}'.format(temp)) + '/potential_samples_GLE_' + str('{0:0{width}}'.format(length,width=3))
    + 'x'+str('{0:0{width}}'.format(length,width=3)) + '_temp'+str('{0:.2f}'.format(temp)) + '.dat','r')

    # DEFINE ARRAYS AND FILL WITH SAMPLES
    Esum = [[] for i in range(2)]
    Potential = []

    for line in input1:
        Esum[0].append(float(line))
    input1.close()
    for line in input2:
        Esum[1].append(float(line))
    input2.close()
    for line in input3:
        Potential.append(float(line))
    input3.close()

    # CONVERT TO numpy ARRAY FOR DATA ANALYTICS IN R
    Esum = np.asarray(Esum)
    Potential = np.asarray(Potential)
    
    # ESTIMATE MEAN AND VARIANCE OF INVERSE ELECTRIC PERMITTIVITY (WITH ERRORS)
    MeanAndVar1 = MCMC_fns.MCMC_MEAN_VAR(Esum[0])
    MeanAndVar2 = MCMC_fns.MCMC_MEAN_VAR(Esum[1])
    R = beta * ( (Esum[0] - float(MeanAndVar1[0])) ** 2 + (Esum[1] - float(MeanAndVar2[0])) ** 2 ) / 2 / volume
    InversePermSamples = 1.0 - R
    MeanAndVar = MCMC_fns.MCMC_MEAN_VAR(InversePermSamples)
    InversePerm = float(MeanAndVar[0])
    err_InversePerm = float(MeanAndVar[1])

    MeanAndVar = MCMC_fns.MCMC_MEAN_VAR(R)
    MeanR = float(MeanAndVar[0])
    DeltaR = R - MeanR
    DeltaROverMean = DeltaR / MeanR
    DeltaRSq = DeltaR * DeltaR
    MeanAndVar = MCMC_fns.MCMC_MEAN_VAR(DeltaRSq)
    RMSDeltaR = (float(MeanAndVar[0])) ** 0.5
    DeltaROverRMS = DeltaR / RMSDeltaR

    fig = plt.hist(DeltaROverRMS, bins = 100)
    plt.title('PDF resistance fluctuations $T = $' + str('{0:.2f}'.format(temp)))
    plt.savefig('PDFresistance_temp' + str('{0:.2f}'.format(temp)) + '.png')
    plt.clf()

    """"
    Inpsect ordered DeltaROverRMS histogram
    
    plt.tight_layout(pad = 0.001)
    pp = pdf('PDFresistance_temp' + str('{0:.2f}'.format(temp)) + '.pdf')
    pp.savefig(fig, bbox_inches = 'tight')
    pp.close()
    """

    # TOPOLOGICAL-SECTOR FLUCTUATIONS
    EsumW = [[] for i in range(2)]
    EsumW[0] = ( ( Esum[0] + TSFunitOverTwo ) // TSFunit ) * TSFunit
    EsumW[1] = ( ( Esum[1] + TSFunitOverTwo ) // TSFunit ) * TSFunit

    # REMOVE TSF CONTRIBUTIONS FROM Ebar = q / 2L (SINCE IN (- q / 2L, q / 2L])
    for i in range(measurements):
        if abs(Esum[0][i] - TSFunitOverTwo) < epsilon:
            EsumW[0][i] = 0.0
        if abs(Esum[1][i] - TSFunitOverTwo) < epsilon:
            EsumW[1][i] = 0.0
            
    # ESTIMATE MEAN AND VARIANCE OF TSFs (WITH ERRORS)
    MeanAndVar1 = MCMC_fns.MCMC_MEAN_VAR(EsumW[0])
    MeanAndVar2 = MCMC_fns.MCMC_MEAN_VAR(EsumW[1])
    TSFsamples = (EsumW[0] - float(MeanAndVar1[0])) ** 2 + (EsumW[1] - float(MeanAndVar2[0])) ** 2
    MeanAndVar = MCMC_fns.MCMC_MEAN_VAR(TSFsamples)
    TSF = beta * float(MeanAndVar[0]) / volume
    err_TSF = beta * float(MeanAndVar[1]) / volume

    # WRITE DATA TO FILE
    if globalTSFon == 1:
        output1.write(str(temp).ljust(5) + '\t' + str(InversePerm).ljust(20) + '\t' + str(err_InversePerm).ljust(20) + '\t' + str(TSF).ljust(20) + '\t' + str(err_TSF).ljust(20) + '\t' + str(auxFieldAcceptanceRate).ljust(20) + '\t' + str(chargeAcceptanceRate).ljust(20) + '\t' + str(twistAcceptanceRate).ljust(20) + '\n')
    else:
        output1.write(str(temp).ljust(5) + '\t' + str(InversePerm).ljust(20) + '\t' + str(err_InversePerm).ljust(20) + '\t' + str(TSF).ljust(20) + '\t' + str(err_TSF).ljust(20) + '\t' + str(auxFieldAcceptanceRate).ljust(20) + '\t' + str(chargeAcceptanceRate).ljust(20) + '\n')

    # SPECIFIC HEAT
    MeanAndVar = MCMC_fns.MCMC_MEAN_VAR(Potential)
    SpecificHeatSamples = (Potential - float(MeanAndVar[0])) ** 2
    MeanAndVar = MCMC_fns.MCMC_MEAN_VAR(SpecificHeatSamples)
    SpecificHeat = beta * beta * float(MeanAndVar[0]) / volume
    err_SpecificHeat = beta * beta * float(MeanAndVar[1]) / volume

    # WRITE DATA TO FILE
    if globalTSFon == 1:
        output2.write(str(temp).ljust(5) + '\t' + str(SpecificHeat).ljust(20) + '\t' + str(err_SpecificHeat).ljust(20) + '\t' + str(auxFieldAcceptanceRate).ljust(20) + '\t' + str(chargeAcceptanceRate).ljust(20) + '\t' + str(twistAcceptanceRate).ljust(20) + '\n')
    else:
        output2.write(str(temp).ljust(5) + '\t' + str(SpecificHeat).ljust(20) + '\t' + str(err_SpecificHeat).ljust(20) + '\t' + str(auxFieldAcceptanceRate).ljust(20) + '\t' + str(chargeAcceptanceRate).ljust(20) + '\n')

    temp += Tincr

output1.close()
output2.close()
