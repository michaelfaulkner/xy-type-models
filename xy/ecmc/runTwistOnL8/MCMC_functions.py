import numpy as np

# CALL Rpy2 LIBRARIES
import rpy2.robjects.numpy2ri as n2ri
import rpy2.robjects.packages as rpackages
n2ri.activate()                                                                 # ACTIVATE numpy - R INTERFACE; MAKES numpy VECTORS SUITABLE FOR R

mcmcse = rpackages.importr('mcmcse')                                            # IMPORT MCMC STANDARD-ERROR PACKAGE FOR ESTIMATING ESTIMATOR WITH A STANDARD ERROR
def MCMC_MEAN_VAR(data):
    ExpAndStdErr = np.asarray(mcmcse.mcse(data))                                # USE STANDARD-ERROR PACKAGE; np.asarray TRANSFORMS R VECTOR TO PYTHON
    mean = float(ExpAndStdErr[0])
    err_mean = float(ExpAndStdErr[1])
    ExpAndStdErr = np.asarray(mcmcse.mcse((data - mean)**2))
    var = float(ExpAndStdErr[0])
    err_var = float(ExpAndStdErr[1])
    output = (mean,err_mean,var,err_var)
    return output

def autocorrelator(data, points=None):
    n = len(data)
    f = np.fft.fft (np.hstack ([data, np.zeros(n)]))
    aco = np.fft.ifft (f * np.conj (f))
    aco = np.real(aco)
    aco = aco[:(n+1) // 2]
    aco /= range (n, n // 2, -1)
    if points is not None and points < len (aco):
        return aco[:points]
    else:
        return aco
