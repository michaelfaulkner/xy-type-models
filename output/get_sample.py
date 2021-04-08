import numpy as np


def magnetisation(sample_directory, temperature_directory):
    return np.linalg.norm(np.loadtxt(sample_directory + temperature_directory + '/magnetisation_sample.dat',
                                     dtype=float, delimiter=','), axis=1)
