import numpy as np
import rpy2.robjects.numpy2ri as n2ri
import rpy2.robjects.packages as r_packages
n2ri.activate()
laplaces_demon_r_package = r_packages.importr('LaplacesDemon')


def get_effective_sample_size(sample):
    return np.array(laplaces_demon_r_package.ESS(sample))[0]


def get_thinned_sample(sample, thinning_level):
    sample_indices_to_keep = np.array([i for i in range(len(sample)) if i % thinning_level == 0])
    return np.take(sample, sample_indices_to_keep)
