import numpy as np

try:
    import rpy2.robjects.numpy2ri as n2ri
    import rpy2.robjects.packages as r_packages
    n2ri.activate()
    laplaces_demon_r_package = r_packages.importr('LaplacesDemon')
    mcmcse_r_package = r_packages.importr('mcmcse')


    def get_effective_sample_size(one_dimensional_sample):
        if len(np.atleast_2d(one_dimensional_sample)) > 1:
            raise Exception("Error: the sample passed to markov_chain_diagnostics.get_effective_sample_size() must be "
                            "one (Cartesian) dimensional.")
        return np.array(laplaces_demon_r_package.ESS(one_dimensional_sample))[0]


    def get_sample_mean_and_error(one_dimensional_sample):
        if len(np.atleast_2d(one_dimensional_sample)) > 1:
            raise Exception("Error: the sample passed to markov_chain_diagnostics.get_sample_mean_and_error() must be "
                            "one (Cartesian) dimensional.")
        sample_mean_and_error = np.array(mcmcse_r_package.mcse(one_dimensional_sample))
        return sample_mean_and_error[0, 0], sample_mean_and_error[1, 0]

except (ModuleNotFoundError, ValueError) as _:
    def get_effective_sample_size(_):
        print("rpy2 not available: get_effective_sample_size() returns None.")
        return None


    def get_sample_mean_and_error(sample):
        """nb, ddof=1 in np.std() uses 1 / (len(sample) - 1) factor"""
        return np.mean(sample), np.std(sample, ddof=1) / len(sample) ** 0.5


def get_thinned_sample(one_dimensional_sample, thinning_level):
    if len(np.atleast_2d(one_dimensional_sample)) > 1:
        raise Exception("Error: the sample passed to markov_chain_diagnostics.get_thinned_sample() must be one "
                        "(Cartesian) dimensional.")
    sample_indices_to_keep = np.array([i for i in range(len(one_dimensional_sample)) if i % thinning_level == 0])
    return np.take(one_dimensional_sample, sample_indices_to_keep)


def get_cumulative_distribution(one_dimensional_sample):
    if len(np.atleast_2d(one_dimensional_sample)) > 1:
        raise Exception("Error: the sample passed to markov_chain_diagnostics.get_cumulative_distribution() must be "
                        "one (Cartesian) dimensional.")
    """alternative calculation commented out, nb, factor of 1 / 10 (in bins) may not be optimal"""
    """count, bins_count = np.histogram(magnetisation_phase, bins=int(len(one_dimensional_sample) / 10))
    cdf = np.array([bins_count[1:], np.cumsum(count / sum(count))])"""
    bin_values = np.arange(1, len(one_dimensional_sample) + 1) / float(len(one_dimensional_sample))
    ordered_sample = np.sort(one_dimensional_sample)
    return [ordered_sample, bin_values]
