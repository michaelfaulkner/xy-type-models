# Note on the topological-susceptibility configuration files

Note that each configuration file whose name is suffixed with a letter is a Markov-process configuration file, where 
the configuration file with no suffix must be used for the sample analysis.  For example, `32x32_hxy_all_a.txt` and 
`32x32_hxy_all_a.txt` are the two Markov-process configuration files that generate the summary-statistic files for 
`32x32_hxy_all.txt`.  This allows us to submit two separate Markov-process runs to the cluster, thus optimising 
simulation time and ensuring larger system sizes fall within the time limit. 
