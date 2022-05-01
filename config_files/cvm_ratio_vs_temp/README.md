# Note on Cramér-von Mises configuration files

Note that each configuration file whose name is suffixed with a letter is a Markov-process configuration file, where 
the configuration file with no suffix must be used for the sample analysis.  For example, `64x64_metrop_local_a.txt`, 
`64x64_metrop_local_b.txt`, `64x64_metrop_local_c.txt` and `64x64_metrop_local_d.txt` are the four Markov-process 
configuration files that generate the samples for `64x64_metrop_local.txt`.  This allows us to submit four separate 
Markov-process jobs to the cluster, thus optimising simulation time and ensuring larger system sizes fall within the 
time limit. 
