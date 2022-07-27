# Note on preliminary Cramér-von Mises configuration files

Note that these are the configuration files for the preliminary Cramér-von Mises simulations.  They were previously 
called by `make_cvm_figs.py` but this is no longer the case as their location no longer matches that called within 
`make_cvm_figs.py` (which now calls the configuration files for the Cramér-von Mises proper simulations).

Note also that each configuration file whose name is suffixed with a letter is a Markov-process configuration file, 
where the configuration file with no suffix must be used for the sample analysis.  For example, `64x64_metrop_a.txt`, 
`64x64_metrop_b.txt`, `64x64_metrop_c.txt`, `64x64_metrop_d.txt` and `64x64_metrop_e.txt` are the five Markov-process 
configuration files that generate the samples for `64x64_metrop.txt`.  This allows us to submit four separate 
Markov-process jobs to the cluster, thus optimising simulation time and ensuring larger system sizes fall within the 
time limit. 
