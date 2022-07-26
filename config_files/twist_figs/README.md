# Note on configuration files for twist figures

Note that each configuration file whose name is suffixed with a letter is a Markov-process configuration file, where 
the configuration file with no suffix must be used for the sample analysis.  For example, `64x64_low_temp_all_a.txt`, 
`64x64_low_temp_all_b.txt`, `64x64_low_temp_all_c.txt`, `64x64_low_temp_all_d.txt` and `64x64_low_temp_all_e.txt` are 
the five Markov-process configuration files that generate the samples for `64x64_low_temp_all.txt`.  This allows us to 
submit four separate Markov-process jobs to the cluster, thus optimising simulation time and ensuring larger system 
sizes fall within the time limit. 
