# graphphylo
visualization of phylogenetic output files

## plot_phylobayes_traces
R script to plot several parameters from the [phylobayes](https://github.com/bayesiancook/pbmpi) trace files. This can be run on a pair of trace files, one from each chain. The output file is a `pdf` automatically named based on the name of the first file. So a hypothetical pair of files `matrix_chain_1.trace` and `matrix_chain_2.trace` will produce a file called `matrix_chain1.trace.pdf`.

`Rscript plot_phylobayes_traces.R matrix_chain_1.trace matrix_chain_2.trace`

By default, this will plot up to the number of iterations (plus 0.1 * the length of the dataset) where there is a value of log(L) that is 99% of the maximum value. There is typically very little change after this point, but the plot can be forced to an arbitrary number of iterations using a third term in the command line:

`Rscript plot_phylobayes_traces.R matrix_chain_1.trace matrix_chain_2.trace 10000`

![Day_CAT_GTR1.trace.png](https://github.com/wrf/graphphylo/blob/master/Day_CAT_GTR1.trace.png)

The [Dayhoff recoded dataset](https://bitbucket.org/bzxdp/feuda_et_al_2017) from [Whelan et al 2017](https://www.nature.com/articles/s41559-017-0331-3): the chains were run for 25000 iterations, though there is essentially no change after 250, meaning it is potentially run 100x longer than it needed to be.
