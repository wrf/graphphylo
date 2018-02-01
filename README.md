# graphphylo
visualization of phylogenetic output files

## plot_phylobayes_traces
R script to plot several parameters from the [phylobayes](https://github.com/bayesiancook/pbmpi) trace files. This can be run on a pair of trace files, one from each chain. The output file is a `pdf` automatically named based on the name of the first file. So a hypothetical pair of files `matrix_chain_1.trace` and `matrix_chain_2.trace` will produce a file called `matrix_chain_1.trace.pdf`.

`Rscript plot_phylobayes_traces.R matrix_chain_1.trace matrix_chain_2.trace`

By default, this will plot up to the number of iterations (plus 0.1 * the length of the dataset) where there is a value of log(L) that is 99% of the maximum value. There is typically very little change after this point, but the plot can be forced to an arbitrary number of iterations using a third term in the command line:

`Rscript plot_phylobayes_traces.R matrix_chain_1.trace matrix_chain_2.trace 10000`

![Day_CAT_GTR1.trace.png](https://github.com/wrf/graphphylo/blob/master/Day_CAT_GTR1.trace.png)

The [Dayhoff recoded dataset](https://bitbucket.org/bzxdp/feuda_et_al_2017) from [Whelan et al 2017](https://www.nature.com/articles/s41559-017-0331-3): the chains were run for 25000 iterations, though there is essentially no change after 250, meaning it is potentially run 100x longer than it needed to be.

Running `bpcomp` on the both chains gives exactly the same values reported in the paper:

```
$ ./bpcomp -x 7000 Day_CAT_GTR1 Day_CAT_GTR2

initialising random
seed was : 916191


Day_CAT_GTR1.treelist : 18853 trees
Day_CAT_GTR2.treelist : 18851 trees

maxdiff     : 0.0760677
meandiff    : 0.00201292
```

Although when `bpcomp` is run on a much narrower set of the chains, the values suggest that the chains have not converged. Nonetheless, it is clear that there is very little activity in the program after 500 iterations. It would seem likely based on the differences in log-likelihood that the only changes after this point are small modifications in branch order of minor clades. Thus, despite the _guidelines_ in the manual of the program, it is unclear what mathematical value there is in running the program any longer.

```
bpcomp -x 250 1 2500 Day_CAT_GTR1 Day_CAT_GTR2

initialising random
seed was : 716406


Day_CAT_GTR1.treelist : 2250 trees
Day_CAT_GTR2.treelist : 2250 trees

maxdiff     : 1
meandiff    : 0.0211733
```
