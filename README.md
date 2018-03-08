# graphphylo
visualization of phylogenetic output files

## plot_phylobayes_traces
R script to plot several parameters from the [phylobayes](https://github.com/bayesiancook/pbmpi) trace files. This can be run on a pair of trace files, one from each chain. The output file is a `pdf` automatically named based on the name of the first file. So a hypothetical pair of files `matrix_chain_1.trace` and `matrix_chain_2.trace` will produce a file called `matrix_chain_1.trace.pdf`.

`Rscript plot_phylobayes_traces.R matrix_chain_1.trace matrix_chain_2.trace`

By default, this will plot up to the number of iterations (plus 0.1 * the length of the dataset) where there is a value of log(L) that is 99% of the maximum value. There is typically very little change after this point, but the plot can be forced to an arbitrary number of iterations using a third term in the command line:

`Rscript plot_phylobayes_traces.R matrix_chain_1.trace matrix_chain_2.trace 10000`

![Day_CAT_GTR1.trace.png](https://github.com/wrf/graphphylo/blob/master/examples/Day_CAT_GTR1.trace.png)

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

## scan_treelist
Despite the above statistical uncertainty, the list of trees itself can be converted into an animated gif, showing how the tree topology becomes fixed after a short period. Most of the changes after this point are within clades, or adjustments of branch lengths. Here is another example from the re-analyses by [Feuda et al](https://bitbucket.org/bzxdp/feuda_et_al_2017), using the D20 dataset from [Whelan et al 2015](https://figshare.com/articles/Error_signal_and_the_placement_of_Ctenophora_sister_to_all_other_animals/1334306) under the CAT-GTR model in `phylobayes`. Sponges and ctenophores are both monophyletic after 50 iterations, and the entire topology is basically fixed at 100. As this is visualizing only one chain, this may imply that the chain is stuck in a local maximum.

**Note: this script requires** `Bio` **and** `matplotlib` **Python libraries**

Several steps are needed to make a similar animation. First, the `.treelist` file needs to be parsed by the `scan_treelist.py` script. This generates `.png` files for each tree in the file. As this might be thousands of trees, it is recommended to specify a maximum with `-m`, maybe a few hundred. The output images are of a rooted tree, so a rooting taxon must be given with `-r`. In this case, the fungi *Spizellomyces punctatus* is used as an outgroup to root animals (though other fungi will appear as the next branch up).

`scan_treelist.py -c examples/whelan_D20_color_scheme.tab -m 200 -r Spizellomyces_punctatus -t CAT_GTR1.treelist --title Whelan-D20-CAT-GTR-1 -d temp_figs`

![whelan_D20_CAT-GTR-1_animated_reduced.gif](https://github.com/wrf/graphphylo/blob/master/examples/whelan_D20_CAT-GTR-1_animated_reduced.gif)

An optional color scheme can be specified for any or all taxa in a file given by `-c`. This should be a tab-delimited text file, and colors given as `#RRGGBB`. Any species not specified are colored black. This is ultimately passed to `matplotlib` as a dictionary of colors by species.

```
Spizellomyces_punctatus	#cf003a
Monosiga_ovata	#9a601b
```

Next, all of the `.png` files are cropped and converted to `.gif`s using `convert` in [ImageMagick](https://www.imagemagick.org/script/index.php). This is applied to all images as a for-loop-in-the-shell. Resizing below 70% makes it difficult to read tip labels.

`for FILE in *.png; do convert $FILE -crop 640x440+80+20 +repage -resize 70% $FILE.gif ; done`

The `.gif`s are then combined into an animated tree, also using `convert`.

`convert -delay 8 -loop 0 *.gif animated_tree.gif`

A [useful post about data visualization](http://viewshed.matinic.us/2018/01/13/1139/) suggested using a program called [gifsicle](https://www.lcdf.org/gifsicle/) to further compress the `.gif`. As only 5 different colors are used (red for the outgroup, brown for choanoflagellates, green for ctenophores, and purple for sponges, and the rest black), this can be reduced in the final image. Specifying fewer than 16 colors starts to have problems with rendering the main ones.

`gifsicle -O3 animated_tree.gif --colors 16 -o animated_tree_reduced_color.gif`
