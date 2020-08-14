# graphphylo
visualization of phylogenetic output files, analyses, and substitution matrices

* [plot_phylobayes_traces](https://github.com/wrf/graphphylo#plot_phylobayes_traces) - plot parameters from the trace file
* [scan_treelist](https://github.com/wrf/graphphylo#scan_treelist) - make animation of the treelist
* [RF distance](https://github.com/wrf/graphphylo#rf_distance) - make plot of progressive RF distance
* [model visualization](https://github.com/wrf/graphphylo#model-visualization) - reference plots of model parameters
* [ancestral state](https://github.com/wrf/graphphylo#ancestral_state) - converting ancestral state probability into fasta sequence

## plot_phylobayes_traces
R script to plot several parameters from the [phylobayes](https://github.com/bayesiancook/pbmpi) trace files. This can be run on a pair of trace files, one from each chain. The output file is a `pdf` automatically named based on the name of the first file. So a hypothetical pair of files `matrix_chain_1.trace` and `matrix_chain_2.trace` will produce a file called `matrix_chain_1.trace.pdf`.

`Rscript plot_phylobayes_traces.R matrix_chain_1.trace matrix_chain_2.trace`

By default, this will plot up to the number of iterations (plus 0.1 * the length of the dataset) where there is a value of log(L) that is 99% of the maximum value. There is typically very little change after this point, but the plot can be forced to an arbitrary number of iterations using a third term in the command line:

`Rscript plot_phylobayes_traces.R matrix_chain_1.trace matrix_chain_2.trace 10000`

![Day_CAT_GTR1.trace.png](https://github.com/wrf/graphphylo/blob/master/examples/Day_CAT_GTR1.trace.png)

The [Dayhoff recoded dataset](https://bitbucket.org/bzxdp/feuda_et_al_2017) from [Whelan et al 2017](https://www.nature.com/articles/s41559-017-0331-3): the chains were run for 25000 iterations, though there is essentially no change after 250, meaning it is potentially run 100x longer than it needed to be.

Running `bpcomp` on both chains gives exactly the same values reported in the paper:

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

**Note: this script requires** `Bio` **and** `matplotlib` **Python libraries, and appears to require python-tk (TK interface)**. **This can be installed on Ubuntu with** `sudo apt install python3-tk`

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

## rf_distance
Pairwise [Robinson-Foulds distances](https://en.wikipedia.org/wiki/Robinson%E2%80%93Foulds_metric) can be quickly calculated for all trees in a file of unrooted trees using the `-f r` option in [RAxML](https://github.com/stamatak/standard-RAxML). This can work directly on the `.treelist` output file of `phylobayes`.

As this is an **all-versus-all** pairwise analysis, the computational time may increase exponentially with more trees, thus it may be useful to subset the trees, here taking the first 500 (i.e. the "normal" *burn-in* period, which is where all activity actually is):

`head -n 500 Day_CAT_GTR2.treelist > Day_CAT_GTR2.treelist_first_500`

Then specify the tree list with the option `-z`, and a model `-m` must be given (even though it appears nothing is calculated with it).

`raxmlHPC-PTHREADS-SSE3 -f r -z Day_CAT_GTR2.treelist_first_500 -n Day_CAT_GTR2_500t -m PROTGAMMALG`

The RAxML output is pairwise RF distance, one per line. It is clear from the output that a lot of changes happen in the first iterations, but not many after several hundred iterations. The data refer to tree 0 against tree 1, the RF distance (raw, 96) and the normalized distance (relative to total nodes, 0.65).

```
0 1: 96 0.657534
0 2: 114 0.780822
0 3: 116 0.794521
0 4: 124 0.849315
0 5: 130 0.890411
0 6: 130 0.890411
0 7: 132 0.904110
0 8: 132 0.904110
0 9: 138 0.945205
0 10: 138 0.945205
0 11: 144 0.986301
0 12: 144 0.986301
0 13: 146 1.000000
0 14: 146 1.000000
0 15: 146 1.000000
...
0 495: 146 1.000000
0 496: 146 1.000000
0 497: 146 1.000000
0 498: 146 1.000000
0 499: 146 1.000000
```

In the interest of directly examining the chain (i.e. seeing the approach to the optimal tree), only those trees in sequence (in the order of the chain) should be considered. That is, tree 0 against tree 1, then 1 against 2, and so forth. RF distances of trees that are in sequence can be extracted using the `filter_rfd.py` script.

`filter_rfd.py RAxML_RF-Distances.Day_CAT_GTR2_500t > filtered_Day_CAT_GTR2_RFd.txt`

```
0 1: 96 0.657534
1 2: 68 0.465753
2 3: 74 0.506849
3 4: 78 0.534247
4 5: 78 0.534247
...
```

Then plot the filtered RF distances in R, which will automatically generate a PDF with the same basename as the input file:

`Rscript graph_filtered_rfd.R filtered_Day_CAT_GTR2_RFd.txt`

![filtered_Day_CAT_GTR2_RFd.png](https://github.com/wrf/graphphylo/blob/master/examples/filtered_Day_CAT_GTR2_RFd.png)

## model visualization
The [CAT model](https://doi.org/10.1093/molbev/msh112) normally will dynamic generate substitution matrices for categories of sites, but usage of fixed substitution matrices (named C10 through C60) can also be selected. AA frequency data is included in the [phylobayes source code](https://github.com/bayesiancook/pbmpi/blob/master/sources/BiologicalSequences.h), though is easier to extract the model frequencies from the [IQtree source code](https://github.com/Cibiv/IQ-TREE/blob/master/model/modelprotein.cpp)

![cat_c10_diagrams_raw.png](https://github.com/wrf/graphphylo/blob/master/mixture_model/cat_c10_diagrams_raw.png)

The default order appears to be arbitrary. Here, they are rearranged in two ways. First, the sets with the lowest variance are put first, showing that many categories can be similar, but ultimately differ based on the inclusion of certain amino acids. Secondly, the amino acids themselves are reordered by a rough measure hydrophobicity (scores [from here](https://www.sigmaaldrich.com/life-science/metabolomics/learning-center/amino-acid-reference-chart.html)).

![cat_c10_diagrams_var.png](https://github.com/wrf/graphphylo/blob/master/mixture_model/cat_c10_diagrams_var.png)

Other mixture models have been created. Here, a plot is made showing the [two-category mixture model EX2](https://doi.org/10.1098/rstb.2008.0180), where the two substitution matrices differ by solvent exposure. The diagonal shows the difference in overall frequency, whereby polar/charged AAs are commonly exposed (in orange), and nonpolar ones are usually buried (in purple). Intersections of two AAs show the relative preference for that transition when in one of the two states. That is, when R is buried, the substitution to K is strongly preferred over any other substitution. It is likely that R-K is preserving an internal salt bridge. This also implies that charge amino acids on the surface may readily be substituted into many other polar AAs, thus making the relative preference less pronounced. The same applies to nonpolar amino acids, that generally change to only certain other nonpolar AAs when exposed. Polar AAs rarely change to nonpolar ones, regardless of solvent exposure.

![EX2_model_bur-exp_difference.png](https://github.com/wrf/graphphylo/blob/master/mixture_model/EX2_model_bur-exp_difference.png)

## ancestral state
This is done in two stages, and must use the [older version 4.1](https://megasun.bch.umontreal.ca/People/lartillot/www/download.html), not the [MPI version](https://github.com/bayesiancook/pbmpi) (as it is not implemented there). First, make a chain of sufficient length:

`~/phylobayes/phylobayes4.1c/data/pb -d Whelan_D16_Opisthokonta_reduced.phy -T d16opi_c1.con.tre -cat -gtr -x 1 200 -s d16_opi_pb_catgtr_for_anc_c1`

Then, generate the ancestral states. This creates a folder of with a lot of files named `.ancstatepostprob`, one for each leaf and each internal node.

`~/phylobayes/phylobayes4.1c/data/ancestral -x 100 1 d16_opi_pb_catgtr_for_anc_c1`

Then, generate the fasta sequence with `ancestral_probability_to_fasta.py`:

`ancestral_probability_to_fasta.py -p  d16_opi_pb_catgtr_for_anc_c1_sample_104_Aque_Scoa.ancstatepostprob > d16_opi_node_104.fasta`

```
# reading probabilities d16_opi_pb_catgtr_for_anc_c1_sample_104_Aque_Scoa.ancstatepostprob
WARNING site 244 maxP 0.29 is tied among 2 states
...
WARNING site 21204 maxP 0.28 is tied among 2 states
# counted 23677 lines
# found probabilities for 23676 sites, target was 23676
# 52 sites multiple states tied for maxP
```

Most of the probabilities are quite high, there is little ambiguity in the state calls. In this dataset for this node, 52 out of 23676 have an ambiguous state, so `ancestral_probability_to_fasta.py` will give them `X` in the fasta.

However, amino acid states are generated for ALL positions, including positions that are mostly missing data. If data are missing, then it will copy the state of the closest neighbor. If many are missing, it will still use the state of the next closest. For instance, if you have a clade where only 1 of 10 species has the gene (say 1 genome vs 9 txomes), then all nodes in that clade will copy the sequence of that lone species. The more missing data, the more that the implied state is biased by fewer species.

## understanding phylobayes output
Since this does not appear to be in the [manual](https://github.com/bayesiancook/pbmpi), this is what the chain appears to contain as best as I can tell:

`head -n 20 simion2017_bottom_10k.chain`

```
(((((Salpingoeca_punica_03:0.12065181307487,Ministeria_vibrans:0.0352483505492743)...;
1
10
1
5000
1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	

0.752456	1.10324	1.6655	0.163412	0.354872	0.248254	0.657502	0.630945	0.728152	0.17976	0.50081	0.69962	0.598637	0.0861455	0.0876559	0.892558	0.0311917	1.3371	0.386137	0.29771	1.39336	0.328663	0.615405	1.2522	1.13879	1.445	0.842922	0.906381	0.754849	0.318021	0.227907	1.5612	0.214874	0.0225296	0.472869	1.00689	0.139297	0.883171	1.67927	1.70409	0.664364	2.06539	0.559632	0.695682	0.801497	3.25483	0.0707738	1.5771	0.466999	1.14675	0.0776192	0.704972	0.0253761	2.69223	0.164571	3.93672	4.63184	0.48477	0.56836	2.92229	0.577695	2.01158	0.771927	0.0566555	0.248625	1.0663	2.96381	1.42662	1.09149	1.90944	2.19236	0.355058	0.703829	0.431295	1.14303	0.483646	3.71787	0.683371	2.13353	1.41832	0.543885	0.754606	1.31183	0.784504	2.52064	0.564413	0.494879	0.338373	1.54069	0.803089	3.09483	0.174367	0.718524	1.21767	0.129853	0.195495	0.0314307	0.280545	0.694194	2.21928	0.431252	0.194706	0.9947	0.0866573	1.37641	0.476275	0.434498	2.62462	0.804932	0.636178	0.306858	0.243217	1.25363	0.439173	0.998232	0.0338243	3.21103	1.41567	0.489619	0.596672	0.182142	0.432279	0.273841	1.25503	1.50839	1.88983	0.567022	3.82803	0.852682	0.359292	0.809722	0.130428	0.157673	0.775841	1.73623	5.32759	1.55121	0.580744	0.87162	0.229947	0.785182	0.617069	0.883117	0.777978	1.01713	1.23474	0.13288	0.178909	0.317019	1.03443	0.331836	0.072567	0.594149	1.62566	1.11257	0.571428	0.342329	0.258291	0.76213	0.546536	0.474314	0.518533	0.705117	1.0516	0.0813968	1.30125	0.440674	0.198017	2.60856	0.0214799	0.70356	1.13102	0.666204	1.11166	2.0335	2.07065	1.10783	0.862025	1.09893	0.837753	1.46409	0.102767	0.583503	0.496023	0.286137	1.11899	0.757297	1.38018	0.148231	2.41132	

0.05092	0.308805	0.0273121	0.21894	0.00244624	0.00380969	0.0402974	0.00654627	0.0333059	0.00462414	0.0543895	0.0144712	0.0187349	0.000976424	0.0364146	0.0354139	0.100391	0.012741	0.0209127	0.00854789	
0.0282293	0.00439629	0.0105243	0.102116	5.72272e-05	0.0559579	0.0761836	0.187137	0.0280789	0.055285	0.0139803	0.0778611	0.0335942	0.00918408	0.0168477	0.0469259	0.00559582	0.0693726	0.103094	0.0755784	
0.00982036	0.0324977	0.038993	0.124121	0.0128107	0.0232892	0.155201	0.0223525	0.193735	0.0273551	0.0118491	0.10893	4.91702e-05	0.0313182	0.0681411	0.0128091	0.0430973	0.0203015	0.0207741	0.0425543	
...
```

This likely corresponds to:

```
1 tree at the start (?) of that iteration
2 ?, value 1 appears constant between iterations
3 ? could be NSPR kappa or kmax, not in any other file
4 appears to be alpha parameter, found in trace file
5 probably refnmodemax, value 5000 appears constant between iterations
6 ?
7 20 numbers, possibly some starting values for each amino acid (last tab is blank)
8 blank
9 190 floats, probably 19*20/2 substitution frequencies (last tab is blank)
10 blank
begin 5000x 20-float lines, probably standing probabilities of each category
5011 in this dataset 10384 tab separated integers, for 10383 sites (last tab is blank)
```

Line 2 is some constant, maybe run status:

`grep -P ";$1" -c simion2017_bottom_10k.chain`

Line 5, the `5000` only appears once in the code, hence is probably that:

`grep 5000 ~/pbmpi/sources/*`

`SBDPProfileProcess.h:const int refnmodemax = 5000;`



