#!/usr/bin/env python
# filter_rfd.py by WRF created 2018-10-22

'''filter pairwise Robinson-Foulds distances to only those in order,
  i.e. 0-1, 1-2, 2-3, etc.

    run with single file in the command line, prints to standard output

filter_rfd.py RAxML_RF-Distances.chain_1_treelist > chain_1_rf_distances

    assumes RF distances are in format from RAxML, using the option -f r
    each line would appear as:
0 1: 96 0.657534
'''

import sys

if len(sys.argv) < 2:
	print >> sys.stderr, __doc__
else:
	for line in open(sys.argv[1],'r'):
		# 0 1: 2 0.014925
		vals = map(int, line.split(":")[0].split(" "))
		if (vals[0]+1) == vals[1]:
			print >> sys.stdout, line.strip()
