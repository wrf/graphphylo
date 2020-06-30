#!/usr/bin/env python
#
# ancestral_probability_to_fasta.py

'''ancestral_probability_to_fasta.py  last modified 2020-06-30

    convert the ancestral state probability from phylobayes to fasta file of amino acids

ancestral_probability_to_fasta.py -p -p sample_87.ancstatepostprob sample_86.ancstatepostprob sample_85.ancstatepostprob sample_95.ancstatepostprob sample_104.ancstatepostprob > ancestral_states.fasta

'''

import sys
import argparse
import gzip
import time
from Bio import Seq

threesites='''
23676	20	A	C	D	E	F	G	H	I	K	L	M	N	P	Q	R	S	T	V	W	Y
1		0.01	0	0	0	0	0	0	0	0	0	0.99	0	0	0	0	0	0	0	0	0
2		1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
3		0	0	0	0	0	0	0.01	0.02	0.01	0.03	0	0	0	0	0	0	0	0.93	0	0
'''

def probabilities_to_aa(probsfile, report_problems=True):
	'''from the file of probabilities, return a string of the most probable AA for each site'''
	if probsfile.rsplit('.',1)[-1]=="gz": # autodetect gzip format
		opentype = gzip.open
		sys.stderr.write("# reading probabilities {} as gzipped  {}\n".format(probsfile, time.asctime() ) )
	else: # otherwise assume normal open
		opentype = open
		sys.stderr.write("# reading probabilities {}  {}\n".format(probsfile, time.asctime() ) )

	linecounter = 0
	ancestral_states = [] # list of strings of AAs
	problem_equal_sites = [] # list of positions that have equal maxP
	problem_low_sites = [] # list of positions that have very low probabilities
	site_target_count = 0 # will be assigned below as sanity check

	for line in opentype(probsfile,'rt'):
		line = line.strip()
		if line and line[0]!="#":
			linecounter += 1
			lsplits = line.split("\t")
			if linecounter == 1: # the first line should be the number of sites
				site_target_count = lsplits[0]
				num_site_types = float(lsplits[1])
				aa_order = lsplits[2:]
			else: # meaning all other sites
				site_num = lsplits[0]
				aa_probabilities = map(float,lsplits[2:])
				if len(aa_probabilities) > num_site_types:
					raise ValueError("ERROR: line has {} elements, should have {}".format( len(aa_probabilities),num_site_types) )
				prob_max = max(aa_probabilities)
				prob_max_index = aa_probabilities.index(prob_max)
				prob_max_aa = aa_order[prob_max_index]

				# track number of sites with low probabilities
				if prob_max <= 1.0/num_site_types: # this should never happen
					problem_low_sites.append(site_num)
					sys.stderr.write("WARNING site {} has maxP {}\n".format(site_num, prob_max) )

				# check for sites with ties for max
				sites_w_max = aa_probabilities.count(prob_max)
				if sites_w_max > 1:
					problem_equal_sites.append(site_num)
					if report_problems:
						sys.stderr.write("WARNING site {} maxP {} is tied among {} states\n".format(site_num, prob_max, sites_w_max) )
					ancestral_states.append("X")
				else:
					ancestral_states.append(prob_max_aa)
	sys.stderr.write("# counted {} lines  {}\n".format(linecounter, time.asctime() ) )
	sys.stderr.write("# found probabilities for {} sites, target was {}\n".format( len(ancestral_states), site_target_count ) )
	if problem_equal_sites:
		sys.stderr.write("# {} sites multiple states tied for maxP\n".format( len(problem_equal_sites) ) )
	if problem_low_sites:
		sys.stderr.write("# {} sites had maxP less than {:.02f}\n".format( len(problem_low_sites), 1.0/num_site_types ) )

	return ancestral_states

def make_fasta_output(filename, ancestral_aas):
	output_header = filename.replace(".ancstatepostprob","")
	sys.stdout.write(">{}\n".format(output_header))
	padded_seq = ""
	for i in xrange(0,len(ancestral_aas),60):
		padded_seq += "".join(ancestral_aas[i:i+60]) + "\n"
	sys.stdout.write("{}".format(padded_seq))

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-p','--probabilities', metavar="FILE", nargs='*', help=".ancstatepostprob file of ancestral state posterior probabilities, can be multiple files, can be .gz")
	parser.add_argument('-q','--quiet', action="store_false", help="suppress warnings of sites with problems, equal likelihoods, etc")
	args = parser.parse_args(argv)

	for probsfile in args.probabilities:
		ancestral_aa_string = probabilities_to_aa(probsfile, args.quiet)
		make_fasta_output(probsfile, ancestral_aa_string)

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
