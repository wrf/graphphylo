#!/usr/bin/env python

'''
change_model_order.py -i bur_hel_model.txt > bur_hel_model.hpi.txt

    normal PAML AA order is alphabetical:
  ARNDCQEGHILKMFPSTWYV
    default new order is by hydrophobicity index:
  DPENKRQSGHTACYMVWLIF
'''

import sys
import argparse
from collections import defaultdict

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-i','--input', help="substitution matrix, in PAML format")
	parser.add_argument('-n','--new-order', default="DPENKRQSGHTACYMVWLIF", help="new amino acid order")
	parser.add_argument('--default-order', default="ARNDCQEGHILKMFPSTWYV", help="original amino acid order [ARNDCQEGHILKMFPSTWYV]")
	args = parser.parse_args(argv)

	aavalues = defaultdict(dict) # dict of dicts where key is AA, value is dict of AA to score
	sys.stderr.write("# reading substitution matrix from {}, using AA order as: {}\n".format( args.input , args.default_order ) )
	for i,line in enumerate(open(args.input, 'r')):
		line = line.strip()
		lsplits = line.split(' ') # assume space separated values
		for j,score in enumerate(lsplits):
			if i < 19:
				# first should be A to R / R to A meaning 0-1 / 1-0
				# then A to N, R to N / N to A, N to R meaning 0-2 1-2 / 2-0 2-1
				aavalues[ args.default_order[i+1] ][ args.default_order[j] ] = score
				aavalues[ args.default_order[j] ][ args.default_order[i+1] ] = score
			else: # last line is overall frequency with 20 values
				aavalues[ args.default_order[j] ][ args.default_order[j] ] = score

	sys.stderr.write("# writing new substitution matrix in order of {}\n".format( args.new_order ) )
	for i,let1 in enumerate(args.new_order):
		scorelist = []
		for j,let2 in enumerate(args.new_order[:i]):
			scorelist.append(aavalues[let1][let2])
		if scorelist:
			sys.stdout.write("{}\n".format( " ".join(scorelist) ) )
	basefreq = [aavalues[letter][letter] for letter in args.new_order]
	sys.stdout.write("{}\n".format( " ".join(basefreq) ) )

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
