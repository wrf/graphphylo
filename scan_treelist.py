#!/usr/bin/env python
#
# scan_treelist.py created by WRF 2018-03-07

'''scan_treelist.py  last modified 2018-03-08

scan_treelist.py -t CAT_GTR1.treelist -m 500 -r Homo_sapiens

    treelist file can be gzipped

  REQUIRES Bio and matplotlib

    then convert to .gif using ImageMagick:

for FILE in *.png; do convert $FILE -crop 640x440+80+20 +repage -resize 70% $FILE.gif ; done

convert -delay 8 -loop 0 *.gif animated_tree.gif

gifsicle -O3 animated_tree.gif --colors 16 -o animated_tree_reduced_color.gif
'''

import sys
import argparse
import os
import time
import gzip
import cStringIO
import matplotlib
from collections import defaultdict
from Bio import Phylo

def get_names_to_colors(colorfile):
	'''read color index file and return a dict where key is species name and value is color'''
	colordict = {}
	print >> sys.stderr, "# reading color vector {}".format(colorfile), time.asctime()
	for line in open(colorfile,'r'):
		line = line.strip()
		if line:
			lsplits = line.split("\t")
			# should be as Genus_species   #aa2255
			colordict[lsplits[0]] = lsplits[1]
	print >> sys.stderr, "# found colors for {} species".format( len(colordict) ), time.asctime()
	return colordict

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-c','--color', help="optional name-to-color file")
	parser.add_argument('-d','--directory', default="./", help="optional directory for temporary figures")
	parser.add_argument('-m','--maximum', type=int, help="stop after this many trees")
	parser.add_argument('-r','--root-taxon', help="species name to root each iteration")
	parser.add_argument('-t','--treelist', help="treelist output file from phylobayes", required=True)
	parser.add_argument('--title', help="title for output graphs")
	args = parser.parse_args(argv)

	colordict = get_names_to_colors(args.color) if args.color else {}

	if not os.path.exists(args.directory):
		os.mkdir(args.directory)
		print >> sys.stderr, "# Creating directory {}".format(args.directory), time.asctime()
	elif os.path.isdir(args.directory):
		print >> sys.stderr, "# Using directory {}".format(args.directory), time.asctime()

	if args.treelist.rsplit('.',1)[1]=="gz": # autodetect gzip format
		opentype = gzip.open
		print >> sys.stderr, "# reading treelist {} as gzipped".format(args.treelist), time.asctime()
	else: # otherwise assume normal open
		opentype = open
		print >> sys.stderr, "# reading treelist {}".format(args.treelist), time.asctime()

	zerolength = 6
	if args.maximum is not None:
		print >> sys.stderr, "# stopping after {} trees".format(args.maximum), time.asctime()
		zerolength = len(str(args.maximum))
	zerostring = "{}.tree-{:0" + str(zerolength) + "}.png"

	linecounter = 0
	treecounter = 0
	eol_counter = defaultdict(int)
	for line in opentype(args.treelist,'r'):
		linecounter += 1
		eol_counter[line[-1]] += 1
		line = line.strip()
		if line:
			treecounter += 1
			treestring = cStringIO.StringIO(line)
			tree = Phylo.read(treestring,"newick")
			tree.root_with_outgroup(args.root_taxon)
			tree.ladderize()
			outtreename = os.path.join(args.directory, zerostring.format(args.treelist,treecounter))
			Phylo.draw(tree, label_colors=colordict, do_show=False)
			matplotlib.pyplot.axis("off")
			matplotlib.pyplot.text(0,len(tree.get_terminals()),"Iter {}".format(treecounter), fontsize=12)
			if args.title is not None:
				matplotlib.pyplot.title(args.title)
			matplotlib.pyplot.savefig(outtreename)
			matplotlib.pyplot.close()
			print >> sys.stdout, "tree", treecounter, "as", outtreename
			if (args.maximum is not None) and (treecounter >= args.maximum):
				print >> sys.stderr, "# encountered final tree {}".format(treecounter), time.asctime()
				break
		else:
			print >> sys.stderr, "WARNING: line {} is empty".format(linecounter)
	print >> sys.stderr, "# {} lines for {} trees".format( linecounter, treecounter )
	#print >> sys.stderr, eol_counter

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
