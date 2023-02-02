#! /usr/bin/env python3

import re
import sys

sys.stderr.write("""This is correct-block coordinates. Takes as input the "blocks/synt" file ouputted by the microsynteny pipeline,
comma-separated list of chrom files, and outputs a .synt file with the correct block coordinates.
Usage: correct_blocks_coordinates.py myblocks.synt mychrom1.chrom,mychrom2.chrom...,mychromn.chrom > corrected_blocks.synt
""")

f_blocks = sys.argv[1]
f_chrom = sys.argv[2].split(',')

d_start = {}
d_stop = {}

for file in f_chrom:
	with open(file, 'r') as f:
		for line in f:
			line = line.strip()
			line = line.split("\t")
			d_start[line[1]] = int(line[4])
			d_stop[line[1]] = int(line[5])

search_string_scaffold = '([^:]+):'


with open(f_blocks, 'r') as g:
	for line in g:
		ls_coordinates = []
		line = line.rstrip()
		line = line.split("\t")
		ls_genes = line[9].split(',')
		scaffold = re.search(search_string_scaffold,line[7]).group(1)
		for gene in ls_genes:
			ls_coordinates.append(d_start[gene])
			ls_coordinates.append(d_stop[gene])
		max_coord = max(ls_coordinates)
		min_coord = min(ls_coordinates)
		output = '\t'.join(line[0:7]) + '\t' +\
					scaffold + ':' + str(min_coord) + '..' + str(max_coord) + '\t' +\
					str(abs(max_coord - min_coord)) + '\t' +\
					line[9]
		print(output)

