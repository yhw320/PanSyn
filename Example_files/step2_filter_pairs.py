#!/usr/bin/env python3

import ete3
import pathlib
import csv
import synphoni.utils as su
import textwrap
from synphoni.logo import logo_ASCII
import argparse


parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter,
                                 description = textwrap.dedent(f"""\
            
    {logo_ASCII()}
    Step 2 of the SYNPHONI (detection of ancestral SYNteny based on PHylogeny and Ortholog Network Inference) pipeline: 
    Determine ancestral distance in node(s) of interest n specified in a species tree t
    The orthogroup pairs are filtered using the following criterias:
        * Either 2+ children, or 1 children, 1+ sister group our 1 children and outgroup needs to be populated.
        * Group is said populated if at least m species (species/threshold) are found in said group (default m is 2)
    Ancestral distance is inferred for every node
    One edgelist generated by node specified.
    """))
parser.add_argument("-s", "--species_tree",
                    help = "Tree of the PREFIX names extracted from chrom files in step 1, newick format with node names\
                            use a tree where only the names of species in your dist file are present\
                            (e.g. (((D:0.723274,F:0.567784)E:0.067192,(B:0.279326,H:0.756049)B:0.807788); ).",
                    required = True)
parser.add_argument("-n", "--node_names", 
                    help = "space-delimited list of nodes which will be used to filter the OG pairs\
                        every node is run sequentially, so launch multiple commands with only one node in order to parallelize",
                    required = True,
                    nargs = "+")
parser.add_argument("-m", "--species_threshold",
                    help = "How many species of a group should posess the syntenic pair to consider group as populated",
                    type = int,
                    default = 3)
parser.add_argument("-d", "--dist",
                    help = "path to the dist file created by step 1",
                    type = str,
                    required = True)
parser.add_argument("-o", "--output_path",
                    help = "path to the outputfiles",
                    type = str,
                    required = True)
	
args = parser.parse_args()
speciestree = ete3.Tree(args.species_tree, format=1)
su.prune_tree(speciestree, args.dist, args.species_tree)
ls_node_names = [su.clade(species, speciestree) for species in args.node_names]
inputpath = pathlib.Path(args.dist)
output_prefix, output_extension = inputpath.stem, inputpath.suffix
for node in ls_node_names:
    output_filename = args.output_path + f"{output_prefix}.m_{args.species_threshold}.{node.name}{output_extension}"
    with open(args.dist, "r") as f, open(output_filename, "w") as g:
        inputcsv = csv.reader(f)
        outputcsv = csv.writer(g)
        outputcsv.writerow(["OGx", "OGy", node.name])
        _,_, *species_ls = next(inputcsv)
        for row in inputcsv:
            OGx, OGy, *dist_ls = row
            all_dist = {sp:int(dist) for sp, dist in zip(species_ls, dist_ls) if dist != ""}
            dist = su.get_ancestral_distance(all_dist, node, args.species_threshold)
            if dist != None:
                outputcsv.writerow([OGx, OGy, dist])
