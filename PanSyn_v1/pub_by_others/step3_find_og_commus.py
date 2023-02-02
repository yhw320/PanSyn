#!/usr/bin/env python3

import argparse
import csv
import textwrap
import networkx as nx
import synphoni.graph_analysis as sg
from synphoni.logo import logo_ASCII


parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter,
                                 description = textwrap.dedent(f"""\
             
    {logo_ASCII()}
    Step 3 of the SYNPHONI (detection of ancestral SYNteny based on PHylogeny and Ortholog Network Inference) pipeline: 
    Infer ancestral microsyntenic orthogroup sets from network built in step 2.
    Only edges that were microsyntenic in the node will be kept (default: nmax = 30), 
    use tools/step2.5_optimal_nmax.py to determine optimal nmax.
    Only connected components that consist of cliques are kept.
    """))
parser.add_argument("node_dist",
                    help = "Name of the node-specific dist file generated with step 2")
parser.add_argument("-n", "--nmax", 
                    help = "distance threshold above which syntenic links should not be discarded",
                    default = 30,
                    type = float)
parser.add_argument("-o",
                    "--output", 
                    help = "prefix of the output file. gpickle extension is a graph object saved for step 4\
                        csv extension contains the orthogroup communities",
                    type = str,
                    required = True)
args = parser.parse_args()

edgelist, edgelist_filt = sg.load_edgelist(args.node_dist, args.nmax)

G = sg.weighted_edgelist_to_graph(edgelist)
G_filt = sg.weighted_edgelist_to_graph(edgelist_filt)

nx.write_gpickle(G_filt, f"{args.output}.gpickle")

with open(f"{args.output}.csv", 'w') as f:
    outputcsv = csv.writer(f)
    cliques_refined = sg.refine_component_to_cliques(raw_graph = G, filtered_graph = G_filt)
    for refined_og_sets in cliques_refined:
        outputcsv.writerow(refined_og_sets)
