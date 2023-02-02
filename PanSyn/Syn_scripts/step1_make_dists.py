#! /usr/bin/env python3

import csv
import itertools
import pathlib
import sys
import textwrap
import synphoni.parsers as sp
import synphoni.utils as su
from synphoni.logo import logo_ASCII
import argparse

parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter,
                                 description = textwrap.dedent(f"""\
    
    {logo_ASCII()}
    Step 1 of the SYNPHONI (detection of ancestral SYNteny based on PHylogeny and Ortholog Network Inference) pipeline: 
    Detection of conserved syntenic pairs of genes, inter- and intra- OGs.\n
    Takes chrom files (PREFIX.chrom) and orthology files as input, Creates one folder containing two files
        * one csv with 2+n columns (n number of species), describing the minimal number of intervening genes
        between 2 syntenic OGs (1-based, i.e. adjacent genes have a distance of 1).
        * one pickle file where coordinates and user-provided data is saved for step 4\n
        """))
parser.add_argument("-c",
                    "--clusfile",
                    help = "Orthology file, clus format.",
                    type = str,
                    required = True)
parser.add_argument("chromfiles",
                    help = "List of chrom files, space separated.",
                    nargs = "+")
parser.add_argument("-o",
                    "--output",
                    help = "Name of the output folder\
                        (default:synphoni_results)",
                    type = str,
                    default = "synphoni_results")
parser.add_argument("-m",
                    "--maxpara",
                    help = "Max number of paralogs allowed for any species in a given OG.\
                        OGs with one species above this threshold will be discarded\
                        (default: 100)",
                    type = int,
                    default = 100)
parser.add_argument("-s",
                    "--species_threshold",
                    help = "minimum number of species where orthoghroup pair needs to be syntenic for being retained\
                        since this step is the longest, this parameter should be set at 2 and stringency \
                            on block detection should rely on parameters ar other steps\
                        (default: 2)",
                    type = int,
                    default = 2)
args = parser.parse_args()
su.make_output_folder(args.output)
output_abspath = pathlib.Path(args.output).absolute()
output_dist = f"{str(output_abspath)}/{output_abspath.name}.dist"
chrom_data_pic = f"{args.output}/chromdata.pickle"
min_sp = su.get_threshold(species_threshold = args.species_threshold,
                          chromfiles_list = args.chromfiles)
ortho = sp.read_orthology(args.clusfile)
chrom, species_ls = sp.read_chromfiles(chromfiles_list = args.chromfiles,
                                       orthology_dict = ortho,
                                       min_species = args.species_threshold)
su.filter_OGs(chrom_data = chrom,
              max_paralogs = args.maxpara)
og_ls = list(chrom.keys())
nb_OG = len(og_ls)
sys.stderr.write(f"{nb_OG} OGs retained\n")
su.save_chrom_data(chrom_data = chrom,
                   filepath = chrom_data_pic)
with open(output_dist, "w") as csv_dist_f:
    dist_writer = csv.writer(csv_dist_f)
    header = ["OG1", "OG2"] + species_ls
    dist_writer.writerow(header)
    sys.stderr.write(f"Looking for syntenic OG pairs...\n")
    og_pairs = itertools.combinations_with_replacement(og_ls, 2)
    for pair in og_pairs:
        row = su.get_min_distances_pair(pair, chrom, min_sp, species_ls)
        if row is not None:
            dist_writer.writerow(row)
