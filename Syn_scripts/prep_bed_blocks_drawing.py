#!/usr/bin/env python3

import argparse
import collections
import pandas as pd
import sys
import pathlib
import csv

parser = argparse.ArgumentParser(
    description = "prepares a bed file for drawing blocks with daw_blocks.rR bed file will be called {block id}_color.bed ")
parser.add_argument("-c", "--multi_species",
                    help = "multi-species blocks file",
                    required = True)
parser.add_argument("-g", "--chrom_folder",
                    help = "folder where the chrom files are located",
                    required = True)
parser.add_argument("-s", "--synt_file", 
                    help = "name of the synt block file",
                    required = True)
parser.add_argument("multi_species_id",
                    help = "id of the multi-species block of interest",
                    type = str)
parser.add_argument("-og", "--ortho",
                    help = "name of the orthology file",
                    type = str,
                    required = True)
args = parser.parse_args()


OKABE_ITO_HEX = ("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
SASHA_TRUBETSKOY_HEX = ('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075',"#000000")
GAINSBORO = "#DCDCDC"

"""
tmp_args = collections.namedtuple("tmp_args", "multi_species chrom_folder synt_file multi_species_id ortho")
args = tmp_args("Planu_m3.len3.ol0.5.phylochecked.clusters",
                "/scratch/robert/2020_08_master_orthology/03_microsynteny/chrom/",
                "Planu_m3.len3.ol0.5.phylochecked.synt",
                "65",
                "../../../03_microsynteny/Orthofinder_n0_hogs_updated_w_singletons.clus")
"""

chromfile_folder = pathlib.Path(args.chrom_folder)
chromfile_list = chromfile_folder.glob("*.chrom")


def block_ids(filepath:str, multi_species_block_id:str) -> list:
    """ load the block ids of a specified multi-species
    raises ValueError if block id not found
    :param filepath: *clusters file, SYNPHONI or MicroSynteny output
    :type filepath: str
    :param block_id: multi species id
    :type block_id: str
    :return: blocks_list
    :rtype: tuple
    """ 
    with open(filepath, 'r') as f:
        for line in f:
            multi_sp, *block_ls = line.rstrip().split()
            if multi_sp == multi_species_block_id:
                return block_ls
    raise ValueError(f"the block with the id {multi_species_block_id} cannot be found in {filepath}")

def load_blocks(filepath:str, block_ls:list) -> list:
    """ loads blocks of interest into dataframe

    :param filepath: *synt file SYNPHONI or MicroSynteny output
    :type filepath: str
    :param block_ls: list of blocks
    :type block_ls: list
    :return: nested list with 3 items block ids, chromosome, species, list of accessions
    :rtype: list
    """
    output_list = []
    with open(filepath, 'r') as f:
        for line in f:
            block_id, species, _, _, _, _, _, coords, _, *acc = line.rstrip().split()
            chromosome = coords.split(":")[0]
            acc_ls = acc[0].split(",")
            if block_id in block_ls:
                output_list.append([block_id, chromosome, species, acc_ls])
    return output_list

def load_chrom(chrom_list:collections.abc.Iterable[str], species_list:list) -> dict[str, pd.DataFrame]:
    """load a list of chrom files into a dict
    keeps only chrom in provided list

    :param chrom_list: list of chrom files
    :type chrom_list: list
    :param species_list: list of species where block is found
    :type species_list: list
    :return: mapping species to loaded chrom file
    :rtype: dict[str, pd.DataFrame]
    """    
    output_dict = {}
    for chromfile in chrom_list:
        mychrom = pd.read_csv(chromfile,
                              names = ["species",
                                       "acc",
                                       "chromosome",
                                       "strand",
                                       "start",
                                       "end"],
                              sep = "\t")
        species = mychrom["species"][0]
        if species in species_list:
            output_dict[species] = mychrom
    return output_dict


def load_ortho(filepath:str) -> dict[str, list]:
    """loads orthology file

    :param filepath: tsv file, first field is OG name, second nc of accessions in OG
    third field and onward are accessions
    :type filepath: str
    :return: accessions to orthogroup mapping
    :rtype: dict[str, list]
    """
    output_dict = {}
    with open(filepath, 'r') as f:
        for line in f:
            og,_, *acc_ls = line.rstrip().split()
            tmp_dict = {acc: og for acc in acc_ls}        
            output_dict.update(tmp_dict)
    return output_dict


def og_palette(og_list: set[str]) -> dict[str, str]:
    """build an orthgroup palette depending on the number of dsitinct OGs the block comprises

    :param og_list: set of orthogroup names
    :type og_list: set[str]
    :return: mapping of OG names to color palette
    :rtype: dict[str, str]
    """
    colormap = {"intervening": GAINSBORO}
    og_list = og_list - {"intervening"}
    if len(og_list) <= len(OKABE_ITO_HEX):
        palette = OKABE_ITO_HEX
    elif len(og_list) <= len(SASHA_TRUBETSKOY_HEX):
        palette = SASHA_TRUBETSKOY_HEX
    elif len(og_list) > len(SASHA_TRUBETSKOY_HEX):
        palette = ["#000000"] * len(og_list)
    return colormap | {og: color for og, color in zip(og_list, palette)}


def get_intervening_coords(indexlist: list, chromosome_data: pd.DataFrame) -> list:
    """get start/end coords of intervening genes

    :param indexlist: indices of genes part of the block
    :type indexlist: list
    :param chromosome_data: loaded chromfile of chromosome of interest
    :type chromosome_data: pd.DataFrame
    :return: list  of coords + number of intervening genes
    :rtype: list
    """
    output = []
    for i in range(len(indexlist) - 1):
        idx1, idx2 = indexlist[i: i + 2]
        nb_intervening = idx2 - idx1 - 1
        if nb_intervening > 0:
            start_i = chromosome_data.loc[idx1]["start"] + 1
            end_i = chromosome_data.loc[idx2]["start"] - 1
            output.append((start_i, end_i, nb_intervening))
    return output

def summarize_block(sp_block:list, gene_location:dict[str, pd.DataFrame], orthology: dict[str, str], spacer = 0.2) -> list:
    """summarize blocks and intervening genes, color by OG

    :param sp_block: nested list with block_id, chromosome, species, acc_ls of a block
    :type sp_blocks: list
    :param gene_location: summary of gene locations
    :type gene_location: dict[str, pd.DataFrame]
    :param spacer: how big spacer is relative to gene
    :return: nested list with 
    molecule gene start end and color (HEX code)
    :rtype: pd.DataFrame
    """
    output = []
    block_id, chromosome, species, acc_ls = sp_block
    molecule = f"{species}_{block_id}"
    df_chr = gene_location[species].query("chromosome == @chromosome").sort_values(by = "start")
    df_chr.index = range(0, len(df_chr.index))
    df_block = df_chr.query("acc in @acc_ls")
    df_block["og"] = df_block['acc'].map(orthology)
    colormap  =  og_palette(set(df_block["og"]))
    idx = df_block.index
    intervening_genes = get_intervening_coords(idx, df_chr)
    for _, gene, _, _, start, end, og in df_block.values.tolist():
        output.append([molecule, gene, start, end, og])
    for start, end, gene in intervening_genes:
        output.append([molecule, gene, start, start+1, "intervening"])
    output = sorted(output, key = lambda x: x[2])
    gene_size, spacer, gene_start = 4, 1, 1
    final_output = []
    for mol, acc, _,_, og in output:
        if og == "intervening":
            gene_end = gene_start + spacer
        else:
            gene_end = gene_start + gene_size
        final_output.append([mol, acc, gene_start, gene_end, og])
        gene_start = gene_end + spacer
    return final_output


my_block_ids = block_ids(args.multi_species,
                         args.multi_species_id)

my_blocks = load_blocks(args.synt_file, my_block_ids)

species_ls = [sp for *_, sp,_ in my_blocks]

gene_loc = load_chrom(chromfile_list, species_ls)

ortho = load_ortho(args.ortho)


output_table = []
for block in my_blocks:
    output_table.extend(summarize_block(block, gene_loc, ortho))

all_ogs = set([og for *_, og in output_table])
colormap = og_palette(all_ogs)

output_table = pd.DataFrame(output_table, columns = ["molecule", "gene", "start", "end", "og"])
output_table["color"] = output_table["og"].map(colormap)



OUTFILE = f"{args.multi_species_id}_color.tsv"
output_table.to_csv(OUTFILE, index = False, sep = "\t")