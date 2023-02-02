#!/usr/bin/env python3

import argparse
import sys
import ete3

parser = argparse.ArgumentParser(description = 'Provides reports of block content of specified nodes')
parser.add_argument('-c',
                    '--clusters_id',
                    help = 'Tsv file, the multi-species clusters of microsyntenic blocks \
                        output of microsynteny pipeline (makeClusters3.pl).',
                    required = True)
parser.add_argument('-b',
                    '--block_list',
                    help = 'Tsv file, microsyntenic blocks details, output of microsynteny pipeline\
                    (makeClusters3.pl).',
                    required = True)
parser.add_argument('-s',
                    '--species_tree',
                    help = 'Tree of the PREFIX in the second column of the block list,\
                        newick format with node names (e.g.\
                        (((D:0.723274,F:0.567784)E:0.067192,(B:0.279326,H:0.756049)B:0.807788); ).',
                    required = True)
parser.add_argument('-n',
                    '--node_names',
                    help = 'Space separated list of names of nodes of interest from in the species tree\
                    (e.g. "E B A C".',
                    nargs = '+',
                    required = True)
parser.add_argument('-m',
                    '--species_threshold',
                    help = "Minimum number n of species per clade for node inference\
                    (e.g. for novel blocks, requires n species in at least two ingroups,\
                    and no species in outgroup\
                    for ancestral blocks, n species in one ingroup an n in one outgroup).\
                    By default, n=2. For ingroup/outgroup of size < n, all species of said\
                    ingroup/outgroup are required to posess the block\
                    (e.g. if ingroup size of 1, the species is required to posess.",
                    type = int,
                    default = 2)
parser.add_argument('-r',
                    '--report',
                    help = 'Type of report printed to STDOUT [Default: short] : \n\
                            - "short": number of blocks per node. \n\
                            - "clusters_list": filters one multi-species block per line,\
                                cluster IDs (field 1), ancestral/novel nodes of the block (field 2),\
                                species list (field 3), block_ids (field 4+).\
                                For getting a *.clusters file with only a subset,\
                                pipe SyntByNode output to `cut -f1,4-`.\n\
                            - "blocks_list": blocks within the filtered multi-species blocks\n\
                            - "tree_ASCII": block count per specified node on an ASCII tree\n\
                            - "tree_NH": block count per specified node on a newick tree',
                    choices = ['short', 'clusters_list', 'blocks_list', 'tree_ASCII', 'tree_NH'],
                    default = 'short')
parser.add_argument('-t',
                    '--block_type',
                    help = 'Specify whether you want to report all blocks ("total"),\
                        only the ones inherited from older nodes ("ancestral") or only "novel" ones.\
                        Multiple options can be specified [Default: total]',
                    choices = ['total','ancestral','novel'],
                    nargs = '+',
                    default = 'total')
args = parser.parse_args()

#Exit the script if you ask for list with more than one type of blocks to report
if len(args.block_type) > 1 and 'list' in args.report:
    sys.stderr.write(f'the {args.report} option can only be used when searching\
        for only one type of blocks! (ONLY total, ONLY ancestral or ONLY novel.)\n')
    sys.exit()



#Makes a generator of lists, each list being an ingroup.
def get_ingroups(taxonomic_node, speciestree):
    """
    :param taxonomic node: of whcih the ingroups and outgroups must be determined
    :param species tree: tree to use for determining outgroups and ingroups
    :return: nested list of ingroups, one list per children node
    """
    output_list = []
    for node in speciestree.traverse("preorder"):
        if node.name == taxonomic_node:
            for ingroup in node.children:
                tmp_ig_list = []
                for TaxonName in ingroup.iter_leaf_names(is_leaf_fn=None):
                    tmp_ig_list.append(TaxonName)
                output_list.append(tmp_ig_list)
    return output_list

def get_outgroups(taxonomic_node, speciestree):
    """
    :param taxonomic node: of whcih the ingroups and outgroups must be determined
    :param species tree: tree to use for determining outgroups and ingroups
    :return: list of outgroups
    """
    output_list = []
    cached_tree = speciestree.copy()
    for node in cached_tree.traverse("postorder"):
        if node.name == taxonomic_node and node.is_root() is False:
            node.detach()
        elif node.name == taxonomic_node and node.is_root() is True:
            return output_list
    for taxon_name in cached_tree.iter_leaf_names(is_leaf_fn = None):
        output_list.append(taxon_name)
    return output_list


def blockspeciespairs(blockinfofile):
    """
    :param blockinfofile: file with block info
    :returns: a dict with blocks as keys, species as values
    """
    output_dict = {}
    with open(blockinfofile, 'r') as f:
        for line in f:
            block_id, species, *_ = line.rstrip().split('\t')
            output_dict[block_id] = species
    return output_dict


def blockclusterpairs(clusterinfofile):
    """
    :param clusterinfofile: file with multi-species block info
    :returns: a dict with block ids as keys, multisp ids as values
    """
    output_dict = {}
    with open(clusterinfofile, 'r') as f:
        for line in f:
            cluster_id, *block_id_ls = line.rstrip().split('\t')
            for block_id in block_id_ls:
                output_dict[block_id] = cluster_id
    return output_dict


class clade: 
    """
    A clade instance. As the class is created, we automatically isolate ingroups and outgroups using the user-provided speciestree.
    'blocks' attributes are lists of the cluster_ids (i.e. multi-species blocks)
    """
    def __init__(clade, name):
        clade.name = name
        clade.ingroups = get_ingroups(clade.name, speciestree)
        clade.outgroups = get_outgroups(clade.name, speciestree)
        clade.total_blocks = []
        clade.novel_blocks = []
        clade.ancestral_blocks = []


def get_block_type(clade, specieslist, n):
    """
    Determines whether a block is ancestral or novel at a node of interest
    :param clade:  clade instance, node of interest
    :param specieslist: list of species posessing the block
    :param n: species threshold. If block novel, needs to be found in 2 ingroups in at least n species
        if ancestral, block needs to be found in n species in ingroup, n species in outgroup
    :returns: a list of the block type states, if it's novel or ancestral it's also part of total
    """
    nb_species_per_ingroup = []
    for ingroup in clade.ingroups:
        species_IG = [x in ingroup for x in specieslist] #bool array, sum of it is nb of positives since True counts as one.
        if len(ingroup) >= n:
            nb_species = sum(x in ingroup for x in specieslist)
        elif len(ingroup) < n: #in the event that an ingroup is smaller than the specified threshold
            if sum(x in ingroup for x in specieslist) == len(ingroup): 
                nb_species = n # if all the species of the ingroup posess the pair, threshold is satisfied
            else:
                nb_species = sum(x in ingroup for x in specieslist)
        nb_species_per_ingroup.append(nb_species)
    nb_populated_ingroups = sum(x >= n for x in nb_species_per_ingroup)
    nb_species_outgroups = sum(x in clade.outgroups for x in specieslist)
    nb_species_ingroups = sum(nb_species_per_ingroup)
    pair_is_novel = [nb_populated_ingroups >= 2 and nb_species_outgroups == 0]
    pair_is_ancestral = [nb_populated_ingroups >= 2 and nb_species_outgroups > 0,
                         nb_species_ingroups >= n and nb_species_outgroups >= n,
                         nb_populated_ingroups > 0 and nb_species_outgroups == len(clade.outgroups)]
    output_ls = []
    if any(pair_is_ancestral):
        return ['ancestral', 'total']
    elif any(pair_is_novel):
        return ['novel', 'total']
    else:
        return []


def get_node_states(clustersIDfile, list_clade_instances, dictspeciesblock):
    with open(clustersIDfile,"r") as f:
        for line in f:
            cluster_id, *block_id_ls = line.rstrip().split('\t')
            species_list = []
            for block_id in block_id_ls:
                species_list.append(dictspeciesblock[block_id])
                set_species_list = set(species_list)
            for i in range(len(list_clade_instances)):
                current_block_type = get_block_type(list_clade_instances[i], set_species_list, args.species_threshold)
                if 'total' in args.block_type and 'total' in current_block_type:
                     list_clade_instances[i].total_blocks.append(cluster_id)
                if 'novel' in args.block_type and 'novel' in current_block_type:
                     list_clade_instances[i].novel_blocks.append(cluster_id)
                if 'ancestral' in args.block_type and 'ancestral' in current_block_type:
                     list_clade_instances[i].ancestral_blocks.append(cluster_id)


def print_clusters_list(clustersIDfile,list_clade_instances, dictspeciesblock):
    with open(clustersIDfile,"r") as f:
        for line in f:
            cluster_id, *block_id_ls = line.rstrip().split('\t')
            block_id_str = "\t".join(block_id_ls) 
            species_list, node_list = [], []
            for taxon in list_clade_instances:
                keep_block = False
                conditions_to_keep_block = [
                    'total' in args.block_type and cluster_id in taxon.total_blocks,
                    'ancestral' in args.block_type and cluster_id in taxon.ancestral_blocks,
                    'novel' in args.block_type and cluster_id in taxon.novel_blocks]
                if any(conditions_to_keep_block): keep_block = True
                else: keep_block = False
                if keep_block is True:
                    for block_id in block_id_ls:
                        species_list.append(dictspeciesblock[block_id])
                        set_species_list = set(species_list)
                        node_list.append(taxon.name)
            if node_list != []:
                node_ls_str = ','.join(set(node_list))
                print(f'{cluster_id}\t{",".join(set(node_list))}\t{",".join(set_species_list)}\t{block_id_str}')


def print_blocks_list(blocksIDfile,list_clade_instances, dictblockclusters):
    with open(blocksIDfile, 'r') as f:
        for line in f:
            block_id, line_rest = line.rstrip().split('\t', maxsplit = 1)
            node_list = []
            try:
                cluster_id = dictblockclusters[block_id]
                for taxon in list_clade_instances:
                    if 'total'in args.block_type and cluster_id in taxon.total_blocks:
                        node_list.append(taxon.name)
                    if 'ancestral' in args.block_type and cluster_id in taxon.ancestral_blocks:
                        node_list.append(taxon.name)
                    if 'novel' in args.block_type and cluster_id in taxon.novel_blocks:
                        node_list.append(taxon.name)
                if node_list != []:
                    print(f'{cluster_id}\t{block_id}\t{line_rest}')
            except KeyError:
                pass


def print_short(list_clade_instances):
    ls_taxons = []
    for taxon in list_clade_instances:
        ls_taxons.append(taxon.name)
    print('\t'.join(['taxon', 'blocktype', 'count']))
    for taxon in list_clade_instances:
        if 'total' in args.block_type:
            print(f'{taxon.name}\ttotal\t{len(taxon.total_blocks)}')
        if 'ancestral' in args.block_type:
            print(f'{taxon.name}\tancestral\t{len(taxon.ancestral_blocks)}')
        if 'novel' in args.block_type:
            print(f'{taxon.name}\tnovel\t{len(taxon.novel_blocks)}')


def print_tree (list_clade_instances, species_tree):
    cached_tree = speciestree.copy() 
    ls_taxons = []
    for taxon in list_clade_instances:
        ls_taxons.append(taxon.name)
    for node in cached_tree.traverse("preorder"):
        if node.name in ls_taxons:
            i = getattr(ls_taxons, 'index')(node.name) #ls_taxons and ls_clade_instances are essentially the same list order. We get the indes from ls_taxons, and use this index to extract class info from ls_clade_instances 
            node_clade = list_clade_instances[i]
            if 'total' in args.block_type:
                node.name = f'{node.name}_{len(node_clade.total_blocks)}'
            if 'ancestral' in args.block_type:
                node.name = f'{node.name}_{len(node_clade.ancestral_blocks)}'
            if 'novel' in args.block_type:
                node.name = f'{node.name}_{len(node_clade.novel_blocks)}'
    if args.report == 'tree_ASCII':
        print(cached_tree.get_ascii(show_internal=True))
    if args.report == 'tree_NH':
        print(cached_tree.write(format = 1))


#import the speciestree
speciestree = ete3.Tree(args.species_tree, format=1)

#argument check. if one of the node names is in root and asks for novel blocks, stop the script
if speciestree.get_tree_root().name in args.node_names:
    raise AttributeError(f'It is not possible to determine novel blocks\
        of {speciestree.get_tree_root().name }, as it is the root node')



#create useful dictionaries for the following functions
d_block_species = blockspeciespairs(args.block_list)
d_block_clusters = blockclusterpairs(args.clusters_id)

#modify args.node_names into a list of clade instances.
for x in range(0,len(args.node_names)):
    args.node_names[x] = clade(args.node_names[x])

# This modifes clade instances within the list. Adding lists of cluster to the attributes created as empty lists at the __init__ of the clade instance
get_node_states(args.clusters_id, args.node_names, d_block_species)

if 'tree' in args.report:
    print_tree(args.node_names, speciestree)
    sys.stderr.write(f'Nodes are as follows: {"_".join(args.block_type)}\n')
elif args.report == 'short':
    print_short(args.node_names)
elif 'list' in args.report:
    sys.stderr.write(f'The {", ".join(args.block_type)} blocks are printed to stdout.\n')
    if args.report == 'clusters_list':
        print_clusters_list(args.clusters_id, args.node_names, d_block_species)
    else:
        print_blocks_list(args.block_list, args.node_names, d_block_clusters)