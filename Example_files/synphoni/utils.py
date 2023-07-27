import numpy as np
import scipy
import csv
import sys
import pathlib
import collections
import itertools
import math
import pickle


class clade:
    """
    A clade object. ingroups, sistergroup(s) and outgroups are isolated using the user-provided speciestree at instantiation.
    """
    def __init__(self, name, speciestree):
        """
        creates a clade object, from a species tree and a node name
        :param name: name of the noce where we define ingroups/outgroups
        :param speciestree: an ete3.Tree object, with nodenames
        :return: clade
        """
        self.name = name
        self.tree = speciestree
        self.ingroups = self.get_ingroups()
        self.outgroups = self.get_outgroups()
        self.sistergroups = self.get_sistergroups()
    def get_ingroups(self):
        """
        :return: nested list of ingroups of the clade
        """
        output_list = []
        for node in self.tree.traverse("preorder"):
            if node.name == self.name:
                for ingroup in node.children:
                    tmp_ig_list = []
                    for TaxonName in ingroup.iter_leaf_names(is_leaf_fn=None):
                        tmp_ig_list.append(TaxonName)
                    output_list.append(tmp_ig_list)
        return output_list
    def get_outgroups(self):
        """
        :return: flat list of outgroup species
        if the node is root, an empty list is returned
        """ 
        cached_tree = self.tree.copy()#create a copy in a different memory space of the species tree. node.detach() will affect cached_tree, but not the initial import of speciestree. For each new outgroup search, we'll bring back the copy of cached tree to its original state.
        for node in cached_tree.traverse("postorder"):
            if node.name == self.name:
                if node.is_root() is False:
                    node.detach()
                else:
                    return []
        return [name for name in cached_tree.iter_leaf_names(is_leaf_fn = None)]
    def get_sistergroups(self):
        """
        :return: nested list of sister clades of the node
        if the node is root, an empty list is returned
        """
        output_list = []
        for node in self.tree.traverse("preorder"):
            if any([node.name == self.name for node in node.children]):
                for sistergroup in [x for x in node.children if x.name != self.name]:
                    tmp_sis_list = []
                    for TaxonName in sistergroup.iter_leaf_names(is_leaf_fn = None):
                        tmp_sis_list.append(TaxonName)
                    output_list.append(tmp_sis_list)
                return output_list
        return output_list


def estimate_maxima(data, no_samples = 1000):
    """
    estimates value where the global maxima of the gaussian kernel density estimate
    :param data: a list
    :param no_samples: number of samples to use
    :return: a float, estimated maxima
    """
    if len(data) <= 2:
        return np.median(data)
    elif data.count(data[0]) == len(data):#all elements are equal
        return data[0]
    else:
        samples = np.linspace(np.min(data), np.max(data), no_samples)
        probs = scipy.stats.gaussian_kde(data).evaluate(samples)
        maxima_index = probs.argmax()
        maxima = samples[maxima_index]
        return maxima


def get_ancestral_distance(dict_distances, myclade, species_threshold):
    """
    calculates ancestral distance between two OG pairs, relies on estimate_maxima
    :param dict_distances: keys = species PREFIX, values = dists (NA included)
    :param myclade: clade instance
    :return: if the pair was present in the node, returns inferred dist
    if pair not found in ingroup, or is a novelty from only one child returns None
    """
    child_dists = []
    sister_dists = []
    outgroup_dists = []
    mean_outside_dist = None
    for child_ingroup in myclade.ingroups:
        child_ingroup_dist = [dict_distances[sp] for sp in child_ingroup
                                 if dict_distances.get(sp) != None]
        if len(child_ingroup_dist) >= species_threshold:
            child_dists.append(child_ingroup_dist)
        elif child_ingroup_dist != []: #if the threshold is larger than the clade, all species required to have the pair
            if len(child_ingroup_dist) == len(child_ingroup):
                child_dists.append(child_ingroup_dist)
    if len(child_dists) == 0:
        return None
    elif len(child_dists) == 1:#calculate dists within sister or outgroup only if needed
        for sistergroup in myclade.sistergroups:
            sistergroup_dist = [dict_distances[sp] for sp in sistergroup
                                  if dict_distances.get(sp) != None]
            if len(sistergroup_dist) >= species_threshold:
                sister_dists.append(sistergroup_dist)
            elif len(sistergroup_dist) == len(sistergroup):
                sister_dists.append(sistergroup_dist)
        if sister_dists == []: #if we don't find info in the sister group, we go look into all the outgroups
            outgroup_dists = [dict_distances[sp] for sp in myclade.outgroups
                               if dict_distances.get(sp) != None]
            if outgroup_dists == []: #check, if we're looking at the root node, outgroup size and outgroup dists would have the same length (empty lists of len 0)
                return None
            elif len(outgroup_dists) < species_threshold:
                if len(outgroup_dists) < len(myclade.outgroups):
                    return None
                elif len(outgroup_dists) > species_threshold:
                    mean_outside_dist =  estimate_maxima(outgroup_dists)
            else:
                mean_outside_dist = estimate_maxima(outgroup_dists)
        else:
            mean_outside_dist = np.mean([estimate_maxima(x) for x in sister_dists])
        if mean_outside_dist is None:
            return None
    dist_children = []
    for x in child_dists:
        dist_children.append(estimate_maxima(x))
    mean_dist_children = np.mean(dist_children)
    if len(child_dists) > 1:
        return mean_dist_children
    else:
        mean_ingroup_outside = np.mean([mean_dist_children, mean_outside_dist])
        return mean_ingroup_outside


def nested_dict():
        """ Create a nested dict. Function is picklable """
        return collections.defaultdict(nested_dict)


def flatten(nested_list):
    """
    flattens a nested list
    :param nested_list: a nested list
    :return: a flattened list
    """
    return list(itertools.chain.from_iterable(nested_list))


def prune_tree(mytree, distfile, input_filename):
    """
    edits inplace the input tree so only leaves in input dist file are kept
    saves to a newick file the pruned tree
    if a taxon from the distfile isn't in the tree, stops the script
    :param mytree: ete3.Tree object
    :param distfile: output of syntkernel step1
    :return: nothing, edits mytree inplace
    """
    leaves_input_tree = mytree.get_leaves()
    with open(distfile, 'r') as f:
        csvreader = csv.DictReader(f)
        _,_, *species_ls = csvreader.fieldnames
    try:
        mytree.prune(species_ls)
        leaves_pruned_tree = mytree.get_leaves()
    except ValueError:
        taxon_not_in_dist = [sp for sp in species_ls if sp not in leaves_input_tree]
        sys.stderr.write(f"ERROR: the following taxons are in the dist file but not in the provided tree: {', '.join(taxon_not_in_dist)}")
        sys.exit()
    if len(leaves_input_tree) != len(leaves_pruned_tree):
        inputpath = pathlib.Path(input_filename)
        pruned_tree_file = f"{inputpath.stem}.pruned.{inputpath.suffix}"
        mytree.write(format = 1, outfile = pruned_tree_file)
        sys.stderr.write(f"""
        WARNING: the tree you provided had taxons absent from the dist file.
        The pruned tree will be saved as {pruned_tree_file}
        """)
        
        
def get_threshold(species_threshold, chromfiles_list):
    """
    adjusts threshold of species
    :param species_threshold: an integer
    :param chromfiles_list: a list
    :return: an integer
    """
    nb_species = len(chromfiles_list)
    if 1 > species_threshold:
        sys.stderr.write(f'WARNING: species_threshold < 2. threshold will be changed to 2\n')
        species_threshold = 1
    elif species_threshold > nb_species:
        sys.stderr.write(f'WARNING: species_threshold > total number of species. threshold will be changed to {nb_species}\n')
        species_threshold = nb_species
    return species_threshold


def make_output_folder(name):
    """
    stop script if name exists
    :param name: folder name
    """
    name_path = pathlib.Path(name)
    if name_path.is_dir():
        sys.stderr.write("ERROR: folder already exists")
        sys.exit()
    else:
        name_path.mkdir()
        sys.stderr.write(f"Files will be saved in folder {name_path.absolute()}\n")


def save_chrom_data(chrom_data, filepath):
    """
    pickles the chrom_data to a filepath
    """
    output_path = pathlib.Path(filepath)
    sys.stderr.write(f'Chrom data will be pickled to {output_path.absolute()}\n')
    with open(output_path, "wb") as f:
        pickle.dump(chrom_data,
                    f,
                    pickle.HIGHEST_PROTOCOL)


def load_chrom_data(filepath):
    """
    unpickles the chrom_data returns a nested_dict()
    """
    with open(filepath, "rb") as f:
        return pickle.load(f)


def filter_OGs(chrom_data, max_paralogs):
    """
    Edits in-place the nested dictionary
    species-specific OGs and OGs where at least one species exhibits more than max_paralogs
    :param chrom_data: output of parsers.read_chromfiles()
    :max_paralogs: an integer
    :return: nothing, edits chrom_data in-place
    """
    og_ls = list(chrom_data.keys())
    for og in og_ls:
        species_ls = chrom_data[og].keys()
        paralog_counts = []
        for species in species_ls:
            paralog_count = 0
            for scaffold in chrom_data[og][species].keys():
                paralog_count += len(chrom_data[og][species][scaffold].keys())
            paralog_counts.append(paralog_count)
        if max(paralog_counts) > max_paralogs:
            del chrom_data[og]


def get_min_distances_pair(og_pair, chrom_data, min_threshold, chromfiles_ls):
    """
    Generates distance matrix of pairwise OGs distances
    :param og_pair: tuple of ogs to check synteny for
    :param chrom_data: nested dict. returned by parsers.read_chromfiles
    :param min_threshold: minimum number of species in which a syntenic pair needs to be found to be retained
    :param chromfiles_ls: list of species
    :return: list with two og names and distances, matching order of chromfiles_ls
    """
    ogx, ogy = og_pair
    dist_ls = []
    for species in chromfiles_ls:
        chr1, chr2 = chrom_data[ogx][species].keys(), chrom_data[ogy][species].keys() #all chromosomes within OG1 and OG2 in species
        shared_chrom_ls = set(chr1) & set(chr2) #OG1 and OG2 are syntenic in species
        dist = math.inf
        for chromosome in shared_chrom_ls:
            ogx_genes_ls = chrom_data[ogx][species][chromosome].keys()
            ogy_genes_ls = chrom_data[ogy][species][chromosome].keys()
            for gene1, gene2 in itertools.product(ogx_genes_ls, ogy_genes_ls):
                pos1, *_ = chrom_data[ogx][species][chromosome][gene1]
                pos2, *_ = chrom_data[ogy][species][chromosome][gene2]
                pos1, pos2 = map(int, (pos1, pos2))
                d = abs(pos2 - pos1)
                overlap_genes = pos1 == pos2
                if overlap_genes is False:
                    if d < dist:
                        dist = d
        if math.isinf(dist):
            dist = None
        dist_ls.append(dist)
    nb_syntenic_sp = sum(x is not None for x in dist_ls)
    if nb_syntenic_sp >= min_threshold:
        return [ogx, ogy] + dist_ls


def window(seq, n = 2):
    """Returns a sliding window (of width n) over data from the iterable
    s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ..."""
    it = iter(seq)
    result = tuple(itertools.islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result