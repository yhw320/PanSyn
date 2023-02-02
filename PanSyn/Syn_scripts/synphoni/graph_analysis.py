import networkx as nx
import statistics
import itertools
import csv
import synphoni.utils as su
from networkx.algorithms.community import k_clique_communities


def load_edgelist(filepath, nmax):
    """load_edgelist for a file
    :param filepath: csv file with three named fields, which will be renamed source, target and weight
    original field names will be ignored
    :param nmax: edges above this threshold will be filetered out
    :return: two edgelists with data in a list, parseable by nx.parse_edgelist
    """
    edge_list = []
    edge_list_filt = []
    with open(filepath, "r") as f:
        myreader = csv.reader(f)
        _ = next(myreader) #skip header
        for row in myreader:
            weight = float(row[2])
            row_str = " ".join(row)
            edge_list.append(row_str)
            if weight <= nmax:
                edge_list_filt.append(row_str)
    return edge_list, edge_list_filt


def weighted_edgelist_to_graph(edge_list):
    """weighted_edgelist_to_graph reads edgelists with data in alist
    :param edge_list: list of strings, edge data space separated
    :return: nx.Graph object
    """
    return nx.parse_edgelist(edge_list, data=(("weight", float),))


def refine_component_to_cliques(raw_graph, filtered_graph):
    """refine_component_to_cliques 
    :param raw_graph: unfiltered graph, all edges are there
    :param filtered_graph: graph where edges > nmax are filtered
    :return: a generator, each item is a set of OGs
    """
    components = sorted(list(nx.connected_components(filtered_graph)), key = len, reverse = True)
    for comp in components:
        cliques_to_keep = []
        G_comp = nx.subgraph(raw_graph, comp)#all links, assume that a path goes through
        max_cliques = sorted(list(nx.find_cliques(G_comp)),
                             key = len,
                             reverse = True)
        max_clique_size = len(max_cliques[0])        
        longest_cliques = [x for x in max_cliques if len(x) > 2 and len(x) == max_clique_size] 
        if len(longest_cliques) == 1: #only one longest string
            cliques_to_keep.append(set(longest_cliques[0]))
        elif len(longest_cliques) > 1: #multiple longest clique, keep the most compact
            most_compact_longest_clique, median_dist_to_keep = None, float("inf")
            for clique in longest_cliques:
                G_clique = nx.subgraph(G_comp, clique)
                G_clique_median_dist = statistics.median([x[2]["weight"] for x in G_clique.edges.data()])
                if median_dist_to_keep > G_clique_median_dist:
                    median_dist_to_keep = G_clique_median_dist
                    most_compact_longest_clique = clique
            if most_compact_longest_clique != None:
                cliques_to_keep.append(set(most_compact_longest_clique))               
        sorted_cliques_w_o_longest = [c for c in max_cliques if set(c) not in cliques_to_keep]
        for clique in sorted_cliques_w_o_longest:
            clique = set(clique)
            shared_ogs_w_kept = set()
            if len(clique) > 2:
                for kept_clique in cliques_to_keep:
                    shared_ogs_w_kept.update(kept_clique & set(clique))
                if len(shared_ogs_w_kept) in [0,1]:
                    cliques_to_keep.append(clique)
            else:
                G_clique_microsynt = nx.subgraph(filtered_graph, clique)
                small_clique_comp = list(nx.connected_components(G_clique_microsynt))#2-clique could be artifact from 2 distant points that we join with a link above nmax
                self_edges = any([source == target for source, target in filtered_graph.edges])
                for comp_clique in small_clique_comp: # in the event that there's two nodes with self edges, that end up connected by a spurious link above nmax
                    for kept_clique in cliques_to_keep:
                        shared_ogs_w_kept.update(kept_clique & set(comp_clique))
                    if self_edges and len(shared_ogs_w_kept) == 0:
                        cliques_to_keep.append(comp_clique)
        for clique in cliques_to_keep:
            yield clique


def genome_location_ogs(og_community, chrom_data, species_list, orthology, min_og_commu = 0.3):
    """
    loads the genome locations of the ogs of the specified community into a dict
    :param og_community: community of ogs
    :param chrom_data: nested defaultdict, output of scripts.parsers.read_chromfiles
    :param species_ls: list/tuple of prefixes of species included.
    :param orthology: dict where keys are accessions, values are ogs
    :param min_og_commu: 
    :returns: a dict. Keys are (species, chrom_name) tuples, values are dicts with following keys:
    acc_ls => list of accessions
    pos_ls => list of start/ends of genes on the chromosome
    og_ls => list of the orthogroups each protein is assigned to
    s_ls => list of strands each gene is encoded by
    same gene share indices across these lists
    """
    tmp_dict = dict()
    #creates a dict, keys  are (species, scaffold) tuples, values are dicts pointing to accessions pos og and strands
    for my_og in og_community:
        species_with_og = [sp for sp in species_list if sp in chrom_data[my_og].keys()]
        for species in species_with_og:
            for scaffold in chrom_data[my_og][species].keys():
                sp_scaff = (species, scaffold)
                accessions_ls = list(chrom_data[my_og][species][scaffold].keys())
                orthogroup_ls = [orthology[acc] for acc in accessions_ls]
                positions_ls = su.flatten([x[1:3] for x in chrom_data[my_og][species][scaffold].values()])
                strand_ls = [x[3] for x in chrom_data[my_og][species][scaffold].values()]
                if tmp_dict.get(sp_scaff) is None:
                    tmp_dict[sp_scaff] = dict(acc_ls = accessions_ls,
                                                 pos_ls = positions_ls,
                                                 og_ls = orthogroup_ls,
                                                 s_ls = strand_ls)
                else:
                    tmp_dict[sp_scaff]["acc_ls"].extend(accessions_ls)
                    tmp_dict[sp_scaff]["pos_ls"].extend(positions_ls)
                    tmp_dict[sp_scaff]["og_ls"].extend(orthogroup_ls)
                    tmp_dict[sp_scaff]["s_ls"].extend(strand_ls)
    output_dict = {k: v for k, v in tmp_dict.items() if
                    (len(v["og_ls"]) / len(og_community)) > min_og_commu}
    return output_dict


def og_info_to_graph(genome_location_orthogroups, fullgraph_ogs_filt, min_len = 3, min_shared = .5):
    """
    Builds a networkx.Graph object where vertices are scaffolds,
        edges are inter-species shared og occurences (>= min len).
        overlap coefficient should be of 0.5 at least
    :param genome_location_orthogroups: nested dict output of genome_location_ogs
    :param fullgraph_ogs_filt: graph loaded from edgelists created by step2 (only).
    used only to see if there's self edges for ogs of interest. If yes, multiple shared paralogs
    will be considered as multiple shared occurences. If not, paralogs will be considered as single shared og occurence
    :param min_len: minimum number of shared og occurences to make an inter-scaffold edge
    :param min shared: minimum overlap coefficient of shared ogs
    :returns: a networkx Graph object if there's at least 3 scaffolds in 3 species
    """
    G = nx.Graph()
    #avoid intraspecies link, replace connected components by 3-clique percolation
    for scaff1, scaff2 in itertools.combinations(genome_location_orthogroups.keys(), 2):
        shared1 = {og: genome_location_orthogroups[scaff1]['og_ls'].count(og) for og in genome_location_orthogroups[scaff1]['og_ls']}
        shared2 = {og: genome_location_orthogroups[scaff2]['og_ls'].count(og) for og in genome_location_orthogroups[scaff2]['og_ls']}
        shared_ogs = set(shared1.keys()) & set(shared2.keys())
        protoblock_shares = {og:min(shared1[og], shared2[og]) for og in shared_ogs}
        for og in protoblock_shares.keys(): #check for self edges. If there's no self edge in the dist network, the paralog shared will be counted as one.
            G_single_og = nx.subgraph(fullgraph_ogs_filt, og)
            og_has_no_self_edge = len(G_single_og.edges()) == 0
            #add a calculation to limit overlap coefficient. Count coefficient length based on value in of protoblock shares
            if og_has_no_self_edge:
                protoblock_shares[og] = 1 #this way paralogous blocks are only counted if there was a self edge in the OG graph
        smallest_scaffold = min(len(shared1), len(shared2))
        len_above_min_len = sum(protoblock_shares.values()) >= min_len
        overlap_above_min_shared = sum(protoblock_shares.values()) >= (min_shared * smallest_scaffold)
        if len_above_min_len and overlap_above_min_shared:
            G.add_edge(scaff1, scaff2)
    if len(G.nodes()) < 3:
        return None
    else:
        return G


def write_blocks(blocks_writer,
                 multi_sp_writer,
                 genome_location_ogs_dict,
                 og_info_graph,
                 k_perco,
                 known_dict):
    """
    write block information specified to a file
    takes as input a scaffold graph (og_info_graph) and the info about the gene it bears
    (genome_location_ogs_dict)
    Multi species blocks are found by 3-clique percolation
    (each block has at least 2 links to two other speciesm since graph has no
    intra-species links)
    :multi_sp_id and  blocks_id should be global variable where current multi_sp id and current block id is stored are stored
    :param blocks_writer: csv writer object writing blocks information
    one line per syntenic unit of each species (multiple occurences by species allowed)
    :param multi_sp_writer: csv writer object writing to multi sp information
    :param genome_location_ogs_dict: dictionary, output of genome_location_ogs function
    :param og_info_graph: nx.Graph object, output of og_info_to_graph
    :param k_perco: size of teh clique to use for percolation
    :returns: block_id dict
    """
    block_id_dict = {}        
    block_clusters = list(k_clique_communities(og_info_graph, k_perco))
    block_clusters = sorted(block_clusters, key = lambda x: len(x), reverse = True)
    for block_cluster in k_clique_communities(og_info_graph, k_perco):
        if len(known_dict.keys()) == 0 and len(block_id_dict.keys()) == 0:
            multi_sp_id = 1
            block_id = 0
        elif len(block_id_dict.keys()) == 0:
            multi_sp_id = max(known_dict.keys()) + 1
            block_id = max(itertools.chain.from_iterable(known_dict.values()))
        elif len(block_id_dict.keys()) != 0:
            multi_sp_id = max(block_id_dict.keys()) + 1
            block_id = max(itertools.chain.from_iterable(block_id_dict.values()))
        last_block_id = block_id + len(block_cluster) + 1
        block_id_ls = list(range(block_id + 1, last_block_id))
        block_id_dict[multi_sp_id] = block_id_ls
        cluster_sp_ls = set([chrom[0] for chrom in block_cluster])
        multi_sp_writer.writerow([multi_sp_id] + block_id_ls)
        for species, chromosome in block_cluster:
            pb = genome_location_ogs_dict[(species, chromosome)]
            block_id += 1
            num_links = len(block_cluster)
            linked_blocks = set(block_id_ls) - {str(block_id)}
            linked_blocks = ','.join(map(str, linked_blocks))
            num_linked_species = len(cluster_sp_ls)
            linked_species = ','.join(cluster_sp_ls)
            strand_ls = pb['s_ls']
            if len(set(strand_ls)) > 1:
                strand = '.'
            else:
                strand = strand_ls[0]
            coordinates_block = [int(x) for x in pb['pos_ls']]
            start_block = min(coordinates_block)
            end_block = max(coordinates_block)
            location = f'{chromosome}:{start_block}..{end_block}'
            length_block = end_block - start_block
            joined_acc = ','.join(pb['acc_ls'])
            joined_og = ','.join(pb['og_ls'])
            output_list = [block_id,
                           species,
                           num_links,
                           linked_blocks,
                           num_linked_species,
                           linked_species,
                           strand,
                           location,
                           length_block,
                           joined_acc,
                           joined_og]
            blocks_writer.writerow(output_list)
    return block_id_dict