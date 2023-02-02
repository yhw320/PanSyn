import synphoni.utils as su
import sys


def read_orthology(clus_file):
    """
    Reads an orthology file in the .clus format.
    Returns a dict such as key:value is accession: OG_name
    """
    output_dict = {}
    sys.stderr.write(f'Reading file {clus_file}\n')
    with open (clus_file, 'r') as f:
        for line in f:
            og,_, *acc_ls = line.rstrip().split('\t')
            output_dict |= {acc: og for acc in acc_ls}
    sys.stderr.write(f'Orthology file loaded, {len(set(output_dict.values()))} OGs found.\n')
    return output_dict


def read_chromfiles(chromfiles_list, orthology_dict, min_species):
    """
    parses chromfiles into nested dictionary
    :param chromfiles_list: list of files, chrom format, prefix of file should be the species PREFIX ()
    :param orthology_dict: output of parsers.read_orthology
    :return: a nested_dict qs follows
        [OG] [Species] [chromosome/scaffold] [accession] [coordinates]
        Coordinates is a tuple of len 3:
            (position on the chrom (1-based) , start, stop)
    """
    output_dict = su.nested_dict()
    species_ls = [name.split('/')[-1].split('.')[0] for name in chromfiles_list]
    for species, file in zip(species_ls, chromfiles_list):
        sys.stderr.write(f'Reading file {file}\n')
        tmp_coords = su.nested_dict()
        with open(file, 'r') as chrom_file:
            for line in chrom_file:
                _, accession, scaffold, strand, start, stop = line.rstrip().split('\t')
                if accession in orthology_dict.keys():
                    tmp_coords[scaffold][accession] = [start, stop, strand]
        for chromosome in tmp_coords.keys():
            chrsort = sorted((tmp_coords[chromosome]).keys(),key = lambda x: int((tmp_coords[chromosome])[x][0]))
            pos = 0
            for gene in chrsort:
                pos = pos + 1
                output_dict[orthology_dict[gene]][species][chromosome][gene]= [pos] + tmp_coords[chromosome][gene]
    sys.stderr.write('Chroms loaded!\n')
    output_dict = {k: v for k,v in output_dict.items() if len(v)>= min_species}
    return output_dict, species_ls
    