#! /usr/bin/env python3 
import numpy as np
import argparse
import csv
import networkx as nx
import pandas as pd
from pathlib import Path

parser = argparse.ArgumentParser(description ="""
Takes as input a folder with step2 results of SYNPHONI pipeline
folder should be nested with node/m_n/*dist with only one dist file per folder (n being a int)
""")
parser.add_argument("input",
                    help = "path where the results are stored. Dist files will be globbed according to folder structure",
                    type = str)
parser.add_argument("--output",
                    "-o",
                    help = "prefix of the output tables, One tab",
                    type = str,
                    required = True)
parser.add_argument("--nmax",
                    "-n",
                    help = "values for the nmax range should be start, stop, step, default is 1,100,1",
                    type = int,
                    nargs = 3,
                    default = [1, 100, 1])
parser.add_argument("--node_name",
                    "-nn",
                    help = "node_names: space-separated names of the nodes for which you want to build ancestral OGs networks",
                    type = str,
                    required=True)
args = parser.parse_args()


def inflection_ede(x,y, index):
    """python implementation of ede function from R package inflection
    :param x: values (should be ordered in increasing order)
    :param y: values
    :param index: curve passing through points is s-shaped and increasing (0) or decreasing (1)
    :return: tuple with x-left, x-right and inflexion point of the curve
        - xleft and x right are xvalues that are unique unconstrained extremes:
    """
    def diff(x, lag = 1):
        """
        implementation of the diff function of R
        Lagged differences of order 1
        :param x: a np.array
        :param lag: lag of the differences
        :return: vector of differences
        """
        n = len(x)
        return x[lag:n] - x[0:(n-lag)]
    def inflection_lin2(x1, y1, x2, y2, x):
        """python implementation of lin2 function from R package inflection 
        """
        return y1 + (y2 - y1) * (x - x1)/(x2 - x1)
    n = len(x)
    if (index == 1):
        y = -y 
    if n >= 4:
        LF = y - inflection_lin2(x[0], y[0], x[n - 1], y[n - 1], x)
        jf1 = LF.argmin()
        xf1 = x[jf1]
        jf2 = LF.argmax()
        xf2 = x[jf2]
        if jf2 < jf1: xfx = np.nan
        else: xfx = (xf1 + xf2) / 2
    else:
        xf1, xf2, xfx = np.nan, np.nan, np.nan
    return xf1, xf2, xfx


def get_dist(input_dir, node_dir, m_dir):
    """
    check whether a filepath is ok
    one dist file per folder
    dist file should have three fields
    """
    file_ls = list((input_dir/node_dir/m_dir).glob("*.dist"))
    if len(file_ls) == 0:
        return None
    elif len(file_ls) > 1:
        parent = file_ls[0].parent
        raise ValueError(f"the folder {parent} contains more than one file, please check")
    data = pd.read_csv(file_ls[0], nrows = 1)
    if len(data.columns) > 3:
        raise ValueError(f"the file {file_ls[0]} contains more than 3 fields, please check format")
    return file_ls[0]
        

def files_to_process(input_dir):
    """
    list files to process from a string popintting to an input directory
    """
    output = []
    folder = Path(input_dir)
    node_ls = [x.stem for x in folder.glob("*")]
    for node in node_ls:
        m_ls = [x.stem for x in (folder / args.node_name).glob("*")]
        for m_folder in m_ls:
            m = m_folder.split("_")[1]
            filepath = get_dist(folder, node, m_folder)
            if filepath != None:
                output.append((filepath, node, m))
    return output


def slice_network(graph, threshold, copy = True):
    """
    Remove all edges with weight > threshold from graph
    if copy is True, returns a copy, else, graph is edited inplace
    """
    trimmed_graph = graph.copy() if copy else graph
    trimmed_graph.remove_edges_from([(n1, n2) for n1, n2, w in trimmed_graph.edges(data = "weight") if w > threshold])
    if copy: return trimmed_graph


def reverse_range(start, stop, step):
    """build_range object to iterate drecreasingly
    makes it so that range includes last value specified as stop
    """
    start = start - 1
    step = -step
    return range(stop, start, step)


def largest_comp_analysis(file_ls, start,stop, step):
    mynmax = [start,stop, step]
    component_values = [["node", "m", "nmax", "proportion"]]
    inflection_points = [["node", "m", "x-left"]]
    for file, node, m in file_ls:
        print(f"processing {file}")
        df = pd.read_csv(file, header = 0,
                               names = ["source", "target", "weight"],
                               dtype = {"source": "str", "target": "str", "weight": "float64"})
        G = nx.convert_matrix.from_pandas_edgelist(df, edge_attr = "weight")
        nmax_ls = list(reverse_range(*mynmax))
        proportion_ls = []
        for nmax in nmax_ls:
            slice_network(graph = G,
                          threshold = nmax,
                          copy = False)
            c_comp = sorted(list(nx.connected_components(G)), key = len, reverse = True)
            if len(c_comp) > 0:
                proportion = len(c_comp[0]) /len(G.nodes())
            else:
                proportion = 0
            proportion_ls.insert(0, proportion)
            component_values.append([node, m, nmax, proportion])
        nmax_ls.reverse()
        inflection_points.append([node,
                                  m,
                                  inflection_ede(x = np.array(nmax_ls),
                                                 y = np.array(proportion_ls),
                                                 index = 0)[0]])
    return component_values, inflection_points


file_list = files_to_process(args.input)
component_values, inflection_points = largest_comp_analysis(file_list, *args.nmax)

with open(f"{args.output}.component.csv", "w") as f:
    component_writer = csv.writer(f)
    component_writer.writerows(component_values)

with open(f"{args.output}.inflection.csv", "w") as f:
    inflection_writer = csv.writer(f)
    inflection_writer.writerows(inflection_points)