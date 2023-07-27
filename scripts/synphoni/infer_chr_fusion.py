# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@File       :  screen.py
@Time       :  2022/12/24
@Author     :  Yuanting Ma
@Version    :  2.1
@Site       :  https://github.com/YuantingMaSC
@Contact    :  yuantingma@189.cn
"""
from os import mkdir
from shutil import rmtree
from os.path import exists
import argparse
from copy import deepcopy


def power_sets_binary(items, min_len=1, max_len=0):
    res = []
    if max_len == 0:
        max_len = len(items)
    N = len(items)
    for i in range(2 ** N):
        zj = []
        for j in range(N):
            if (i >> j) % 2 == 1:
                zj.append(items[j])
        res.append(zj)
    return [i for i in res if max_len >= len(i) >= min_len]


def remove_sub_set(set_list):
    set_list_ = []
    for item in set_list:
        set_list_.append(sorted(item, reverse=False))
    filter_l = list(filter(lambda f: not any(set(f) < set(g) for g in set_list_), set_list_))
    res = []
    [res.append(_) for _ in filter_l if _ not in res]
    return res


class Node(object):
    def __init__(self, name):
        self.name = str(name).replace("(", "").replace(")", "").replace("'", "").replace(",", "_").replace(' ', '')
        self.parent = None
        self.property = 'root'
        self.gene_status = 'unknown'
        # self.left_leaf = None
        # self.right_leaf = None
        self.sons = []

    def create_leaf(self, name):
        if len(self.sons) <= 2:
            self.sons.append(name)
            self.sons[-1].parent = self
        else:
            raise ValueError("leaf of this node {} is more than 2 !".format(self.name))


def is_gene_exist(gene_list: list, data: dict):
    for refer_gene_list in list(data.values()):
        if set(gene_list).issubset(refer_gene_list):
            return True
    else:
        return False


class CreatureTree(object):
    def __init__(self, tree_path, work_space):
        self.tree_path = tree_path
        self.work_folder = work_space
        self.tree_list = eval(open(tree_path).readlines()[0].replace(";", ""))
        self.depth_dict = dict()
        self.max_depth = 0
        self.node_list = []
        self.root = Node("root")
        self.leafs_of_node = {}
        self.data = dict()

        self.root.parent = Node('root_parent')
        self.root.parent.property = 'rootp'
        self.root.parent.sons = [self.root]

    def count_max_depth(self):
        self.max_depth = max(list(self.depth_dict.values()))

    def decode_tree(self, node, tree_tuple, depth):
        if not isinstance(tree_tuple, tuple):
            node.property = "leaf"
            node.gene_status = 'known'
            self.node_list.append(tree_tuple)
            self.depth_dict[tree_tuple] = depth
            return tree_tuple, depth
        else:
            for item in tree_tuple:
                node.create_leaf(Node(item))
                if depth > 0:
                    node.property = 'parent'
                self.decode_tree(node.sons[-1], item, depth=depth + 1)

    def read_tree(self):
        path = self.tree_path
        tree = eval(open(path).readlines()[0].replace(";", ""))
        self.decode_tree(self.root, tree, 0)
        # print(tree)

    def process_raw_data(self, res_dir):
        temp_dict = dict()
        file_rows = open(res_dir, "r").readlines()
        for row in file_rows:
            row = row.replace("\n", "")
            name = row.split("\t")[0]
            ances = row.split("\t")[1:]
            temp_dict[name] = ances
        return temp_dict

    def prepare_data(self):
        self.read_tree()
        for creature_name in self.node_list:
            self.data[creature_name] = self.process_raw_data(
                "{}/{}_chr_breakage_fusion.result".format(self.work_folder, creature_name))

    def find_leafs_of_target_node(self, node, target_node):
        if len(node.sons) == 1:
            if node.sons[0].property == 'root':
                self.find_leafs_of_target_node(node.sons[0], target_node)
        else:
            if node.sons[0].property == "leaf":
                if target_node.name not in self.leafs_of_node.keys():
                    self.leafs_of_node[target_node.name] = []
                self.leafs_of_node[target_node.name].append(node.sons[0])
            elif node.sons[0].property != "leaf":
                self.find_leafs_of_target_node(node.sons[0], target_node)
            if node.sons[1].property == "leaf":
                if target_node.name not in self.leafs_of_node.keys():
                    self.leafs_of_node[target_node.name] = []
                self.leafs_of_node[target_node.name].append(node.sons[1])
            elif node.sons[1].property != "leaf":
                self.find_leafs_of_target_node(node.sons[1], target_node)

    def find_outgroup_nodes(self, node1, node2, inclued_target_nodes=False):
        if node1.parent.name != node2.parent.name:
            raise ValueError("{} and {} belongs to different parents !".format(node1.name, node2.name))
        if node1.parent.parent.name not in self.leafs_of_node.keys():
            self.find_leafs_of_target_node(node1.parent.parent, node1.parent.parent)
        if node1.parent.name not in self.leafs_of_node.keys():
            self.find_leafs_of_target_node(node1.parent, node1.parent)
        # 外群是父节点的父节点的另外一支的节点
        if not inclued_target_nodes:
            return list(set(self.leafs_of_node[node1.parent.parent.name]) - set(self.leafs_of_node[node1.parent.name]))
        else:
            return self.leafs_of_node[node1.parent.parent.name]

    def store_data(self, creature1, creature2, over_lap, mark):
        ances_name = "{}_{}".format(creature1, creature2)
        if ances_name not in self.data.keys():
            self.data[ances_name] = dict()
        self.data[ances_name]["chr_{}".format(mark)] = over_lap

    def store_data_v2(self, creature1, creature2, over_lap, mark):
        ances_name = "{}_{}".format(creature1, creature2)
        if ances_name not in self.data.keys():
            self.data[ances_name + "_v2"] = dict()
        self.data[ances_name + "_v2"]["chr_{}".format(mark)] = over_lap

    def recover_nodes(self, node):
        if node.property == "leaf":
            node.gene_status = 'known'
        else:
            node.gene_status = 'unknown'
            for son in node.sons:
                self.recover_nodes(son)

    def fuzzy_overlap(self, list1: list, list2: list):
        """
        subitem of item should be seperate by "_"
        :param list1:
        :param list2:
        :return:
        """
        fully_overlap = []
        semi_overlap = []
        pass

    def write_res(self, res, creature1, creature2, name):
        mark = 1
        for _ in res:
            self.store_data(creature1, creature2, _, mark)
            open("{}/{}/{}_{}_anc.result".format(self.work_folder, name, creature1, creature2), 'a').write(
                "{}\t{}\n".format(mark, ",".join(_)))
            mark += 1

    def match(self, creature1_node, creature2_node, name):
        creature1, creature2 = creature1_node.name, creature2_node.name
        c1_data = self.data[creature1]
        c2_data = self.data[creature2]

        temp_res = []
        for c1_anc in c1_data.values():
            for c2_anc in c2_data.values():
                over_lap = list(set(c1_anc) & set(c2_anc))
                # 如果正好能够找到完全重合的序列，就直接把完全重合的部分写到祖先节点的结果中
                if len(over_lap) <= 1:
                    # self.store_data(creature1, creature2, over_lap, mark)
                    continue
                elif len(over_lap) == max(len(c1_anc), len(c2_anc)):
                    # self.store_data(creature1, creature2, over_lap, mark)
                    temp_res.append(over_lap)
                    # open("{}/ancestor/{}_{}_anc.result".format(self.work_folder, creature1, creature2), 'a').write(
                    #     "chr_{}\t{}\n".format(mark, "\t".join(over_lap)))
                    # mark += 1
                # 否则挨个到外群查找不匹配的祖先基因，超过设定比例外群物种
                elif len(over_lap) >= 2:
                    temp_over_l = over_lap.copy()
                    dif_set_left = list(set(c1_anc) - (set(c1_anc) & set(c2_anc)))
                    dif_set_right = list(set(c2_anc) - (set(c1_anc) & set(c2_anc)))
                    outfroup_nodes = self.find_outgroup_nodes(creature1_node, creature2_node)
                    if len(outfroup_nodes) == 0:
                        temp_res.append(temp_over_l)
                        continue
                    cmb_subset_form = power_sets_binary(dif_set_left,min_len=0)+power_sets_binary(dif_set_right)
                    for anc_gene in cmb_subset_form:
                        anc_gene_express = 0
                        for outgroup_node in outfroup_nodes:
                            outgroup_node_data = self.data[outgroup_node.name]
                            if is_gene_exist(over_lap + anc_gene, outgroup_node_data):
                                anc_gene_express += 1
                        express_rate = anc_gene_express / len(outfroup_nodes)
                        if express_rate > 0.1:
                            temp_res.append(temp_over_l+anc_gene)
                    else:
                        temp_res.append(temp_over_l)

                    # self.store_data(creature1, creature2, over_lap, mark)
                    # open("{}/{}/{}_{}_anc.result".format(self.work_folder, name, creature1, creature2), 'a').write(
                    #     "chr_{}\t{}\n".format(mark, "\t".join(over_lap)))
                    # mark += 1
        # print(temp_res)
        res = remove_sub_set(temp_res)
        self.write_res(res, creature1, creature2, name)
        creature1_node.parent.gene_status = "known"

    def infer_ancestor(self, node_, depth=0):
        name = 'ancestor'
        if depth == 0:
            if not exists("{}/ancestor/".format(self.work_folder)):
                mkdir("{}/ancestor/".format(self.work_folder))
            else:
                rmtree("{}/{}/".format(self.work_folder, name))
                mkdir("{}/{}/".format(self.work_folder, name))
        if depth == 0 and len(self.data) == 0:
            self.prepare_data()
        if node_.sons[0].gene_status == "known" and node_.sons[1].gene_status == "known":
            print("Inferring the chromosomal rearrangements of {} and {}".format(node_.sons[0].name, node_.sons[1].name))
            self.match(node_.sons[0], node_.sons[1], name=name)
        else:
            for son_node in node_.sons:
                if son_node.property == 'leaf':
                    continue
                self.infer_ancestor(son_node, depth + 1)
            print("Inferring the chromosomal rearrangements of {} and {}".format(node_.sons[0].name, node_.sons[1].name))
            self.match(node_.sons[0], node_.sons[1], name=name)

    def match_v2(self, creature1_node, creature2_node, v2_name):
        if creature1_node.parent.name != creature2_node.parent.name:
            raise ValueError(
                "{} and {} belongs to different parents !".format(creature1_node.name, creature2_node.name))
        group_nodes = self.find_outgroup_nodes(creature1_node, creature2_node, inclued_target_nodes=True)
        combine_form = power_sets_binary(group_nodes, 2, 2)
        binary_combine_res_list = [[] for i in range(len(combine_form))]
        for form_st, cmb_index in zip(combine_form, range(len(combine_form))):
            creature1, creature2 = form_st[0].name, form_st[1].name
            c1_data = self.data[creature1]
            c2_data = self.data[creature2]
            # mark = 1
            for c1_anc in c1_data.values():
                for c2_anc in c2_data.values():
                    # 所有两两组合的形式下，每两个物种间可以得到很多的overlap，最后求所有组合形式的并集
                    over_lap = list(set(c1_anc) & set(c2_anc))
                    if len(over_lap) >= 2:
                        binary_combine_res_list[cmb_index].append(over_lap)
        temp = []
        for one_form_res in binary_combine_res_list:
            temp.extend(one_form_res)
        # print(temp)
        res = remove_sub_set(temp)
        self.write_res(res, creature1_node.name, creature2_node.name, v2_name)
        # mark = 1
        # for over_l in res:
        #     open("{}/{}/{}_{}_anc_v2.result".format(self.work_folder, v2_name, creature1_node.name, creature2_node.name), 'a').write(
        #         "{}\t{}\n".format(mark, "\t".join(over_l)))
        #     mark += 1

    def infer_ancestor_v2(self, node_, depth=0):
        v2_name = 'ancestor_v2'
        if depth == 0:
            if not exists("{}/{}/".format(self.work_folder, v2_name)):
                mkdir("{}/{}/".format(self.work_folder, v2_name))
            else:
                rmtree("{}/{}/".format(self.work_folder, v2_name))
                mkdir("{}/{}/".format(self.work_folder, v2_name))
        if depth == 0 and len(self.data) == 0:
            self.prepare_data()
        if node_.sons[0].gene_status == "known" and node_.sons[1].gene_status == "known":
            print("<<V2>>--Inferring the chromosomal rearrangements of {} and {}".format(node_.sons[0].name, node_.sons[1].name))
            self.match_v2(node_.sons[0], node_.sons[1], v2_name=v2_name)
        else:
            for son_node in node_.sons:
                if son_node.property == 'leaf':
                    continue
                self.infer_ancestor_v2(son_node, depth + 1)
            print("<<V2>>--Inferring the chromosomal rearrangements of {} and {}".format(node_.sons[0].name, node_.sons[1].name))
            self.match_v2(node_.sons[0], node_.sons[1], v2_name=v2_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree_dir", default="./data/tree.nwk", help="the path of tree.nwk")
    parser.add_argument("--work_space", default="./d", help="Chr folder path")

    args = parser.parse_args()

    a = CreatureTree(args.tree_dir, args.work_space)
    a.infer_ancestor(a.root)
    a.recover_nodes(a.root)
    a.infer_ancestor_v2(a.root)

