#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Remove the characters RAxML adds to internal node labels, root the tree on a
given outgroup leaf, then detach and write out the clade defined by the
common ancestor of two reference leaves (i.e. the opsin clade).

usage: detach_opsin_clade.py <tree_in> <outgroup_leaf> <anc_leaf_1> <anc_leaf_2> <tree_out>
"""
import sys
from ete3 import Tree

tree_in, outgroup, anc_leaf_1, anc_leaf_2, tree_out = sys.argv[1:6]

tree = Tree(tree_in, format=5)
tree.set_outgroup(tree & outgroup)
opsin_ancestor = tree.get_common_ancestor(anc_leaf_1, anc_leaf_2)
opsin_clade = opsin_ancestor.detach()
opsin_clade.write(format=5, outfile=tree_out)