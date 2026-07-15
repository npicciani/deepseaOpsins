#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
List the query sequence names (leaves prefixed with 'QUERY___') found in a
detached opsin clade tree.

usage: list_opsin_queries.py <tree_in> <list_out>
"""
import re
import sys
from ete3 import Tree

tree_in, list_out = sys.argv[1:3]

tree = Tree(tree_in, format=5)
with open(list_out, "w") as out:
    for leaf in tree:
        if leaf.name.startswith("QUERY"):
            match = re.search(r"QUERY___(.+)", leaf.name)
            if match:
                out.write(match.group(1) + "\n")