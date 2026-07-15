#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert an aligned fasta file to relaxed phylip format, as required by RAxML.

usage: fasta_to_phylip.py <input.aln> <output.phy>
"""
import sys
from Bio import AlignIO

aln_in, phy_out = sys.argv[1:3]

alignment = AlignIO.read(aln_in, "fasta")
n_seqs = len(alignment)
length = alignment.get_alignment_length()

with open(phy_out, "w") as out:
    out.write(f"{n_seqs}\t{length}\n")
    for record in alignment:
        out.write(f"{record.id}\t{record.seq}\n")