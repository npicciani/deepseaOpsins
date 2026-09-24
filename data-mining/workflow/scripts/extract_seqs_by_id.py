#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract fasta records whose IDs appear in a text file (one ID per line).

Used twice in the opsin-identification pipeline: once to pull blastp hit
sequences out of the ORF peptide file, and once to pull the final query
sequences that fell in the opsin clade out of the same file.

usage: extract_seqs_by_id.py <id_list.txt> <input.fasta> <output.fasta>
"""
import sys
from Bio import SeqIO

id_list_file, fasta_in, fasta_out = sys.argv[1:4]

wanted = set()
with open(id_list_file) as f:
    for line in f:
        line = line.strip()
        if line:
            wanted.add(line)

records = [rec for rec in SeqIO.parse(fasta_in, "fasta") if rec.id in wanted]
SeqIO.write(records, fasta_out, "fasta")