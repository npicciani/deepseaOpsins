# deepseaOpsins

A project to profile sequence variation associated with depth in opsins from ctenophores and siphonophores.

**July 14:** created repo to add Snakemake rules for pulling opsins from translated transcriptomes in `moredata/isoseq`.

## Setup

Create a conda environment with the tools needed to run the pipeline:

```
conda env create -f pia_snakemake.yaml -n pia_snakemake
conda activate pia_snakemake
```

## Running the pipeline

```
snakemake --cores all
```

Snakemake will build every `results/opsins/{transcriptome}.opsins.fasta` requested for the translated transcriptomes (`{transcriptome}.pep`) placed in the folder `resources`.

## Pipeline overview

The pipeline identifies opsin sequences in a translated transcriptome by combining a BLAST search against known opsin sequences with phylogenetic placement on a reference opsin tree (adapted from the PIA approach in Speiser et al. 2010; reference tree from Picciani et al. 2018). This gives higher-confidence opsin calls than a BLAST hit alone, since a candidate has to also land inside the opsin clade of the reference tree.

Each step below is its own Snakemake rule, chained in this order:

`make_blastdb → blastp_opsins → extract_blast_hits → extract_candidate_opsins → align_candidate_opsins → aln_to_phylip → raxml_epa → clean_raxml_tree → detach_opsin_clade → list_opsin_queries → get_opsins`

1. **make_blastdb** — Builds a BLAST protein database from the transcriptome's translated ORFs (`{transcriptome}.pep`).

2. **blastp_opsins** — Searches that database with a curated set of opsin "bait" sequences (`resources/baits0206.fasta`) to find candidate opsin hits (blastp, e-value cutoff 1e-10).

3. **extract_blast_hits** — Parses the blastp results table and produces a deduplicated list of candidate opsin sequence IDs.

4. **extract_candidate_opsins** — Pulls the full peptide sequences for those candidate IDs out of the transcriptome's ORF file.

5. **align_candidate_opsins** — Adds the candidate opsin sequences to a pre-built reference opsin alignment using MAFFT.

6. **aln_to_phylip** — Converts the resulting alignment from FASTA to relaxed PHYLIP format, as required by RAxML.

7. **raxml_epa** — Places the candidate sequences onto a fixed reference opsin phylogeny (`resources/rep17.treefile`) using RAxML's Evolutionary Placement Algorithm, producing a labelled tree.

8. **clean_raxml_tree** — Strips the internal node labels RAxML adds to the placement tree, so it can be parsed cleanly downstream.

9. **detach_opsin_clade** — Roots the tree on a designated outgroup and detaches the clade containing the known reference opsins, isolating which placements fall within it.

10. **list_opsin_queries** — Lists the names of the candidate sequences that were placed inside the detached opsin clade — i.e., the confirmed opsins.

11. **get_opsins** — Pulls those confirmed opsin sequences (by name) out of the original transcriptome, producing the final `{transcriptome}.opsins.fasta` output.

## Reference resources

Fixed across all transcriptomes (`resources`):

- `baits0206.fasta` — opsin bait sequences for the BLAST search
- `all_0512_gt1_rs50_s65.aln` — reference opsin alignment
- `rep17.treefile` — reference opsin phylogeny used for placement