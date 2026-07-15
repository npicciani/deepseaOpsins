# deepseaOpsins


A project to profile sequence variation associated with depth in opsins from ctenophores and siphonophores

July 14: created repo to add snakemake rules for pulling opsins from translated transcriptomes in moredata/isoseq


## First step: 

### Create a conda environment with the tools for running the main script:

``` conda create -n pia_snakemake -f pia_snakemake.yaml ```


Rule chain: make_blastdb → blastp_opsins → extract_blast_hits → extract_candidate_opsins → align_candidate_opsins → aln_to_phylip → raxml_epa → clean_raxml_tree → detach_opsin_clade → list_opsin_queries → get_opsins