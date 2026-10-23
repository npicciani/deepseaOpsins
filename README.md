# deepseaOpsins

A project to profile sequence variation associated with depth in opsins from ctenophores and siphonophores.

## Layout

Each stage of the analysis is a self-contained Snakemake workflow in its own folder, run from inside that folder.

| folder | stage |
|---|---|
| [`data-mining/`](data-mining/) | pull opsins from translated transcriptomes (BLAST + phylogenetic placement, PIA approach) |
| [`scoring-sequences/`](scoring-sequences/) | keep opsins that span all seven TM helices (profile-HMM coverage + DeepTMHMM topology) |
