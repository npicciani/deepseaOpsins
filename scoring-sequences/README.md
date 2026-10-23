# scoring-sequences: score opsins for seven-TM completeness

Part of [deepseaOpsins](../README.md).

This workflow **measures**; it does not decide. For every sequence it writes a
report of the quantities the seven-TM filter is decided on -- profile coverage
of each TM window, DeepTMHMM's helix count and where the helices fall, the
N-terminal length -- and applies no threshold.
[`../filtering-sequences`](../filtering-sequences/) reads the report and makes
the keep / review / drop calls. Everything costly (hmmalign, DeepTMHMM on the
GPU) is here, so changing a threshold there never re-runs anything here.

## Resources

`resources/opsin_profile/` holds the clade-inclusive opsin profile HMM that
the filter measures TM coverage against, and everything needed to rebuild it:

| file | |
|---|---|
| `animalOpsins_090426.hmm` | the profile -- 23 seeds, 336 match states |
| `baits_090426.fasta` | the 23 seed sequences |
| `baits_090426.aln` | their E-INS-i alignment, which the profile is built from |
| `build_hmm.py` | alignment -> profile (pyhmmer, default Builder) |

The profile is a **versioned resource, not a pipeline product**. Editing the
seeds changes which alignment columns become match states, which can shift TM1
and TM7 placement, so a new profile needs revalidating on held-out sequences
before it is trusted. To rebuild:

```bash
mafft --ep 0 --genafpair --maxiterate 1000 baits_NEW.fasta > baits_NEW.aln   # E-INS-i
python build_hmm.py baits_NEW.aln -o animalOpsins_NEW.hmm
```

Keep dated filenames; never overwrite a profile in place.

## Test input

`input/pseudopsins_melatonin_trichoplax_outgroups_devivo2023.fasta` -- 30
sequences (pseudopsins, melatonin receptors, *Trichoplax* and other outgroups
from de Vivo et al. 2023). Small enough to run the whole workflow in seconds.

It contains two sequences entered twice under different headers:
`Pseudopsin_MolCephalopodaOctopusvulgaris_2483` = `Pseudopsin_XM029798206OcsiOCTOCE`
and `Pseudopsin_MolCephalopodaOctopusbimaculoides_2422` = `Pseudopsin_XM014924517OcbiOCTOCE`.

## What is measured

`docs/method.md` is the full account -- what each measurement captures, why the
profile replaced Pfam 7tm_1, and the held-out validation. In short:

- **A. Profile-HMM coverage.** Each sequence is aligned to the profile with
  `hmmalign` alongside bovine rhodopsin (UniProt P02699), and for each of
  rhodopsin's seven TM spans the report gives the fraction of its match
  columns covered by real residues.
- **B. DeepTMHMM topology.** The number of predicted helices, and which of the
  seven expected positions a predicted helix overlaps.

## Setup

```bash
conda env create -f envs/filtering.yaml     # snakemake + hmmalign
```

DeepTMHMM 1.0 is academically licensed and is **not** in this repo or the env.
It must be installed on the host with a launcher taking
`deeptmhmm <fasta> <outdir>`; point `deeptmhmm.bin` in `config/config.yaml` at
it. On deepedna2 it is installed system-wide at `/usr/local/bin/deeptmhmm`
(program in `/opt/deeptmhmm`), and uses the GPUs.

## Running

```bash
cd scoring-sequences
conda activate filtering
snakemake --cores 4 --resources gpu=2 -n      # dry run
snakemake --cores 4 --resources gpu=2
```

`gpu=2` lets two DeepTMHMM batches run at once, one per device listed under
`deeptmhmm.gpus`; each batch claims a free device through a lock file, so two
jobs never share a card. Snakemake treats an undeclared resource as
unlimited, so leave the flag off and batches run `--cores` wide.

Outputs per sample, under `results/<sample>/`:

| file | |
|---|---|
| `report.csv` | one row per scored sequence -- measurements only, no verdicts |
| `sequences.faa` | exactly the sequences in the report, gaps stripped |
| `report.json` | what the report was measured against: profile path and sha1, match columns per TM window, counts in and out |
| `work_7tm/` | de-gapped input and the `hmmalign` Stockholm |
| `dtm_out/TMRs.gff3` | this dataset's topology, assembled from the store |
| `dtm_plan/` | which sequences this run needed to predict (`plan.json`) |
| `dtm_runs/batch_XXX/` | raw DeepTMHMM output per batch, kept as a log |

### Report columns

| column | |
|---|---|
| `name`, `length` | FASTA header and de-gapped length |
| `cov_TM1` .. `cov_TM7`, `min_cov` | fraction of each window's profile match columns covered by residues |
| `has_topology` | 1 if DeepTMHMM has a prediction for this sequence; 0 (and the DeepTMHMM columns blank) if not, e.g. over `max_len` |
| `n_tm_helices` | DeepTMHMM's raw helix count |
| `windows_with_helix`, `helix_pattern` | how many / which of the seven windows a predicted helix overlaps, e.g. `1111110` |
| `n_before_TM1`, `n_after_TM7` | residues outside the 7TM core |
| `extra_helix_location` | where helices beyond seven fall: `N-terminal`, `interior` (the fusion signature), `C-terminal` |
| `has_signal_peptide` | DeepTMHMM calls a signal peptide |
| `res_at_K296` | residue in rhodopsin's retinal-binding lysine column; reported, never filtered on |

Coverage is written to six decimals. Each value is k/n with n of about 20-25,
so a threshold applied to the report decides exactly as on the unrounded value.

Adding a dataset is one line under `samples:` in `config/config.yaml`, plus the
matching entry in `../filtering-sequences/config/config.yaml`.

## DeepTMHMM: one prediction per sequence

```
input FASTA ──> dtm_plan (checkpoint: what is NOT already predicted?)
                   └──> dtm_batch × N (DeepTMHMM, 1 GPU each) ──> resources/dtm_cache/
input FASTA ──────────────────────> dtm_collect <────────────────────┘
                                         └──> dtm_out/TMRs.gff3 ──> score_7tm
```

DeepTMHMM predicts each sequence independently; a batch only amortises loading
ESM-1b. So predictions are cached **per sequence**, in `resources/dtm_cache/`,
keyed by a sha1 of the sequence itself (uppercased, gaps stripped) rather than
its header. That key is what DeepTMHMM actually saw, so:

- adding five sequences to a dataset predicts five, not a whole batch;
- a renamed header reuses its prediction, while a header reused for a changed
  sequence misses the cache -- where a name-keyed cache would serve a stale
  topology;
- the same sequence under two headers is predicted once and reported under
  both (the test input has two such pairs);
- the store is shared across samples, so overlap between datasets is
  predicted once.

`dtm_plan` is a checkpoint because the work set depends on the store at run
time. The cost is that `-n` shows the checkpoint first and the batch jobs only
after it has run.

Every record carries a `source` tag (`deeptmhmm.source` in the config), so a
store holding predictions from more than one install says where each came
from.

**Long sequences.** ESM-1b takes at most 1,022 residues; DeepTMHMM embeds longer
sequences in consecutive 1,022-residue windows rather than truncating them, so
the whole sequence is predicted. Anything over `max_len` is skipped: it is
still scored for profile coverage, with `has_topology = 0`, and the filter
classes it `no_prediction`.

### Recovery

`dtm_batch.py` exits non-zero on failure so snakemake sees it. Each record is
written atomically and the batch sentinel is written last, so a batch killed
midway keeps the sequences it finished, and the next run predicts only the
rest.

A `dtm_store.FormatError` means batch output did not look like DeepTMHMM
wrote it. It is fatal on purpose: a silently mis-parsed topology would enter
the store and be trusted by every later run.

### Importing existing DeepTMHMM output

`dtm_ingest_legacy.py` loads batch directories produced outside the workflow
-- including earlier runs on the BioLib web service -- so existing predictions
are reused. No FASTA is needed; `predicted_topologies.3line` carries the
sequences.

```bash
python workflow/scripts/dtm_ingest_legacy.py -n <batch dirs...> --cache-dir resources/dtm_cache
```

Imported records are tagged `biolib:DTU/DeepTMHMM` unless `--source` says
otherwise. It never overwrites a stored record: a sequence already stored with
a *different* prediction is reported as a conflict (non-zero exit), since two
runs disagreeing means the model differed between them.

On the 30-sequence test input, the local install on deepedna2 and the BioLib
service produced byte-identical `TMRs.gff3` output.
