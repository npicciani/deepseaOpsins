# filtering-sequences: keep opsins that span all seven TM helices

Part of [deepseaOpsins](../README.md).

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

## How the filter works

`docs/method.md` is the full account -- what each criterion measures, why the
profile replaced Pfam 7tm_1, and the held-out validation. In short, a sequence
is kept only if **both** hold:

- **A. Profile-HMM coverage.** Aligned to the profile with `hmmalign` alongside
  bovine rhodopsin (UniProt P02699), each of rhodopsin's seven TM spans must be
  covered by real residues in >= `hmm_cov` of its match columns.
- **B. DeepTMHMM topology.** A predicted helix must overlap each of the seven
  expected positions.

## Setup

```bash
conda env create -f envs/filtering.yaml     # snakemake + hmmalign
```

## Running

```bash
cd filtering-sequences
conda activate filtering
snakemake --cores 4 -n      # dry run
snakemake --cores 4
```

Outputs per sample, under `results/<sample>/`:

| file | |
|---|---|
| `opsins_7tm.faa` | keep |
| `opsins_7tm_review.faa` | criteria disagree |
| `opsins_7tm_dropped.faa` | removed; `drop_reason` in the report says why |
| `opsins_7tm_no_prediction.faa` | no topology record |
| `opsins_7tm_report.csv` | per-sequence coverage, helix pattern, verdict |
| `work_7tm/` | de-gapped input and the `hmmalign` Stockholm |

Adding a dataset is one line under `samples:` in `config/config.yaml`.
