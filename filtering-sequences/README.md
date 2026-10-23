# filtering-sequences: keep opsins that span all seven TM helices

Part of [deepseaOpsins](../README.md).

This workflow **decides**. It reads the report written by
[`../scoring-sequences`](../scoring-sequences/) -- per-sequence profile
coverage and DeepTMHMM topology, with no threshold applied -- and makes the
keep / review / drop calls. It runs no alignment and no DeepTMHMM, needs no
GPU and no licensed software, and finishes in seconds, so thresholds can be
explored freely: change one in `config/config.yaml` and only this workflow
re-runs.

## The rules

A sequence is kept only if **both** criteria pass:

- **A. Profile coverage.** `min_cov >= hmm_cov`: every one of rhodopsin's
  seven TM windows is covered by residues in at least that fraction of its
  profile match columns.
- **B. DeepTMHMM topology.** With `dtm_mode: position` (default), a predicted
  helix overlaps all seven windows (`helix_pattern` is `1111111`). With
  `dtm_mode: count`, at least `min_tm` helices are predicted -- cheaper, but it
  cannot tell a missing TM1 from a complete core.

and none of these drop reasons applies:

| `drop_reason` | |
|---|---|
| `<7_TM` | both criteria call the 7TM core incomplete |
| `>N_TM` | DeepTMHMM predicts more than `max_tm` helices (`0` disables) |
| `>Naa_before_TM1` | more than `max_nterm` residues precede TM1 -- the fusion / mis-assembly signature (`0` disables) |

A sequence tripping several rules gets every reason, joined by `+`.

| verdict | when |
|---|---|
| `no_prediction` | no DeepTMHMM record (`has_topology = 0`) -- a data gap, never a failure |
| `drop` | any drop reason |
| `keep` | both criteria pass |
| `review` | the criteria disagree |

The scientific rationale for both criteria, and the validation, is in
[`../scoring-sequences/docs/method.md`](../scoring-sequences/docs/method.md).

## Running

Score first, then filter:

```bash
conda activate filtering           # any env with snakemake; scripts are stdlib-only
cd scoring-sequences  && snakemake --cores 4 --resources gpu=2
cd ../filtering-sequences && snakemake --cores 1
```

Each sample in `config/config.yaml` names the three files one scoring run
produced for it -- `report.csv`, `sequences.faa`, `report.json`. Adding a
dataset means adding it to both workflows' configs.

Outputs per sample, under `results/<sample>/`:

| file | |
|---|---|
| `opsins_7tm.faa` | keep |
| `opsins_7tm_review.faa` | criteria disagree |
| `opsins_7tm_dropped.faa` | removed; `drop_reason` says why |
| `opsins_7tm_no_prediction.faa` | no topology record |
| `opsins_7tm_verdicts.csv` | the scoring report plus `pass_hmm`, `pass_deeptmhmm`, `drop_reason`, `verdict` |
| `opsins_7tm_params.json` | the thresholds applied, and the scoring run's metadata (profile sha1, match columns per window) |
| `filter_summary.tsv` | counts by verdict, by drop reason, and by criterion |

Every FASTA is always written; an empty one means that class is empty.

In `filter_summary.tsv` the `drop_reason` block partitions the drops
exclusively; the `criterion_met` block counts sequences meeting each rule and
therefore **overlaps** -- a sequence tripping two rules is counted under both.
