> **Record of the method as developed, before it became a workflow.** File names in
> *Contents* and *Running it* refer to that analysis directory; several
> (`seed_qc_090426.csv`, `combined_081426_nocontam.fasta`, `dtm_out/`) are not in
> this repo. For how to run the workflow, see [`../README.md`](../README.md).

# 7-TM filtering with a clade-inclusive opsin profile

In this version of the filtering pipeline, **profile**: stock
Pfam `7tm_1` (PF00001) is replaced by a profile built from a curated animal
opsin bait set, which fixes a systematic under-call of TM7 in ctenophores.

**DeepTMHMM was not re-run.** `dtm_out/` is a copy of
`../curated/work/dtm_out/` — 13 batches covering all 2,451 sequences. Criterion B
was never the failing half, so there is nothing to redo.

## Contents

| file | role |
|---|---|
| `filter_7tm.py` | the filter; profile and quality exclusions are parameters |
| `build_hmm.py` | rebuild the profile from a bait alignment (pyhmmer, default Builder) |
| **`animalOpsins_090426.hmm`** | **current profile** — 23 seeds, 336 match states, EFFN 2.54 |
| **`baits_090426.fasta`** | **current bait set** — 23 seeds |
| **`baits_090426.aln`** | E-INS-i alignment the current profile was built from |
| `seed_qc_090426.csv` | per-seed diagnostic residues; all 23 carry K296 |
| `ref/P02699.json` | bovine rhodopsin sequence + annotated TM spans |
| `dtm_out/` | existing DeepTMHMM output, copied, not regenerated |
| `combined_081426_nocontam.fasta` | input sequences (2,451) |
| `opsins_7tm*.faa`, `opsins_7tm_report.csv` | outputs of the run below |
| `work_7tm/` | intermediates (de-gapped input, `hmmalign` Stockholm) |

## Running it

```bash
conda activate structural

python filter_7tm.py <input.fasta> -o opsins_7tm \
    --profile animalOpsins_090426.hmm \
    --dtm-dir dtm_out
```

`--max-tm 7` and `--max-nterm 200` are the defaults, so the quality exclusions
below are on unless you pass `0` to disable them.

Outputs `opsins_7tm.faa`, `opsins_7tm_review.faa`, `opsins_7tm_dropped.faa`,
and `opsins_7tm_report.csv`. Omit `--dtm-dir` for a homology-only pass — note
that disables both quality filters, since each needs the topology prediction.

To rebuild the profile after editing the baits:

```bash
mafft --ep 0 --genafpair --maxiterate 1000 baits_NEW.fasta > baits_NEW.aln   # E-INS-i
python build_hmm.py baits_NEW.aln -o animalOpsins_NEW.hmm
python seed_qc.py baits_NEW.fasta -o seed_qc.csv
```

Then **re-run the held-out validation** before trusting it (see below) — a seed
change can quietly remove the ctenophore signal that makes TM7 work.

## Results

On `combined_081426_nocontam.fasta` (2,451 sequences):

| verdict | n | file |
|---|---|---|
| keep | 1,537 | `opsins_7tm.faa` |
| review (criteria disagree) | 70 | `opsins_7tm_review.faa` |
| drop | 844 | `opsins_7tm_dropped.faa` |
| no_prediction | 0 | — |

There is one removal verdict, `drop`, and every dropped sequence carries a
`drop_reason`:

| drop_reason | n |
|---|---|
| `<7_TM` | 829 |
| `>200aa_before_TM1` | 7 |
| `>7_TM` | 6 |
| `>7_TM+>200aa_before_TM1` | 1 |
| `<7_TM+>200aa_before_TM1` | 1 |

Reasons are combined with `+` when a sequence trips more than one rule, so the
per-criterion totals are `<7_TM` 830, `>7_TM` 7, `>200aa_before_TM1` 9. The
three FASTA files partition the input exactly: 1,537 + 70 + 844 = 2,451.

`no_prediction` is empty because `dtm_out/` covers all 2,451 sequences.

Where the kept count comes from:

| | kept |
|---|---|
| PF00001, no exclusions (first release) | 1,310, plus 351 unscored |
| `animalOpsins_090426`, no exclusions | 1,550 |
| `animalOpsins_090426` + exclusions | **1,537** |

So the profile change accounts for +240 kept sequences and the closing of the
351-sequence prediction gap; the quality exclusions remove 13 of those.

## Why sequences are dropped

Three rules remove a sequence. All are recorded in `drop_reason`, joined by `+`
when more than one applies.

**`<7_TM`** — both criteria call the 7TM core incomplete: profile coverage is
below `--hmm-cov` on at least one window *and* no predicted helix sits at that
position. Met by **830** sequences. Where only one criterion fails the verdict
is `review`, not `drop`.

**`>7_TM`** (`--max-tm 7`) — more than seven predicted TM helices. Met by **7**
sequences.

**`>200aa_before_TM1`** (`--max-nterm 200`) — more than 200 residues before the
TM1 window. Met by **9** sequences. This is the fusion / mis-assembly filter.
Bovine rhodopsin has 36 residues before TM1; across the kept set the median is
26 and the p99 is 125. The worst case here, `Aiptasia_pallida_1847`, is 2,836 aa
with **2,471** before TM1 — a complete 7TM core bolted onto unrelated sequence,
with a predicted signal peptide at 1–28 and a **proline** where the
retinal-binding lysine belongs. Six of the nine share that signature (proline at
K296 plus a several-hundred-residue N-terminal region); two carry the lysine
(`Hydra_vulgaris_XP_012559098.1`, `Acyrthosiphon_pisum_J9JXD7_canonical_r`) and
one has a gap there (`Tetranychusurticae_T1K299_canonical_r`).

The last two rules run independently of the 7TM verdict, so a sequence with
extra membrane segments or a large N-terminal region is dropped regardless of
how its core scores. That is why `>7_TM` and `>200aa_before_TM1` appear both
alone and combined with `<7_TM`.

The per-criterion totals above count sequences **meeting** each rule, so they
overlap and do not sum to 844. The `drop_reason` column partitions the 844
exclusively — see the table under Results.

`>7_TM` and `>200aa_before_TM1` are nearly disjoint: they overlap on one
sequence (`Hydra_vulgaris_XP_012559098.1`) and each misses what the other finds.
`CTE_Cydi_redx|D497-D8|17442c4g6i2.1` is the illustrative case: 1,104 aa, 747
residues before TM1, proline at K296 — and exactly 7 predicted helices, so only
the N-terminal rule sees it.

Note on the helix count: it is DeepTMHMM's raw output. Where DeepTMHMM splits a
single helix into two adjacent calls, the sequence really has seven TMs but is
reported as eight and is dropped. `Hydra_vulgaris_XP_012559098.1` (calls at
231–241 and 243–254, one residue apart) and
`Pundamilia_nyererei_XP_005748590_canonical_c` are the two cases here. This is
a deliberate choice: dropping them keeps the retained set unambiguous.

A `>7 TM helices` count filter is **not** a fusion filter. The fusion signature
is an extra helix *interior* to the 7TM core, reported in
`extra_helix_location`, and it currently matches **zero** sequences — all extras
sit before TM1 or after TM7. If you later want to filter on architecture rather
than count, that column is the one to use.

`res_at_K296` and `has_signal_peptide` are reported but **not** filtered on:
the input deliberately carries non-opsin GPCRs as outgroups (514 of the kept
sequences lack the retinal lysine), so either as a drop rule would remove them.
`has_signal_peptide` in particular is unreliable alone — 18 kept sequences with
an SP call have entirely normal N-termini.

## The two criteria

**A. Profile-HMM coverage — supplies register.** Each of rhodopsin's seven TM
spans is mapped onto the model; every window must be ≥70% covered by real
residues. Coverage counts **match columns only**: insert columns are not
comparable across sequences, because `hmmalign` places each sequence's inserted
residues independently. Match columns are read from HMMER's case convention
(uppercase = match), so the script works with or without an `RF` line — stock
Pfam models have `RF yes`, a `hmmbuild`/`pyhmmer` model from a plain alignment
has `RF no`.

**B. DeepTMHMM topology — supplies presence.** DeepTMHMM labels every residue,
so it yields helix *coordinates*, not just a count. Default `--dtm-mode position`
requires a predicted helix overlapping each expected window. `--dtm-mode count`
reproduces the older total-count rule, but that rule conflates "has a TM7" with
"is complete": a sequence missing TM1 has 6 helices *and* an intact TM7.

A sequence is kept only when both criteria pass. Disagreements go to `review`.

## Why the profile was changed

PF00001's TM7 columns carry high information content at the class-A NPxxY motif,
which most ctenophore opsins lack. `AFK83788.1` (*Mnemiopsis* opsin 1) reads
**WIIIY** at those five columns. `hmmalign` then scores delete+insert above
match, dumps the real TM7 into an insert, and `cov_TM7` collapses to ~0 — even
though that helix aligns gaplessly (22/22) to rhodopsin TM7 pairwise and
DeepTMHMM predicts seven helices. `hmmsearch` and `hmmalign` independently stop
at model state 229 of 260, and TM7 is states 241–260, so this is the optimal path
under that model rather than a filtering artifact. PF00001 also begins near the
end of TM1, contributing only 7 match columns to the TM1 window.

Measured on **258 ctenophore sequences held out of the seed set**:

| | PF00001 | this profile |
|---|---|---|
| TM7 covered | 65.5% | **88.0%** |
| all seven covered | 58.1% | **79.5%** |
| TM1 match columns | 7 | **25** |
| under-capture (helix at TM7, profile misses it) | 26 | **0** |
| over-capture (no helix at TM7, profile covers it) | 0 | **0** |

TM1–TM6 coverage is unchanged between the two profiles, so the gain is
TM7-specific rather than general permissiveness. All 17 non-opsin GPCR outgroups
score 100% on all seven windows under both profiles — `hmmalign` is glocal, so a
distant class-A receptor still traverses an opsin-specific model.

An earlier version carried a `--rescue-tm7` flag that kept sequences whose only
failing window was TM7. **It has been retired.** A rescue rule lets a broken
measurement stand and then overrides it; the right profile fixes the measurement.

## Validation

**Held-out design.** `K9LK83` in the bait set is effectively the same sequence as
`AFK83788.1`, so seeds must be excluded before measuring anything. All figures
above are on 258 ctenophore sequences absent from both the current and the
previous seed set.

**Truncation control.** On 205 sequences passing all seven windows *and* carrying
a predicted helix at all seven positions, TM7 was deleted and the sequences
rescored:

| removal mode | false-pass rate |
|---|---|
| C-terminal truncation at TM7 start | 0/205 (0.0%) |
| TM7 excised, C-terminal tail retained | 2/205 (1.0%) |

Both residual cases are the same protein entered twice — *Mnemiopsis* opsin 3
(`AFK83790`, Schnitzler et al. 2012; Feuda et al. 2014). Its intact `cov_TM7` is
0.86 and rises to 1.00 after excision, because a downstream sequence patch slides
into the vacated columns. That protein is also the documented case of a **local
K296 mimic**: a lysine sits at rhodopsin's K296 column in sequence alignment
without being the structurally conserved retinal-binding residue. The duplicate
entry is worth removing from the input.

Caveat: removing TM7 also perturbed TM6 alignment in ~17% of cases, so the
"passes all seven" figures from the control are conservative. The TM7-specific
false-pass rate is the clean number.

## Caveats

- **The residual 12% at TM7 is the data, not the profile.** Independently of any
  profile, DeepTMHMM finds ≥7 helices in only 80.6% of these held-out sequences.
  The profile covers TM7 in 88.0% — *above* that ceiling. Interior helices score
  96–99% because a transcript that lost its termini still contains its middle, so
  96–99% is not a standard TM7 should be held to. Of the 31 residual TM7 failures,
  7 fail TM7 alone with 0–24 residues left after TM6 (rhodopsin needs 74), and 24
  fail TM6 as well with TM1–TM5 at 1.00; DeepTMHMM reports 5–6 helices for all 31.
- **`seed_qc.py` K296 screening is necessary, not sufficient.** A local mimic can
  place a lysine at that column (see `AFK83790` above), so a positive is a screen,
  not a determination.
- **Single-position claims need a better alignment.** The apparent loss of the
  C110/C187 disulfide in ctenophore opsins is consistent across four seeds in
  three genera, but `mafft --add` can misalign locally; confirm against a
  structural or curated alignment before publishing it.
- **Check the bait FASTA ends with a newline.** Without one, concatenating it
  (`cat baits.fasta ref.faa`) silently glues the next header onto the last
  sequence, and the aligner then fails with no obvious cause. `baits_090426.fasta`
  is fine; `seed_qc.py` guards against it, ad-hoc `cat` does not.
- **Re-run the held-out validation after any seed edit.** Changing a single seed
  alters the whole E-INS-i alignment and therefore which columns become match
  states, which can shift TM1 and TM7 placement. The profile is stable for
  sequences carrying the retinal lysine and less so where that column is gapped,
  so a seed change can flip divergent sequences without touching the ones you
  care about — or the reverse.
- The profile is opsin-specific. Non-opsin GPCRs still score well (above), but
  their coverage values are less meaningful than an opsin's.
