#!/usr/bin/env python3
"""Decide what actually needs predicting, and cut it into batches.

This is the checkpoint that makes the store worth having. It hashes every
eligible sequence in the input, asks the store which are already predicted, and
writes batch FASTAs for ONLY the remainder. Add five sequences to a dataset and
this plans one five-sequence batch rather than re-running the full
batch they happened to land in.

The work set is knowable only at run time -- it depends on the contents of the
store -- which is why this is a Snakemake `checkpoint` and not a plain rule.

Two kinds of sequence never reach a batch:

  * already in the store (the point of the exercise);
  * longer than --max-len, which DeepTMHMM is slow on. These are skipped, not
    truncated: a truncated prediction would be silently wrong. They arrive at
    filter_7tm.py with no topology record and are classed `no_prediction`.

Duplicate sequences under different headers collapse to one batch entry, since
they share a hash. They are re-expanded per header at collect time.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dtm_store  # noqa: E402
from fastaio import read_fasta  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta")
    ap.add_argument("--cache-dir", required=True)
    ap.add_argument("--plan-dir", required=True)
    ap.add_argument("-b", "--batch-size", type=int, default=200)
    ap.add_argument("--max-len", type=int, default=4000)
    a = ap.parse_args()

    os.makedirs(a.plan_dir, exist_ok=True)
    os.makedirs(a.cache_dir, exist_ok=True)

    recs = read_fasta(a.fasta)
    eligible, oversize, empty = [], [], 0
    for name, seq in recs:
        if not seq:
            empty += 1
        elif len(seq) > a.max_len:
            oversize.append(name.split()[0])
        else:
            eligible.append((name.split()[0], seq))

    # dedupe on the hash, keeping first occurrence, then drop what's cached
    missing, seen = [], set()
    n_cached = 0
    for name, seq in eligible:
        key = dtm_store.seq_key(seq)
        if key in seen:
            continue
        seen.add(key)
        if dtm_store.has(a.cache_dir, key):
            n_cached += 1
        else:
            missing.append((key, name, seq))

    # a stale plan from an earlier run would be silently re-executed
    for old in os.listdir(a.plan_dir):
        if old.startswith("batch_") and old.endswith(".faa"):
            os.remove(os.path.join(a.plan_dir, old))

    batches = [missing[i:i + a.batch_size]
               for i in range(0, len(missing), a.batch_size)]
    for bi, batch in enumerate(batches):
        with open(os.path.join(a.plan_dir, f"batch_{bi:03d}.faa"), "w") as fh:
            for _, name, seq in batch:
                fh.write(f">{name}\n{seq}\n")

    plan = {
        "input": a.fasta,
        "records": len(recs),
        "empty": empty,
        "oversize": oversize,
        "unique_eligible": len(seen),
        "already_cached": n_cached,
        "to_predict": len(missing),
        "batches": len(batches),
        "batch_size": a.batch_size,
        "max_len": a.max_len,
    }
    with open(os.path.join(a.plan_dir, "plan.json"), "w") as fh:
        json.dump(plan, fh, indent=2, sort_keys=True)
        fh.write("\n")

    print(f"{len(recs)} records -> {len(seen)} unique eligible sequences: "
          f"{n_cached} cached, {len(missing)} to predict in {len(batches)} "
          f"batch(es)")
    if oversize:
        print(f"skipped {len(oversize)} over {a.max_len} aa: "
              f"{', '.join(oversize[:5])}{' ...' if len(oversize) > 5 else ''}")
    if empty:
        print(f"skipped {empty} empty record(s)")


if __name__ == "__main__":
    main()
