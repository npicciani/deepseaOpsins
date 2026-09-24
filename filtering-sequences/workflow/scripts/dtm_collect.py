#!/usr/bin/env python3
"""Reattach headers to stored predictions and write one dataset's topology files.

The store is keyed by sequence, so this is where a dataset's own headers come
back: for every record in the input FASTA, look up the prediction by sequence
hash and emit it under that record's name. Two headers carrying the same
sequence each get a block referring to the same prediction, which is what keeps
filter_7tm.py scoring every input record.

Output is a single dtm_out/TMRs.gff3 rather than the batch_XXX/ tree the
earlier version produced. filter_7tm.py globs `**/TMRs.gff3` under --dtm-dir,
so it reads one file or thirteen without caring which.

Sequences skipped by the plan for being over --max-len are legitimately absent
and are simply omitted, which is what lands them in `no_prediction`. An
ELIGIBLE sequence missing from the store is a bug -- the batches should have
produced it -- so it is an error rather than a silent gap in the input to a
filter whose whole job is deciding what is complete.
"""
import argparse
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
    ap.add_argument("--out-gff3", required=True)
    ap.add_argument("--out-3line", required=True)
    ap.add_argument("--max-len", type=int, default=4000)
    a = ap.parse_args()

    recs = read_fasta(a.fasta)

    names = [n.split()[0] for n, _ in recs]
    dupes = sorted({n for n in names if names.count(n) > 1})
    if dupes:
        sys.exit(f"{a.fasta}: {len(dupes)} duplicated header token(s): "
                 f"{', '.join(dupes[:5])}{' ...' if len(dupes) > 5 else ''}\n"
                 f"filter_7tm.py joins topology to sequences by this token, so "
                 f"duplicates make that join ambiguous and would put two rows "
                 f"under one name in the report. Deduplicate the input.")

    out, seqs_by_name, skipped, missing = [], {}, [], []
    for name, seq in recs:
        key_name = name.split()[0]
        if not seq or len(seq) > a.max_len:
            skipped.append(key_name)
            continue
        key = dtm_store.seq_key(seq)
        if not dtm_store.has(a.cache_dir, key):
            missing.append(key_name)
            continue
        out.append((key_name, dtm_store.load(a.cache_dir, key)))
        seqs_by_name[key_name] = seq

    if missing:
        sys.exit(f"{len(missing)} eligible sequence(s) absent from the store "
                 f"after the batch jobs ran: "
                 f"{', '.join(missing[:5])}{' ...' if len(missing) > 5 else ''}")

    dtm_store.emit_gff3(out, a.out_gff3)
    dtm_store.emit_3line(out, a.out_3line, seqs_by_name)
    print(f"{len(out)} predictions written to {a.out_gff3}"
          + (f"; {len(skipped)} sequence(s) skipped as empty or over "
             f"{a.max_len} aa -> no_prediction" if skipped else ""))


if __name__ == "__main__":
    main()
