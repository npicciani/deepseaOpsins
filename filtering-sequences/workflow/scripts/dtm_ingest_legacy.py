#!/usr/bin/env python3
"""Load DeepTMHMM output that already exists on disk into the store.

One-time migration, run by hand, not part of the DAG. Every batch directory
produced before the store existed -- the 13 batches covering the 2,451-sequence
curated set, and anything run since -- represents cloud compute that has
already been paid for. This reads those directories and writes per-sequence
records, so the first workflow run after the migration plans zero DeepTMHMM runs.

No FASTA is needed: predicted_topologies.3line carries the sequences, so a
batch directory is self-describing.

    python dtm_ingest_legacy.py results/devivo2023/dtm_out ../curated/work/dtm_out \\
        --cache-dir resources/dtm_cache

CONFLICTS are the reason this prints a summary rather than running silently. If
a sequence is already stored with a DIFFERENT prediction, the two runs
disagree -- which on a cloud service usually means the hosted model changed
between them. That is exactly the provenance gap a version number would close,
and it is worth knowing about before the two sets of predictions are mixed in
one analysis. Existing records are never overwritten; conflicts are reported
and the incoming record is dropped.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dtm_store  # noqa: E402

PAIR = ("TMRs.gff3", "predicted_topologies.3line")


def find_batch_dirs(roots):
    found = []
    for root in roots:
        if all(os.path.exists(os.path.join(root, f)) for f in PAIR):
            found.append(root)
            continue
        for dirpath, _, filenames in os.walk(root):
            if all(f in filenames for f in PAIR):
                found.append(dirpath)
    return sorted(set(found))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dirs", nargs="+",
                    help="batch directories, or parents to walk")
    ap.add_argument("--cache-dir", required=True)
    ap.add_argument("-n", "--dry-run", action="store_true")
    ap.add_argument("--source", default="biolib:DTU/DeepTMHMM",
                    help="provenance tag for the imported records (default "
                         "assumes they came from the BioLib web service)")
    a = ap.parse_args()

    batch_dirs = find_batch_dirs(a.dirs)
    if not batch_dirs:
        sys.exit(f"no directory containing both {' and '.join(PAIR)} found "
                 f"under: {', '.join(a.dirs)}")

    os.makedirs(a.cache_dir, exist_ok=True)
    n_new = n_same = 0
    conflicts = []

    for bdir in batch_dirs:
        gff3 = os.path.join(bdir, PAIR[0])
        line3 = os.path.join(bdir, PAIR[1])
        label = os.path.relpath(bdir)
        seqs = dtm_store.parse_3line(line3)
        feats = dtm_store.parse_gff3(gff3)
        fresh = 0
        for name, (seq, topo) in seqs.items():
            if name not in feats:
                raise dtm_store.FormatError(
                    f"{label}: {name} is in the 3line file but not the gff3")
            key = dtm_store.seq_key(seq)
            incoming = [list(f) for f in feats[name]["features"]]
            if dtm_store.has(a.cache_dir, key):
                existing = dtm_store.load(a.cache_dir, key)
                if existing["features"] == incoming \
                        and existing["topology"] == topo:
                    n_same += 1
                else:
                    conflicts.append((name, existing.get("batch", "?"), label))
                continue
            fresh += 1
            if not a.dry_run:
                dtm_store.save(a.cache_dir, {
                    "sha1": key,
                    "length": len(seq),
                    "n_tm_helices": sum(1 for t, _, _ in feats[name]["features"]
                                        if t == "TMhelix"),
                    "comments": feats[name]["comments"],
                    "features": incoming,
                    "topology": topo,
                    "source": a.source,
                    "batch": f"legacy:{label}",
                    "ingested": "migrated",
                })
        n_new += fresh
        print(f"{label}: {len(seqs)} records, {fresh} new")

    print(f"\n{len(batch_dirs)} batch director{'y' if len(batch_dirs) == 1 else 'ies'}: "
          f"{n_new} new record(s), {n_same} already stored and identical, "
          f"{len(conflicts)} conflict(s)")
    if a.dry_run:
        print("dry run -- nothing written")
    if conflicts:
        print("\nCONFLICT: same sequence, different prediction. The hosted "
              "model may have changed between these runs; do not mix them "
              "without checking.")
        for name, was, now in conflicts[:20]:
            print(f"  {name}: stored from {was}, differs in {now}")
        if len(conflicts) > 20:
            print(f"  ... and {len(conflicts) - 20} more")
        sys.exit(1)


if __name__ == "__main__":
    main()
