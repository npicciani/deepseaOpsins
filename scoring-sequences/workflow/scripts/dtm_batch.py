#!/usr/bin/env python3
"""Run DeepTMHMM over one planned batch, then split the result into
per-sequence records in the store.

DeepTMHMM runs from a local install -- the academic-licensed release, not the
BioLib web service -- through its launcher (`deeptmhmm <fasta> <outdir>` by
default). The batch comes from dtm_plan.py and holds only sequences the store
lacks, so this never recomputes cached work.

Each invocation loads ESM-1b (~2.6 GB) before predicting, so batches should
be large: the load is paid once per batch, not per sequence.

DeepTMHMM also writes one ESM-1b embedding per sequence under
<outdir>/embeddings/ (~2.5 MB each). Nothing downstream reads them, so they
are deleted once the batch is stored.

Snakemake cannot declare the store files as outputs -- which sequences a batch
holds is only known after the checkpoint -- so the declared output is a
sentinel written last, after every record has landed. A job killed midway
leaves complete records for whatever it finished (each is written atomically)
and no sentinel; the next plan re-runs only the rest.
"""
import argparse
import contextlib
import fcntl
import os
import shutil
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dtm_store  # noqa: E402

EXPECTED = ("TMRs.gff3", "predicted_topologies.3line")


@contextlib.contextmanager
def claim_gpu(gpus, lock_dir):
    """Hold an exclusive lock on one free GPU for the duration of the batch.

    Snakemake's `--resources gpu=N` caps how many batch jobs run at once but
    cannot say WHICH device each gets, so two concurrent jobs could land on the
    same card. Each job instead takes a non-blocking flock on
    <lock_dir>/gpu<i>.lock for the first free device and holds it until done.
    The kernel releases a flock when its process dies, so a killed job never
    leaves a device stuck.
    """
    if not gpus:
        yield None
        return
    os.makedirs(lock_dir, exist_ok=True)
    while True:
        for g in gpus:
            fh = open(os.path.join(lock_dir, f"gpu{g}.lock"), "w")
            try:
                fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                fh.close()
                continue
            try:
                yield g
            finally:
                fcntl.flock(fh, fcntl.LOCK_UN)
                fh.close()
            return
        time.sleep(5)  # more jobs than devices: wait for one to free


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("batch_fasta", help="a batch_XXX.faa written by dtm_plan.py")
    ap.add_argument("--cache-dir", required=True)
    ap.add_argument("--run-dir", required=True,
                    help="DeepTMHMM's output directory for this batch, kept as a log")
    ap.add_argument("--sentinel", required=True)
    ap.add_argument("--deeptmhmm", default="deeptmhmm",
                    help="the DeepTMHMM launcher (default: `deeptmhmm` on PATH)")
    ap.add_argument("--gpus", default="",
                    help="comma-separated CUDA device indices this batch may "
                         "claim, e.g. '0,1'; empty = the launcher's default")
    ap.add_argument("--lock-dir", default=None,
                    help="where GPU claim locks live (default: beside --run-dir)")
    ap.add_argument("--source", default=dtm_store.SOURCE,
                    help="provenance tag written into every record")
    a = ap.parse_args()

    n_in = sum(1 for line in open(a.batch_fasta) if line.startswith(">"))
    if n_in == 0:
        sys.exit(f"{a.batch_fasta} contains no sequences")
    launcher = shutil.which(a.deeptmhmm)
    if launcher is None:
        sys.exit(f"DeepTMHMM launcher {a.deeptmhmm!r} not found. It is licensed "
                 f"software installed per host; set deeptmhmm.bin in "
                 f"config/config.yaml to its path.")

    # a stale partial run would otherwise be read as this batch's output
    if os.path.isdir(a.run_dir):
        shutil.rmtree(a.run_dir)
    os.makedirs(os.path.dirname(os.path.abspath(a.run_dir)), exist_ok=True)

    gpus = [g.strip() for g in a.gpus.split(",") if g.strip()]
    lock_dir = a.lock_dir or os.path.join(
        os.path.dirname(os.path.abspath(a.run_dir)), ".gpu_locks")
    label = os.path.basename(a.run_dir.rstrip("/"))

    with claim_gpu(gpus, lock_dir) as gpu:
        env = dict(os.environ)
        if gpu is not None:
            env["CUDA_VISIBLE_DEVICES"] = gpu
        where = f"GPU {gpu}" if gpu is not None else "default device"
        print(f"{label}: {n_in} sequences on {where}", flush=True)
        proc = subprocess.run([launcher, a.batch_fasta, a.run_dir], env=env)
    if proc.returncode != 0:
        sys.exit(f"{label}: DeepTMHMM exited {proc.returncode}")

    paths = {f: os.path.join(a.run_dir, f) for f in EXPECTED}
    missing = [f for f, p in paths.items() if not os.path.exists(p)]
    if missing:
        sys.exit(f"{label}: DeepTMHMM exited 0 but did not write "
                 f"{', '.join(missing)} to {a.run_dir}")

    keys = dtm_store.ingest_batch(paths["TMRs.gff3"],
                                  paths["predicted_topologies.3line"],
                                  a.cache_dir, label, source=a.source)
    if len(keys) != n_in:
        sys.exit(f"{label}: submitted {n_in} sequences but DeepTMHMM returned "
                 f"{len(keys)} -- refusing to write the sentinel")

    shutil.rmtree(os.path.join(a.run_dir, "embeddings"), ignore_errors=True)

    os.makedirs(os.path.dirname(os.path.abspath(a.sentinel)), exist_ok=True)
    with open(a.sentinel, "w") as fh:
        fh.write(f"{label}\t{len(keys)} sequences\t{a.source}\n")
    print(f"{label}: {len(keys)} records stored", flush=True)


if __name__ == "__main__":
    main()
