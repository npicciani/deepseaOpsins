#!/usr/bin/env python3
"""Apply the seven-TM filter's thresholds to a scoring report.

The report (../scoring-sequences, score_7tm.py) holds measurements only. This
script makes every decision, so the thresholds can be changed and re-applied in
seconds without re-running hmmalign or DeepTMHMM. It reads nothing but the
report, the matching sequences FASTA, and optionally the report's metadata.

Two criteria, each decided from report columns:

  A) profile coverage   pass if min_cov >= --hmm-cov  (every TM window covered)
  B) DeepTMHMM topology
       --dtm-mode position  (default) pass if helix_pattern == 1111111: a
                            predicted helix overlaps each of the seven expected
                            positions. The like-for-like comparison against A.
       --dtm-mode count     pass if n_tm_helices >= --min-tm. Cheaper, but it
                            conflates "has a TM7" with "is complete": a sequence
                            missing TM1 has 6 helices AND an intact TM7, and a
                            count rule cannot tell those apart.
     A sequence with has_topology = 0 has no B at all -- a data gap, never a
     failure.

Every reason a sequence is removed is recorded in drop_reason, joined by '+'
when more than one applies, and all removals share the verdict "drop":

  <7_TM              both criteria call the 7TM core incomplete
  >N_TM              DeepTMHMM's raw helix count exceeds --max-tm (0 disables).
                     Where DeepTMHMM splits one helix into two adjacent calls
                     the sequence has seven but is counted as eight, and drops.
  >Naa_before_TM1    more than --max-nterm residues precede TM1 (0 disables) --
                     the fusion / mis-assembly signature. Bovine rhodopsin has
                     36; the p99 of a real opsin set is ~125.

Verdicts, in order of precedence:

  no_prediction   has_topology = 0 (drop_reason is still filled in, so a
                  sequence that would have dropped anyway is visible)
  drop            any drop reason
  keep            both criteria pass
  review          the criteria disagree

Outputs, all always written (an empty FASTA means that class is empty):

  <prefix>.faa, <prefix>_review.faa, <prefix>_dropped.faa,
  <prefix>_no_prediction.faa
  <prefix>_verdicts.csv    the report plus pass_hmm, pass_deeptmhmm,
                           drop_reason, verdict
  <prefix>_params.json     the thresholds applied, and the report's metadata
"""
import argparse
import csv
import json
import os
import sys

REQUIRED = ("name", "min_cov", "has_topology", "n_tm_helices",
            "helix_pattern", "n_before_TM1")


def read_fasta(path):
    recs, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def write_fasta(path, recs):
    with open(path, "w") as fh:
        for n, s in recs:
            fh.write(f">{n}\n{s}\n")


def opt_int(v):
    return None if v == "" else int(v)


def decide(row, a):
    """Return (pass_hmm, pass_dtm, reasons, verdict) for one report row."""
    pass_hmm = float(row["min_cov"]) >= a.hmm_cov

    has_topo = row["has_topology"] == "1"
    ntm = opt_int(row["n_tm_helices"]) if has_topo else None
    if not has_topo:
        pass_dtm = None
    elif a.dtm_mode == "count":
        pass_dtm = ntm >= a.min_tm
    else:
        pass_dtm = row["helix_pattern"] == "1" * 7

    n_before = opt_int(row["n_before_TM1"])
    reasons = []
    if pass_dtm is not None and not pass_hmm and not pass_dtm:
        reasons.append("<7_TM")
    if a.max_tm and ntm is not None and ntm > a.max_tm:
        reasons.append(f">{a.max_tm}_TM")
    if a.max_nterm and n_before is not None and n_before > a.max_nterm:
        reasons.append(f">{a.max_nterm}aa_before_TM1")

    if pass_dtm is None:
        verdict = "no_prediction"
    elif reasons:
        verdict = "drop"
    elif pass_hmm and pass_dtm:
        verdict = "keep"
    else:
        verdict = "review"
    return pass_hmm, pass_dtm, reasons, verdict


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("report", help="report.csv from scoring-sequences")
    ap.add_argument("sequences", help="sequences.faa from scoring-sequences")
    ap.add_argument("-o", "--outprefix", required=True)
    ap.add_argument("--meta", default=None, help="report.json, carried into params")
    ap.add_argument("--hmm-cov", type=float, default=0.7)
    ap.add_argument("--dtm-mode", choices=("position", "count"), default="position")
    ap.add_argument("--min-tm", type=int, default=7)
    ap.add_argument("--max-tm", type=int, default=7)
    ap.add_argument("--max-nterm", type=int, default=200)
    a = ap.parse_args()

    with open(a.report, newline="") as fh:
        reader = csv.DictReader(fh)
        missing = [c for c in REQUIRED if c not in (reader.fieldnames or [])]
        if missing:
            sys.exit(f"{a.report} lacks column(s) {', '.join(missing)} -- not a "
                     f"scoring-sequences report, or one from an older version")
        fields = list(reader.fieldnames)
        rows = list(reader)
    if not rows:
        sys.exit(f"{a.report} has no rows")

    seqs = read_fasta(a.sequences)
    by_name = dict(seqs)
    if len(by_name) != len(seqs):
        sys.exit(f"{a.sequences} repeats a header; verdicts could not be joined "
                 f"to sequences unambiguously")
    names = [r["name"] for r in rows]
    if len(set(names)) != len(names) or set(names) != set(by_name):
        only_r = sorted(set(names) - set(by_name))[:5]
        only_s = sorted(set(by_name) - set(names))[:5]
        sys.exit(f"report and sequences describe different sets "
                 f"(report only: {only_r}; FASTA only: {only_s}). They must come "
                 f"from the same scoring run.")

    out = {"keep": [], "review": [], "drop": [], "no_prediction": []}
    for r in rows:
        pass_hmm, pass_dtm, reasons, verdict = decide(r, a)
        r["pass_hmm"] = int(pass_hmm)
        r["pass_deeptmhmm"] = "" if pass_dtm is None else int(pass_dtm)
        r["drop_reason"] = "+".join(reasons)
        r["verdict"] = verdict
        out[verdict].append((r["name"], by_name[r["name"]]))

    os.makedirs(os.path.dirname(os.path.abspath(a.outprefix)), exist_ok=True)
    write_fasta(f"{a.outprefix}.faa", out["keep"])
    write_fasta(f"{a.outprefix}_review.faa", out["review"])
    write_fasta(f"{a.outprefix}_dropped.faa", out["drop"])
    write_fasta(f"{a.outprefix}_no_prediction.faa", out["no_prediction"])
    with open(f"{a.outprefix}_verdicts.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields + ["pass_hmm", "pass_deeptmhmm",
                                                    "drop_reason", "verdict"],
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    params = {
        "thresholds": {"hmm_cov": a.hmm_cov, "dtm_mode": a.dtm_mode,
                       "min_tm": a.min_tm, "max_tm": a.max_tm,
                       "max_nterm": a.max_nterm},
        "report": a.report,
        "counts": {k: len(v) for k, v in out.items()},
    }
    if a.meta:
        with open(a.meta) as fh:
            params["scoring"] = json.load(fh)
    with open(f"{a.outprefix}_params.json", "w") as fh:
        json.dump(params, fh, indent=2)
        fh.write("\n")

    print(f"{len(rows)} scored: " + ", ".join(f"{k}={len(v)}" for k, v in out.items()))


if __name__ == "__main__":
    main()
