#!/usr/bin/env python3
"""Condense opsins_7tm_report.csv into the tables the README reports.

Three blocks, long-format, one TSV:

  verdict       keep / review / drop / no_prediction        counts partition the input
  drop_reason   the exclusive `+`-joined reason strings     counts partition the drops
  criterion     sequences MEETING each rule                 counts OVERLAP, by design

The last block is the one that is easy to misread: a sequence tripping both
`<7_TM` and `>200aa_before_TM1` is counted under each, so criterion counts sum
to more than the number of drops. The drop_reason block is the exclusive one.
"""
import argparse
import csv
from collections import Counter


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("report", help="opsins_7tm_report.csv")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--sample", default="")
    a = ap.parse_args()

    with open(a.report, newline="") as fh:
        rows = list(csv.DictReader(fh))

    verdicts = Counter(r["verdict"] for r in rows)
    reasons = Counter(r["drop_reason"] for r in rows if r["verdict"] == "drop")
    criteria = Counter()
    for r in rows:
        for part in filter(None, r["drop_reason"].split("+")):
            criteria[part] += 1

    with open(a.out, "w", newline="") as fh:
        # csv defaults to CRLF; a TSV read by awk/pandas/R should be LF
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample", "block", "key", "n"])
        w.writerow([a.sample, "total", "sequences_scored", len(rows)])
        for k in ("keep", "review", "drop", "no_prediction"):
            w.writerow([a.sample, "verdict", k, verdicts.get(k, 0)])
        for k, n in sorted(reasons.items(), key=lambda kv: -kv[1]):
            w.writerow([a.sample, "drop_reason", k, n])
        for k, n in sorted(criteria.items(), key=lambda kv: -kv[1]):
            w.writerow([a.sample, "criterion_met", k, n])

    print(f"{a.sample or a.report}: {len(rows)} scored, "
          + ", ".join(f"{k}={verdicts.get(k, 0)}"
                      for k in ("keep", "review", "drop", "no_prediction")))


if __name__ == "__main__":
    main()
