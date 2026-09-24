#!/usr/bin/env python3
"""Keep only opsin / GPCR protein sequences that span all seven transmembrane helices.

Two independent criteria are combined.

  A) Profile-HMM coverage (homology-aware, gives REGISTER).
     Sequences are aligned to a profile HMM with `hmmalign`, together with bovine
     rhodopsin (UniProt P02699) whose TM1-TM7 boundaries are annotated from
     structure. Rhodopsin's TM spans are mapped onto the model, and each of the
     seven windows must be covered by real residues above --hmm-cov.

     Coverage is computed over profile MATCH columns only. Insert columns are not
     comparable across sequences -- hmmalign places each sequence's inserted
     residues independently -- so counting them makes the score partly a measure
     of N-terminal extension length rather than of helix presence. Match columns
     are identified from HMMER's case convention (uppercase = match state,
     lowercase = insert), which works whether or not the profile carries an
     `RF` reference-annotation line. Stock Pfam models have `RF yes`; a model
     built by `hmmbuild`/`pyhmmer` from a plain alignment has `RF no`.

  B) DeepTMHMM topology (homology-free, gives PRESENCE).
     DeepTMHMM predicts per-residue topology, so it yields both the NUMBER of
     transmembrane helices and their COORDINATES. Two modes:

       --dtm-mode position  (default) a predicted helix must overlap each of the
                            seven expected window positions. This is the
                            like-for-like comparison against criterion A.
       --dtm-mode count     total predicted helices >= --min-tm. Cheaper, but it
                            conflates "has a TM7" with "is complete": a sequence
                            missing TM1 has 6 helices AND an intact TM7, and a
                            count rule cannot tell those apart.

WHY THE PROFILE MATTERS (and why --rescue-tm7 is gone)
-----------------------------------------------------
Pfam 7tm_1 (PF00001) systematically under-calls TM7 in ctenophores. Its TM7
columns carry high information content at the class-A NPxxY motif, which most
ctenophore opsins lack (e.g. AFK83788.1 / Mnemiopsis opsin 1 reads WIIIY at
those five columns). For such a sequence hmmalign scores delete+insert above
match, so the real TM7 is dumped into an insert and cov_TM7 collapses -- even
though the helix aligns gaplessly to rhodopsin TM7 pairwise and DeepTMHMM sees
seven helices. PF00001 also begins near the end of TM1, contributing only 7
match columns to the TM1 window.

Measured on 258 ctenophore sequences held out of the seed set, with a
clade-inclusive opsin profile instead of PF00001:

    TM7 covered      : 65.5% -> 88.0%
    all seven covered: 58.1% -> 79.5%
    TM1 match columns:     7 -> 25
    under-capture (helix present at TM7, profile misses it): 26 -> 0
    over-capture  (no helix at TM7, profile covers it)     :  0 -> 0

An earlier version of this script carried a --rescue-tm7 flag that kept
sequences whose only failing window was TM7. A profile that represents the
clade removes the need for it, so the flag has been retired: a rescue rule
lets a broken measurement stand and then overrides it, whereas the right
profile fixes the measurement.

Truncation control (205 sequences passing all seven windows AND carrying a
helix at all seven positions): deleting TM7 gives a false-pass rate of 0/205
for C-terminal truncation and 2/205 for internal excision. Both residual cases
are the same protein entered twice, Mnemiopsis opsin 3 (AFK83790), whose
downstream tail can slide into the vacated TM7 columns.

Usage
-----
    python filter_7tm.py input.fasta -o opsins_7tm \\
        --profile animalOpsins.hmm --dtm-dir dtm_out

Omit --dtm-dir to run criterion A alone (fully local, no uploads).
Omit --profile to fetch and use stock PF00001 (NOT recommended for ctenophores).
"""
import argparse
import csv
import glob
import json
import os
import subprocess
import sys
import urllib.request

P02699_JSON = ("https://rest.uniprot.org/uniprotkb/P02699.json"
               "?fields=ft_transmem,sequence")
PF00001_HMM = ("https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/"
               "PF00001?annotation=hmm")


# ---------------------------------------------------------------- FASTA / IO

def read_fasta(path):
    recs, name, buf = [], None, []
    for line in open(path):
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


def degap(seq):
    return (seq.replace("-", "").replace(".", "").replace("*", "").upper())


# ---------------------------------------------------------------- reference

def fetch_reference(refdir, profile=None):
    """Return (hmm_path, rhodopsin_seq, [(tm_start, tm_end), ...])."""
    os.makedirs(refdir, exist_ok=True)

    if profile:
        hmm = profile
        if not os.path.exists(hmm):
            sys.exit(f"profile not found: {hmm}")
    else:
        hmm = os.path.join(refdir, "PF00001.hmm")
        if not os.path.exists(hmm):
            import gzip
            gz = hmm + ".gz"
            urllib.request.urlretrieve(PF00001_HMM, gz)
            with gzip.open(gz, "rb") as f_in, open(hmm, "wb") as f_out:
                f_out.write(f_in.read())

    js = os.path.join(refdir, "P02699.json")
    if not os.path.exists(js):
        urllib.request.urlretrieve(P02699_JSON, js)
    d = json.load(open(js))
    seq = d["sequence"]["value"]
    spans = [(f["location"]["start"]["value"], f["location"]["end"]["value"])
             for f in d["features"] if f["type"] == "Transmembrane"]
    if len(spans) != 7:
        sys.exit(f"expected 7 annotated TM spans in P02699, got {len(spans)}")
    return hmm, seq, spans


# ---------------------------------------------------------------- alignment

def parse_stockholm(path):
    """Return {name: aligned_string} from an hmmalign Stockholm file."""
    seqs = {}
    for line in open(path):
        line = line.rstrip("\n")
        if not line.strip() or line.startswith("#") or line.startswith("//"):
            continue
        parts = line.split(None, 1)
        if len(parts) == 2:
            seqs[parts[0]] = seqs.get(parts[0], "") + parts[1].strip()
    return seqs


def hmm_align(fasta, hmm, ref_seq, workdir, threads, refname="P02699_BOVIN_RHO"):
    os.makedirs(workdir, exist_ok=True)
    combined = os.path.join(workdir, "with_ref.faa")
    recs = read_fasta(fasta)
    write_fasta(combined, [(refname, ref_seq)] + recs)
    sto = os.path.join(workdir, "aln.sto")
    with open(sto, "w") as out:
        subprocess.run(["hmmalign", "--amino", "--outformat", "Stockholm",
                        hmm, combined], stdout=out, check=True)
    return parse_stockholm(sto), refname


def tm_windows(aln, refname, ref_spans):
    """Alignment columns of each rhodopsin TM span, restricted to MATCH states.

    HMMER writes match-state residues uppercase and insert-state residues
    lowercase, so the reference residue's case identifies the state. This is
    used instead of the `#=GC RF` line because a profile built from a plain
    alignment has no RF annotation.
    """
    ref = aln[refname]
    pos, res2col = 0, {}
    for i, ch in enumerate(ref):
        if ch.isalpha():
            pos += 1
            res2col[pos] = (i, ch.isupper())
    windows = []
    for k, (s, e) in enumerate(ref_spans, 1):
        cols = [c for p in range(s, e + 1) if p in res2col
                for c, is_match in [res2col[p]] if is_match]
        if not cols:
            sys.exit(f"TM{k} maps to zero match columns -- profile does not "
                     f"cover this helix; check the profile/reference pairing")
        windows.append(cols)
    return windows


def column_to_residue(seq):
    """{alignment column -> 1-based residue index} for one aligned sequence."""
    m, k = {}, 0
    for i, ch in enumerate(seq):
        if ch.isalpha():
            k += 1
            m[i] = k
    return m


# ---------------------------------------------------------------- DeepTMHMM

def parse_dtm(dtm_dir):
    """Return (counts, intervals, regions), or three Nones if no output present.

    counts[name]     -> number of predicted TM helices (0 is meaningful)
    intervals[name]  -> [(start, end), ...] 1-based residue coordinates
    regions[name]    -> [(type, start, end), ...] every predicted region,
                        including 'signal', 'inside' and 'outside'

    The per-sequence '# <name> Number of predicted TMRs: N' comment is the
    authoritative record of which sequences were processed, so a sequence with
    genuinely zero helices stays distinguishable from one never run.
    """
    files = sorted(glob.glob(os.path.join(dtm_dir, "**", "TMRs.gff3"),
                             recursive=True))
    if not files:
        return None, None, None
    counts, intervals, regions = {}, {}, {}
    for f in files:
        for line in open(f):
            if line.startswith("# ") and "Number of predicted TMRs" in line:
                nm = line[2:].rsplit(" Number of predicted TMRs:", 1)[0]
                key = nm.split()[0]
                counts.setdefault(key, 0)
                intervals.setdefault(key, [])
            elif line.startswith("#") or line.startswith("//") or not line.strip():
                continue
            else:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 4:
                    key = p[0].split()[0]
                    regions.setdefault(key, []).append(
                        (p[1], int(p[2]), int(p[3])))
                    if p[1] == "TMhelix":
                        counts[key] = counts.get(key, 0) + 1
                        intervals.setdefault(key, []).append(
                            (int(p[2]), int(p[3])))
    return counts, intervals, regions


def helix_at_windows(aligned_seq, windows, helices):
    """Per window, True if a predicted helix overlaps that window's residues."""
    c2r = column_to_residue(aligned_seq)
    out = []
    for cols in windows:
        rs = [c2r[c] for c in cols if c in c2r]
        if not rs:
            out.append(False)
            continue
        lo, hi = min(rs), max(rs)
        out.append(any(not (end < lo or start > hi) for start, end in helices))
    return out


def window_residue_ranges(aligned_seq, windows):
    """Per window, the (first, last) residue index of this sequence in it."""
    c2r = column_to_residue(aligned_seq)
    out = []
    for cols in windows:
        rs = [c2r[c] for c in cols if c in c2r]
        out.append((min(rs), max(rs)) if rs else None)
    return out


def extra_helix_locations(ranges, helices):
    """Where do predicted helices beyond the seven expected windows sit?

    'interior' is the fusion / mis-assembly signature: an extra membrane
    segment *inside* the 7TM core. 'before_TM1' / 'after_TM7' are terminal
    extensions or marginal hydrophobic calls. An empty list with a helix count
    above seven means one window absorbed two adjacent calls, i.e. DeepTMHMM
    split a single helix in two -- the sequence really has seven.
    """
    occ = [r for r in ranges if r]
    if not occ:
        return []
    lo7, hi7 = min(r[0] for r in occ), max(r[1] for r in occ)
    out = []
    for start, end in sorted(helices):
        if any(r and not (end < r[0] or start > r[1]) for r in ranges):
            continue
        out.append("before_TM1" if end < lo7
                   else "after_TM7" if start > hi7 else "interior")
    return out


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta", help="input protein FASTA (gaps are stripped)")
    ap.add_argument("-o", "--outprefix", default="opsins_7tm")
    ap.add_argument("--profile", default=None,
                    help="profile HMM to score against. Strongly recommended: a "
                         "clade-inclusive opsin profile. Default fetches stock "
                         "PF00001, which under-calls TM7 in ctenophores.")
    ap.add_argument("--dtm-dir", default=None,
                    help="directory of existing DeepTMHMM output (searched "
                         "recursively for TMRs.gff3). Omit to use criterion A alone.")
    ap.add_argument("--dtm-mode", choices=("position", "count"), default="position",
                    help="how criterion B is evaluated (default: position)")
    ap.add_argument("--hmm-cov", type=float, default=0.7,
                    help="min fraction of each TM window's match columns covered")
    ap.add_argument("--min-tm", type=int, default=7,
                    help="min predicted helices for --dtm-mode count")
    ap.add_argument("--max-tm", type=int, default=7,
                    help="DROP sequences with more than this many predicted TM "
                         "helices (default 7; drop_reason '>7_TM'). Note the count "
                         "is DeepTMHMM's raw output: where it splits one helix into "
                         "two adjacent calls the sequence really has seven but is "
                         "reported as eight, and it is dropped. Set 0 to disable. "
                         "Requires --dtm-dir.")
    ap.add_argument("--max-nterm", type=int, default=200,
                    help="DROP sequences with more than this many residues before "
                         "the TM1 window (default 200; drop_reason "
                         "'>200aa_before_TM1'). This is the fusion / mis-assembly "
                         "filter: bovine rhodopsin has 36, and the p99 of a real "
                         "opsin set is ~125. Set 0 to disable.")
    ap.add_argument("--refdir", default="ref")
    ap.add_argument("--workdir", default="work_7tm")
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()

    recs = read_fasta(a.fasta)
    clean = [(n, degap(s)) for n, s in recs]
    clean = [(n, s) for n, s in clean if s]
    if not clean:
        sys.exit("no sequences read")
    os.makedirs(a.workdir, exist_ok=True)
    clean_fa = os.path.join(a.workdir, "degapped.faa")
    write_fasta(clean_fa, clean)
    print(f"input     : {len(recs)} records -> {len(clean)} non-empty, gaps stripped")

    hmm, ref_seq, ref_spans = fetch_reference(a.refdir, a.profile)
    print(f"profile   : {hmm}")
    aln, refname = hmm_align(clean_fa, hmm, ref_seq, a.workdir, a.threads)
    windows = tm_windows(aln, refname, ref_spans)
    print(f"match columns per TM window: {[len(w) for w in windows]}")

    counts, intervals, dtm_regions = (parse_dtm(a.dtm_dir) if a.dtm_dir
                                      else (None, None, None))
    if counts is None:
        print("DeepTMHMM : none found -- criterion A only")
    else:
        print(f"DeepTMHMM : records for {len(counts)} sequences "
              f"(mode={a.dtm_mode})")
        print(f"drop rules: >{a.max_tm} TM helices" if a.max_tm else
              "drop rules: helix-count filter disabled", end="")
        print(f", >{a.max_nterm} residues before TM1" if a.max_nterm
              else ", N-terminal filter disabled")

    # the retinal-binding lysine column (rhodopsin K296) -- reported, not filtered
    # on: the input may deliberately carry non-opsin GPCRs as outgroups
    ref_pos, k296_col = 0, None
    for i, ch in enumerate(aln[refname]):
        if ch.isalpha():
            ref_pos += 1
            if ref_pos == 296:
                k296_col = i
                break

    rows, keep, review, nopred, dropped = [], [], [], [], []
    for n, s in clean:
        key = n.split()[0]
        if key not in aln:
            continue
        cov = [sum(1 for c in w if aln[key][c].isalpha()) / len(w)
               for w in windows]
        pass_hmm = min(cov) >= a.hmm_cov

        ntm = counts.get(key) if counts is not None else None
        hits = (helix_at_windows(aln[key], windows, intervals.get(key, []))
                if (counts is not None and key in counts) else None)
        if counts is None:
            pass_dtm = None
        elif ntm is None:
            pass_dtm = None          # no record: a data gap, never a failure
        elif a.dtm_mode == "count":
            pass_dtm = ntm >= a.min_tm
        else:
            pass_dtm = all(hits)

        ranges = window_residue_ranges(aln[key], windows)
        n_before = (ranges[0][0] - 1) if ranges[0] else None
        n_after = (len(s) - ranges[6][1]) if ranges[6] else None
        extras = (extra_helix_locations(ranges, intervals.get(key, []))
                  if counts is not None and key in counts else None)
        has_sp = (int(any(t == "signal" for t, _, _ in dtm_regions.get(key, [])))
                  if dtm_regions is not None else "")

        # Every reason a sequence is removed is recorded, and all removals share
        # the single verdict "drop". A sequence can trip more than one rule --
        # an incomplete core AND a large N-terminal region, say -- and then
        # every applicable reason appears, joined by '+'.
        #
        #   <7_TM              both criteria call the 7TM core incomplete
        #                      (profile coverage below --hmm-cov on at least one
        #                      window AND no predicted helix at that position)
        #   >7_TM              DeepTMHMM's raw helix count exceeds --max-tm
        #   >Naa_before_TM1    more than --max-nterm residues precede the TM1
        #                      window (the fusion / mis-assembly signature)
        reasons = []
        if counts is not None and pass_dtm is not None \
                and not pass_hmm and not pass_dtm:
            reasons.append("<7_TM")
        if a.max_tm and ntm is not None and ntm > a.max_tm:
            reasons.append(f">{a.max_tm}_TM")
        if a.max_nterm and n_before is not None and n_before > a.max_nterm:
            reasons.append(f">{a.max_nterm}aa_before_TM1")

        if counts is None:
            verdict = "keep" if pass_hmm else "drop"
            if not pass_hmm:
                reasons = ["<7_TM"]
        elif pass_dtm is None:
            verdict = "no_prediction"
        elif reasons:
            verdict = "drop"
        elif pass_hmm and pass_dtm:
            verdict = "keep"
        else:
            verdict = "review"

        row = dict(name=n, length=len(s))
        row.update({f"cov_TM{i+1}": round(cov[i], 3) for i in range(7)})
        row["min_cov"] = round(min(cov), 3)
        row["n_tm_helices"] = "" if ntm is None else ntm
        row["windows_with_helix"] = "" if hits is None else sum(hits)
        row["helix_pattern"] = ("" if hits is None
                                else "".join("1" if h else "0" for h in hits))
        row["n_before_TM1"] = "" if n_before is None else n_before
        row["n_after_TM7"] = "" if n_after is None else n_after
        row["extra_helix_location"] = ("" if not extras
                                       else "+".join(sorted(set(extras))))
        row["has_signal_peptide"] = has_sp
        row["res_at_K296"] = (aln[key][k296_col].upper()
                              if k296_col is not None else "")
        row["pass_hmm"] = int(pass_hmm)
        row["pass_deeptmhmm"] = "" if pass_dtm is None else int(pass_dtm)
        row["drop_reason"] = "+".join(reasons)
        row["verdict"] = verdict
        rows.append(row)

        if verdict == "keep":
            keep.append((n, s))
        elif verdict == "review":
            review.append((n, s))
        elif verdict == "no_prediction":
            nopred.append((n, s))
        else:
            # every removal, whatever the rule, lands in one file; the report's
            # drop_reason column says which rule or rules removed each sequence
            dropped.append((n, s))

    write_fasta(f"{a.outprefix}.faa", keep)
    if review:
        write_fasta(f"{a.outprefix}_review.faa", review)
    if nopred:
        write_fasta(f"{a.outprefix}_no_prediction.faa", nopred)
    if dropped:
        write_fasta(f"{a.outprefix}_dropped.faa", dropped)
    with open(f"{a.outprefix}_report.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print(f"\nkept          : {len(keep)}  -> {a.outprefix}.faa")
    print(f"review        : {len(review)}  -> {a.outprefix}_review.faa")
    if dropped:
        print(f"dropped       : {len(dropped)}  -> {a.outprefix}_dropped.faa")
        by_reason = {}
        for r in rows:
            if r["verdict"] == "drop":
                by_reason[r["drop_reason"]] = by_reason.get(r["drop_reason"], 0) + 1
        for reason, cnt in sorted(by_reason.items(), key=lambda z: -z[1]):
            print(f"                {cnt:4d}  {reason}")
    if nopred:
        print(f"no prediction : {len(nopred)}  -> {a.outprefix}_no_prediction.faa")
    print(f"report        : {a.outprefix}_report.csv")
    n_int = sum(1 for r in rows if "interior" in str(r["extra_helix_location"]))
    if n_int:
        print(f"\nNOTE: {n_int} sequences carry an extra TM helix INTERIOR to the "
              f"7TM core -- the fusion signature. See extra_helix_location.")


if __name__ == "__main__":
    main()
