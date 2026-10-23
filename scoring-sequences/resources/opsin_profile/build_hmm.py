#!/usr/bin/env python3
"""Build a profile HMM from a bait alignment (pyhmmer, default Builder settings).

Parameterised version of the original build_hmm.py, which hardcoded its paths.
Defaults reproduce that script exactly: pyhmmer's `Builder` with no arguments,
which mirrors hmmbuild's defaults (--fast construction at symfrac 0.5, --wpb
position-based sequence weighting, --eent relative-entropy weighting).

Note on entropy weighting: --eent reduces the effective sequence count to hit a
target relative entropy, so a 23-sequence seed reports EFFN ~2.6. That flattens
every column, which is what makes the model tolerant of divergent residues. If
you ever need to sharpen or flatten TM7 specifically, that is not reachable from
here -- it is a property of the seed's composition at those columns.

Usage
-----
    python build_hmm.py baits_090326.aln -o animalOpsins_090326.hmm \\
        --name opsins_type_II

The input must be an ALIGNMENT (e.g. mafft --ep 0 --genafpair --maxiterate 1000,
i.e. E-INS-i), not an unaligned FASTA.
"""
import argparse
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("alignment", help="aligned FASTA (or any format easel reads)")
    ap.add_argument("-o", "--out", required=True, help="output .hmm path")
    ap.add_argument("--name", default="opsins_type_II", help="model NAME field")
    a = ap.parse_args()

    try:
        import pyhmmer
    except ImportError:
        sys.exit("pyhmmer not available -- activate the env that has it, or use:\n"
                 f"  hmmbuild -n {a.name} --amino {a.out} {a.alignment}")

    alphabet = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.easel.MSAFile(a.alignment, digital=True,
                               alphabet=alphabet) as handle:
        msa = handle.read()
    if msa is None:
        sys.exit(f"could not read an alignment from {a.alignment}")

    msa.name = a.name.encode() if isinstance(a.name, str) else a.name
    builder = pyhmmer.plan7.Builder(alphabet)
    background = pyhmmer.plan7.Background(alphabet)
    hmm, _, _ = builder.build_msa(msa, background)

    with open(a.out, "wb") as fh:
        hmm.write(fh)

    print(f"built {a.out}")
    print(f"  name  {hmm.name.decode() if hmm.name else '-'}")
    print(f"  M     {hmm.M} match states")
    print(f"  nseq  {msa.sequences.__len__()}")


if __name__ == "__main__":
    main()
