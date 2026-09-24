#!/usr/bin/env python3
"""Cache bovine rhodopsin (UniProt P02699) with its annotated TM spans.

filter_7tm.py will fetch this itself if it is missing, but then every filter
run is a network call and an offline rerun fails inside the filter rather than
in a rule of its own. Making it an explicit one-line rule means the reference
is a tracked input: fetched once, reused by every sample, and visible in the
DAG as the thing that defines the seven TM windows.

The payload is written byte-for-byte as filter_7tm.fetch_reference() expects:
the `ft_transmem,sequence` projection of the UniProtKB JSON record.
"""
import argparse
import json
import os
import sys
import urllib.request

URL = ("https://rest.uniprot.org/uniprotkb/P02699.json"
       "?fields=ft_transmem,sequence")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--out", required=True, help="output P02699.json path")
    a = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    tmp = a.out + ".part"
    urllib.request.urlretrieve(URL, tmp)

    with open(tmp) as fh:
        d = json.load(fh)
    spans = [f for f in d.get("features", []) if f["type"] == "Transmembrane"]
    if len(spans) != 7:
        os.remove(tmp)
        sys.exit(f"expected 7 annotated TM spans in P02699, got {len(spans)}. "
                 f"The UniProt annotation changed; the seven windows the filter "
                 f"is built on are no longer defined by this record.")
    if not d.get("sequence", {}).get("value"):
        os.remove(tmp)
        sys.exit("P02699 record carries no sequence")

    os.replace(tmp, a.out)
    print(f"{a.out}: {len(d['sequence']['value'])} aa, {len(spans)} TM spans")


if __name__ == "__main__":
    main()
