"""Content-addressed store of DeepTMHMM predictions, one record per sequence.

WHY. DeepTMHMM predicts each sequence independently -- a batch exists only to
amortise loading the ESM-1b model. Caching at batch granularity would tie the
unit of recomputation to that: adding five sequences to a dataset would re-run
the whole batch they landed in. This store makes the sequence the unit instead.

KEY. sha1 of the sequence itself (uppercased, gaps stripped), not the FASTA
header. That is the input DeepTMHMM actually saw, so it is the only key under
which a cached prediction is certainly still valid:

  * a header renamed between assemblies reuses the prediction, correctly;
  * a header REUSED for a changed sequence misses the cache, correctly --
    a name-keyed store would silently serve a stale topology here;
  * the same sequence entered twice under different headers (the AFK83790
    double-entry README.md flags) is predicted once and emitted twice.

Headers are therefore not stored in a record. They are reattached at collect
time from whichever FASTA is being processed, which is also what lets one store
serve several datasets.

RECORD. JSON, one file per unique sequence:

    {"sha1": ..., "length": 348, "n_tm_helices": 7,
     "comments":  ["Length: 348", "Number of predicted TMRs: 7"],
     "features":  [["signal", 1, 28], ["outside", 29, 36], ...],
     "topology":  "SSSS...oooMMMM...",
     "source":    "DeepTMHMM 1.0 local", "batch": "batch_003",
     "ingested":  "2026-09-23T17:40:00Z"}

`comments` holds the text of each `# <name> ...` line with the name stripped,
so whatever DeepTMHMM emits round-trips even if it adds fields later; `features`
and `topology` are the two representations the batch output carries.

CONCURRENCY. Distinct sequences hash to distinct filenames, so parallel batch
jobs never write the same path and no locking is needed. Writes go through a
temp file and os.replace, so a killed job leaves either nothing or a complete
record -- never a half-written one a later run would trust.
"""
import hashlib
import json
import os
import time

# Default provenance tag. Records imported from the BioLib web service carry
# "biolib:DTU/DeepTMHMM" instead (dtm_ingest_legacy.py --source).
SOURCE = "DeepTMHMM 1.0 local"


# ------------------------------------------------------------------ keying

def seq_key(seq):
    """sha1 of the sequence as DeepTMHMM saw it: uppercase, no gaps."""
    norm = seq.upper().replace("-", "").replace(".", "").replace("*", "")
    return hashlib.sha1(norm.encode()).hexdigest()


def record_path(cache_dir, key):
    return os.path.join(cache_dir, f"{key}.json")


def has(cache_dir, key):
    return os.path.exists(record_path(cache_dir, key))


def load(cache_dir, key):
    with open(record_path(cache_dir, key)) as fh:
        return json.load(fh)


def save(cache_dir, rec):
    os.makedirs(cache_dir, exist_ok=True)
    path = record_path(cache_dir, rec["sha1"])
    tmp = path + f".part{os.getpid()}"
    with open(tmp, "w") as fh:
        json.dump(rec, fh, sort_keys=True)
        fh.write("\n")
    os.replace(tmp, path)
    return path


# ------------------------------------------- reading DeepTMHMM batch output

class FormatError(RuntimeError):
    """Raised when batch output does not look like DeepTMHMM wrote it.

    Loud on purpose. A silently mis-parsed topology file would populate the
    store with wrong records that every later run would trust without
    re-running anything.
    """


def parse_3line(path):
    """{name: (sequence, topology)} from predicted_topologies.3line.

    The 3line file carries the SEQUENCE, which is what makes the store
    self-sufficient: a batch directory can be ingested without the FASTA it
    came from.
    """
    out, block = {}, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if line.startswith(">"):
                if block:
                    _commit_3line(out, block, path)
                block = [line]
            else:
                block.append(line)
    if block:
        _commit_3line(out, block, path)
    if not out:
        raise FormatError(f"{path}: no records found")
    return out


def _commit_3line(out, block, path):
    if len(block) != 3:
        raise FormatError(
            f"{path}: expected a 3-line record (header, sequence, topology), "
            f"got {len(block)} lines starting {block[0][:60]!r}")
    # ">name | TM" -- the annotation after '|' is DeepTMHMM's class call
    name = block[0][1:].strip().split("|")[0].strip().split()[0]
    seq, topo = block[1].strip(), block[2].strip()
    if len(seq) != len(topo):
        raise FormatError(f"{path}: {name}: sequence is {len(seq)} residues but "
                          f"topology string is {len(topo)}")
    out[name] = (seq, topo)


def parse_gff3(path):
    """{name: {'comments': [...], 'features': [(type, start, end), ...]}}.

    Mirrors filter_7tm.parse_dtm's reading of the same file, with the name
    stripped out of both line kinds so records are header-independent.
    """
    out = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            # "##" lines are GFF pragmas (DeepTMHMM writes "##gff-version 3"
            # first), not per-sequence comments -- read as one, the pragma
            # would become a sequence named "#gff-version".
            if not line.strip() or line.startswith("//") \
                    or line.startswith("##"):
                continue
            if line.startswith("#"):
                body = line[1:].strip()
                name = body.split()[0] if body.split() else None
                if name is None:
                    continue
                rest = body[len(name):].strip()
                out.setdefault(name, {"comments": [], "features": []})
                out[name]["comments"].append(rest)
            else:
                p = line.split("\t")
                if len(p) < 4:
                    raise FormatError(f"{path}: expected 4+ tab-separated "
                                      f"fields, got {len(p)}: {line[:80]!r}")
                name = p[0].split()[0]
                out.setdefault(name, {"comments": [], "features": []})
                out[name]["features"].append((p[1], int(p[2]), int(p[3])))
    if not out:
        raise FormatError(f"{path}: no records found")
    return out


def ingest_batch(gff3, line3, cache_dir, batch_label, source=SOURCE):
    """Split one batch's output into per-sequence records. Returns [keys].

    Every sequence in the 3line file must also appear in the gff3, and vice
    versa -- a mismatch means the two files describe different runs.
    """
    seqs = parse_3line(line3)
    feats = parse_gff3(gff3)
    only_3line = set(seqs) - set(feats)
    only_gff3 = set(feats) - set(seqs)
    if only_3line or only_gff3:
        raise FormatError(
            f"{batch_label}: {os.path.basename(line3)} and "
            f"{os.path.basename(gff3)} disagree on which sequences were run "
            f"(3line only: {sorted(only_3line)[:3]}, "
            f"gff3 only: {sorted(only_gff3)[:3]})")

    now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    keys = []
    for name, (seq, topo) in seqs.items():
        key = seq_key(seq)
        save(cache_dir, {
            "sha1": key,
            "length": len(seq),
            "n_tm_helices": sum(1 for t, _, _ in feats[name]["features"]
                                if t == "TMhelix"),
            "comments": feats[name]["comments"],
            "features": [list(f) for f in feats[name]["features"]],
            "topology": topo,
            "source": source,
            "batch": batch_label,
            "ingested": now,
        })
        keys.append(key)
    return keys


# ------------------------------------------------------- writing back out

def emit_gff3(records, out_path):
    """Reattach headers and write a single TMRs.gff3 filter_7tm.py can read.

    `records` is [(header_name, record), ...]. The same record may appear under
    several names -- that is the duplicate-sequence case, and each header gets
    its own block so the filter scores every input record.
    """
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write("##gff-version 3\n")
        for name, rec in records:
            for c in rec["comments"]:
                fh.write(f"# {name} {c}\n")
            for typ, start, end in rec["features"]:
                fh.write(f"{name}\t{typ}\t{start}\t{end}\n")
            fh.write("//\n")


def emit_3line(records, out_path, seqs_by_name):
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with open(out_path, "w") as fh:
        for name, rec in records:
            fh.write(f">{name}\n{seqs_by_name[name]}\n{rec['topology']}\n")
