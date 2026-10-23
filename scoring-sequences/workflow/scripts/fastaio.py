"""FASTA reading shared by the workflow scripts."""


def read_fasta(path):
    """Return [(header, sequence), ...] in file order. Header excludes '>'."""
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
