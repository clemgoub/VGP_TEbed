#!/usr/bin/env python3
"""Recover the sequence-name mapping for a RepeatMasker .out run on a renamed FASTA.

Some upstream runs are performed on EDTA's `.fa.mod`, which rewrites sequence
names to short opaque ids (`_J0000000`, `_J000001q`, ...). The resulting .out
cannot be joined to the assembly, and no mapping file is published alongside it.

RepeatMasker records `end` and `(left)` for every hit, and their sum is the full
length of the query sequence. That gives an exact length for each renamed
sequence without needing the modified FASTA, which is enough to match against
the assembly's chrom.sizes -- uniquely, wherever the length is unique.

Where several assembly sequences share a length (small unplaced scaffolds with
round sizes are the usual case) the mapping is genuinely undetermined by this
evidence, and this script does NOT guess. Those rows are written with an
`AMBIGUOUS:` marker so the mapping fails loudly if used unedited, and the
candidates are listed so a human can resolve or discard them.

Sequence ids look base-62 sequential in FASTA order, but assembly FASTA order
is not chrom.sizes order, nor NCBI assembly-report order (both were checked and
neither reproduces the confidently-mapped anchors), so interpolating the
ambiguous ids from their neighbours is not sound either.

Usage:
    python scripts/infer_rm_rename.py in.out assembly.chrom.sizes -o rename.tsv
"""

import argparse
import sys
from collections import defaultdict


def seq_lengths(path):
    """{renamed sequence: length} from end + (left), verifying consistency."""
    seen = defaultdict(set)
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 15 or not f[0].isdigit():
                continue
            seen[f[4]].add(int(f[6]) + int(f[7].strip("()")))
    lengths, bad = {}, {}
    for s, vals in seen.items():
        if len(vals) == 1:
            lengths[s] = vals.pop()
        else:
            bad[s] = sorted(vals)
    return lengths, bad


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rmout", help="RepeatMasker .out run on the renamed FASTA")
    ap.add_argument("sizes", help="assembly chrom.sizes (name<TAB>length)")
    ap.add_argument("-o", "--output", default="-", help="output TSV (default stdout)")
    args = ap.parse_args()

    lengths, inconsistent = seq_lengths(args.rmout)
    if inconsistent:
        print(f"[infer_rm_rename] {len(inconsistent)} sequence(s) gave inconsistent "
              "lengths across their hits; refusing to map them", file=sys.stderr)

    by_len = defaultdict(list)
    with open(args.sizes) as fh:
        for line in fh:
            if not line.strip():
                continue
            name, size = line.split()[:2]
            by_len[int(size)].append(name)

    out = sys.stdout if args.output == "-" else open(args.output, "w")
    n_ok = n_amb = n_missing = 0
    try:
        out.write("# renamed\tassembly\t# inferred by sequence length "
                  "(RepeatMasker end + (left))\n")
        for s in sorted(lengths, key=lambda k: (len(k), k)):
            cands = by_len.get(lengths[s], [])
            if len(cands) == 1:
                out.write(f"{s}\t{cands[0]}\n")
                n_ok += 1
            elif not cands:
                out.write(f"# NO_MATCH\t{s}\tlength={lengths[s]} absent from sizes\n")
                n_missing += 1
            else:
                out.write(f"{s}\tAMBIGUOUS:{','.join(cands)}\t"
                          f"# length={lengths[s]} shared by {len(cands)} sequences\n")
                n_amb += 1
    finally:
        if out is not sys.stdout:
            out.close()

    print(f"[infer_rm_rename] {n_ok} unique, {n_amb} ambiguous, {n_missing} unmatched "
          f"of {len(lengths)} renamed sequences", file=sys.stderr)
    if n_amb:
        print("[infer_rm_rename] resolve or delete the AMBIGUOUS rows before use; "
              "rmout2bed.py --drop-unmapped discards records on unmapped sequences",
              file=sys.stderr)


if __name__ == "__main__":
    main()
