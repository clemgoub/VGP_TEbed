#!/usr/bin/env python3
"""Convert Satellome's 5-column BED to the pipeline's 16-column BED.

Satellome reports satellite arrays as BED5:

    chrom  chromStart  chromEnd  family  length

Column 5 is redundant (verified equal to chromEnd - chromStart on every row of
the pilot file) and is dropped. Coordinates are already 0-based half-open BED.

FAMILY AND CLASS. Column 4 is the satellite family (`fGobNig19A`, assembly-
prefixed with a serial and an `A` suffix). It goes to `name` unchanged.
`repeat_class_family` is `Satellite` for every row: detecting satellite arrays
is the tool's entire method, so unlike a period-based inference (see
fastan2bed.py --classify-period) this is the tool's own assertion, not a
converter guess. `Satellite` maps to repeat:tandem:satellite in
config/class_map.tsv for all tools.

NO SCORES OR DIVERGENCE. The native output carries no score, identity, or
divergence of any kind -- every RepeatMasker-style column is emitted NA and the
tool is registered rm_fields=no.

STRAND. Satellite arrays have no meaningful orientation in this output;
emitted as ".".

MINIMUM ARRAY LENGTH. The pilot file (`*.10kb.bed`) contains only arrays
>= 10 kb (verified: minimum observed 10,005 bp). Absence of a Satellome call
therefore does not mean absence of satellite -- shorter arrays were filtered
upstream. Recorded in the manifest notes; the `tandem` scope already keeps the
tool out of the eligibility denominator for non-tandem loci.

HIT_ID. Arrays are independent insertions, not fragments of one, so hit_id is
`<family>_<n>` with a per-family serial -- unique per tool per assembly, and
groupable by family prefix.

Usage:
    python scripts/satellome2bed.py GCA_..._genomic.10kb.bed -o satellome.bed
"""

import argparse
import collections
import gzip
import sys


def opener(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def convert(infile, outfile):
    outfile.write(
        "#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\t"
        "SW_score\tperc_div\tperc_del\tperc_ins\tquery_left\t"
        "repeat_class_family\trepeat_start\trepeat_end\trepeat_left\thit_id\n"
    )
    kept = bad = 0
    serial = collections.Counter()
    for line in infile:
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 4:
            bad += 1
            continue
        chrom, start, end, fam = f[0], int(f[1]), int(f[2]), f[3]
        if end <= start:
            bad += 1
            continue
        serial[fam] += 1
        outfile.write("\t".join([
            chrom, str(start), str(end),
            fam, "0", ".",
            "NA", "NA", "NA", "NA", "NA",   # SW, div, del, ins, query_left
            "Satellite",
            "NA", "NA", "NA",               # repeat_start/end/left
            f"{fam}_{serial[fam]}",
        ]) + "\n")
        kept += 1

    print(f"[satellome2bed] kept {kept} array(s), "
          f"{len(serial)} families", file=sys.stderr)
    if bad:
        print(f"[satellome2bed] skipped {bad} malformed line(s)", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("bed", help="Satellome BED5 (optionally .gz)")
    ap.add_argument("-o", "--out", required=True, help="output BED16")
    args = ap.parse_args()
    with opener(args.bed) as fh, open(args.out, "w") as out:
        convert(fh, out)


if __name__ == "__main__":
    main()