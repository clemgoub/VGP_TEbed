#!/usr/bin/env python3
"""Convert UCSC GenArk's windowMasker (WM + SDust) track to the pipeline's BED16.

Input is the text dump of the bigBed 3 -- bare intervals, nothing else:

    bigBedToBed GCA_..._fGobNig1.1.windowMasker.bb stdout | windowmasker2bed.py -

WHAT THIS TOOL ASSERTS. WindowMasker flags over-represented k-mer windows
(plus low-complexity via the SDust pass merged into this track). That is a
statement that the sequence is REPETITIVE -- not that it is a TE, a tandem
array, or any class at all. The track cannot distinguish its WM intervals from
its SDust intervals, so even `low_complexity` cannot be asserted per row.

CLASS. Every row emits `Unknown` -> bare `repeat`: full existence support,
abstention from every classification vote (the map's UNINFORMATIVE rule).
This makes WindowMasker a pure detector -- exactly what a k-mer masker is.
Register scope=general_homology (it can flag any repetitive sequence) and
rm_fields=no (bigBed 3 carries no score, identity or divergence).

EVIDENCE. hit_id prefix `wm_` -- k-mer statistics, neither homology nor
structure; the manifest's evidence note should say so.

Usage:
    ./bin/bigBedToBed windowMasker.bb stdout | \\
        python scripts/windowmasker2bed.py - -o inputs/windowmasker.bed
"""

import argparse
import sys


def convert(infile, outfile):
    outfile.write(
        "#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\t"
        "SW_score\tperc_div\tperc_del\tperc_ins\tquery_left\t"
        "repeat_class_family\trepeat_start\trepeat_end\trepeat_left\thit_id\n"
    )
    kept = bad = 0
    for i, line in enumerate(infile, 1):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            bad += 1
            continue
        chrom, start, end = f[0], int(f[1]), int(f[2])
        if end <= start:
            bad += 1
            continue
        outfile.write("\t".join([
            chrom, str(start), str(end),
            "wm_masked", "0", ".",
            "NA", "NA", "NA", "NA", "NA",
            "Unknown",
            "NA", "NA", "NA",
            f"wm_{i}",
        ]) + "\n")
        kept += 1
    print(f"[windowmasker2bed] kept {kept} interval(s)", file=sys.stderr)
    if bad:
        print(f"[windowmasker2bed] skipped {bad} malformed line(s)", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("bed", help="bigBedToBed dump of windowMasker.bb, or - for stdin")
    ap.add_argument("-o", "--out", required=True)
    args = ap.parse_args()
    fh = sys.stdin if args.bed == "-" else open(args.bed)
    with open(args.out, "w") as out:
        convert(fh, out)
    if fh is not sys.stdin:
        fh.close()


if __name__ == "__main__":
    main()