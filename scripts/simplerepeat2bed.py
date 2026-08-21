#!/usr/bin/env python3
"""Convert UCSC GenArk's simpleRepeat (TRF) track to the pipeline's BED16.

Input is the text dump of the bigBed:

    bigBedToBed GCA_..._fGobNig1.1.simpleRepeat.bb stdout | simplerepeat2bed.py -

16 columns per the simpleRepeat autoSql: chrom, chromStart, chromEnd, name
(repeat-unit tag), period, copyNum, consensusSize, perMatch, perIndel, score,
A, C, G, T, entropy, sequence.

CLASS. Emits `tandem` -> repeat:tandem for every row (existing `tandem*` rule).
TRF finds tandem arrays of ANY period, from dinucleotide microsatellites to
multi-kb satellite units, so asserting Simple_repeat (or Satellite) from period
alone would be a converter guess, not a tool assertion -- the same reasoning as
fastan2bed.py, whose --classify-period stays opt-in. This file's periods run
1 bp to multi-kb (period > 6 on 88% of rows), so the guess would be wrong often.

FIELDS KEPT. name = `p<period>_x<copyNum>_<unit tag>` so period and copy number
reach the per-tool mouseover (unit tag truncated to 24 bp to keep bigBed names
sane; the full unit sequence is not carried). SW_score = TRF alignment score,
uncapped. perc_div = 100 - perMatch: like FasTAN this is unit-to-unit identity
within the array, NOT divergence from a library consensus -- register the tool
rm_fields=divergence_only. perIndel is a combined indel percentage that maps to
neither perc_del nor perc_ins alone; carried in neither (NA both), noted here.

NOT KEPT. Base composition (A/C/G/T percents), entropy, consensusSize, and the
full unit sequence. All are derivable from the public GenArk track this file
came from; the residual-losses list in INPUT_FORMAT.md section 6 points there.

STRAND. TRF arrays have no orientation; ".".

Usage:
    ./bin/bigBedToBed simpleRepeat.bb stdout | \\
        python scripts/simplerepeat2bed.py - -o inputs/trf.bed
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
        if len(f) < 16:
            bad += 1
            continue
        chrom, start, end = f[0], int(f[1]), int(f[2])
        if end <= start:
            bad += 1
            continue
        tag, period, copynum = f[3], f[4], f[5]
        permatch, score = f[7], f[9]
        div = f"{100 - int(permatch):.1f}" if permatch.isdigit() else "NA"
        name = f"p{period}_x{copynum}_{tag[:24]}"
        outfile.write("\t".join([
            chrom, str(start), str(end),
            name,
            str(min(int(score), 1000)) if score.isdigit() else "0",
            ".",
            score if score.isdigit() else "NA",   # SW_score: raw TRF score
            div, "NA", "NA", "NA",                # div, del, ins, query_left
            "tandem",
            "NA", "NA", "NA",
            f"trf_{i}",
        ]) + "\n")
        kept += 1
    print(f"[simplerepeat2bed] kept {kept} array(s)", file=sys.stderr)
    if bad:
        print(f"[simplerepeat2bed] skipped {bad} malformed line(s)", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("bed", help="bigBedToBed dump of simpleRepeat.bb, or - for stdin")
    ap.add_argument("-o", "--out", required=True)
    args = ap.parse_args()
    fh = sys.stdin if args.bed == "-" else open(args.bed)
    with open(args.out, "w") as out:
        convert(fh, out)
    if fh is not sys.stdin:
        fh.close()


if __name__ == "__main__":
    main()