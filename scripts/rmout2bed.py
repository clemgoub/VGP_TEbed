#!/usr/bin/env python3
"""
Convert a RepeatMasker .out file to BED format.

Mandatory BED6 fields:
    chrom, chromStart (0-based), chromEnd, name, score, strand

Extra columns appended (RepeatMasker fields not captured by BED6):
    SW_score, perc_div, perc_del, perc_ins, query_left,
    repeat_class_family, repeat_start, repeat_end, repeat_left, RM_ID

Usage:
    python rmout2bed.py input.out > output.bed
    python rmout2bed.py input.out -o output.bed
"""

import argparse
import sys


def parse_int_or_paren(s):
    """RepeatMasker writes '(left)' values in parentheses. Strip them."""
    return int(s.strip().lstrip("(").rstrip(")"))


def convert(infile, outfile):
    header = (
        "#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\t"
        "SW_score\tperc_div\tperc_del\tperc_ins\tquery_left\t"
        "repeat_class_family\trepeat_start\trepeat_end\trepeat_left\tRM_ID\n"
    )
    outfile.write(header)

    n_malformed = 0      # lines with too few fields
    n_empty = 0          # records that would yield a zero-length/negative interval

    for raw in infile:
        line = raw.strip()
        if not line:
            continue
        # Skip RepeatMasker header lines
        if line.startswith(("SW", "score")):
            continue

        f = line.split()
        if len(f) < 15:
            # malformed / unexpected line
            n_malformed += 1
            continue

        sw_score    = f[0]
        perc_div    = f[1]
        perc_del    = f[2]
        perc_ins    = f[3]
        chrom       = f[4]
        q_begin     = int(f[5])           # 1-based, inclusive
        q_end       = int(f[6])           # inclusive
        q_left      = parse_int_or_paren(f[7])
        rm_strand   = f[8]                # '+' or 'C'
        repeat_name = f[9]
        repeat_cf   = f[10]

        # Repeat coordinates: order depends on strand.
        # '+' strand: begin end (left)
        # 'C' strand: (left) end begin
        if rm_strand == "+":
            r_begin = parse_int_or_paren(f[11])
            r_end   = parse_int_or_paren(f[12])
            r_left  = parse_int_or_paren(f[13])
            strand  = "+"
        else:  # 'C'
            r_left  = parse_int_or_paren(f[11])
            r_end   = parse_int_or_paren(f[12])
            r_begin = parse_int_or_paren(f[13])
            strand  = "-"

        rm_id = f[14]

        # BED is 0-based, half-open
        chrom_start = q_begin - 1
        chrom_end   = q_end

        # Guard against zero-length / negative intervals (chromEnd <= chromStart).
        # These are invalid BED features: they survive `bedtools merge` untouched
        # and make `bedtools intersect` counts orientation-dependent. Drop them.
        if chrom_end <= chrom_start:
            n_empty += 1
            continue

        # BED score must be 0..1000; cap SW score for that column,
        # but keep the raw value in the extras.
        bed_score = min(int(sw_score), 1000)

        out_fields = [
            chrom,
            str(chrom_start),
            str(chrom_end),
            repeat_name,
            str(bed_score),
            strand,
            sw_score,
            perc_div,
            perc_del,
            perc_ins,
            str(q_left),
            repeat_cf,
            str(r_begin),
            str(r_end),
            str(r_left),
            rm_id,
        ]
        outfile.write("\t".join(out_fields) + "\n")

    if n_malformed:
        print(f"[rmout2bed] skipped {n_malformed} malformed line(s) (<15 fields)",
              file=sys.stderr)
    if n_empty:
        print(f"[rmout2bed] skipped {n_empty} zero-length/negative interval(s) "
              "(chromEnd <= chromStart)", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description="Convert RepeatMasker .out to BED.")
    ap.add_argument("infile", help="RepeatMasker .out file")
    ap.add_argument("-o", "--output", help="Output BED file (default: stdout)")
    args = ap.parse_args()

    with open(args.infile) as fh:
        if args.output:
            with open(args.output, "w") as out:
                convert(fh, out)
        else:
            convert(fh, sys.stdout)


if __name__ == "__main__":
    main()
