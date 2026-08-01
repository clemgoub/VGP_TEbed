#!/usr/bin/env python3
"""Convert an EDTA TEanno GFF3 to the pipeline's 16-column BED.

EDTA publishes coordinates only as GFF3 (`*.fa.mod.EDTA.TEanno.gff3`), so this
is the entry point for rebuilding the hub from the GenomeArk copy.

FEATURE SELECTION. The GFF3 mixes whole elements with their sub-parts. A
structurally-detected LTR element is emitted as a `repeat_region` container
holding `lTSD`/`rTSD` target-site duplications, an `LTRRT` element, and its
`lLTR`/`rLTR` long terminal repeats -- six records describing one insertion.
Keeping all of them would count the same locus up to three times and inflate
both coverage and tool support. This script keeps one record per element:

    TE_homo_*     homology-based hits (the bulk)
    TE_struc_*    structurally-detected non-LTR TEs
    LTRRT_*       structurally-detected LTR retrotransposons (the element,
                  not its container, TSDs, or individual LTRs)

DIVERGENCE. EDTA reports `identity` (0-1 to the library consensus) rather than
RepeatMasker's percent divergence. Homology calls are converted as
perc_div = 100 * (1 - identity), which is the same quantity RepeatMasker's
column reports. Structural calls have no library alignment -- EDTA still writes
an `identity` for them, but it is an LTR-to-LTR identity (a within-element
measure), NOT divergence from a consensus. Those are emitted as NA so they
cannot contaminate the mean-divergence signal; this reproduces the NA pattern
of the RepeatMasker-derived EDTA input.

Usage:
    python scripts/edtagff2bed.py TEanno.gff3 -o edta.bed
"""

import argparse
import gzip
import sys

KEEP_PREFIXES = ("TE_homo", "TE_struc", "LTRRT")


def attrs(field):
    out = {}
    for kv in field.rstrip(";").split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def opener(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def convert(infile, outfile):
    outfile.write(
        "#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\t"
        "SW_score\tperc_div\tperc_del\tperc_ins\tquery_left\t"
        "repeat_class_family\trepeat_start\trepeat_end\trepeat_left\thit_id\n"
    )
    kept = skipped = empty = 0
    for line in infile:
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        a = attrs(f[8])
        hit_id = a.get("ID", "")
        if not hit_id.startswith(KEEP_PREFIXES):
            skipped += 1
            continue

        start, end = int(f[3]) - 1, int(f[4])   # GFF3 is 1-based inclusive
        if end <= start:
            empty += 1
            continue

        structural = a.get("method", "") == "structural"
        if structural:
            perc_div = "NA"
        else:
            try:
                perc_div = f"{100.0 * (1.0 - float(a['identity'])):.1f}"
            except (KeyError, ValueError):
                perc_div = "NA"

        try:
            score = min(int(float(f[5])), 1000)
        except ValueError:
            score = 0
        strand = f[6] if f[6] in ("+", "-") else "."

        outfile.write("\t".join([
            f[0], str(start), str(end),
            a.get("Name", hit_id), str(score), strand,
            f[5] if f[5] != "." else "NA",   # SW_score
            perc_div, "NA", "NA", "NA",      # perc_del, perc_ins, query_left
            a.get("classification", "Unknown"),
            "NA", "NA", "NA",                # repeat_start/end/left
            hit_id,
        ]) + "\n")
        kept += 1

    print(f"[edtagff2bed] kept {kept} element(s); skipped {skipped} sub-feature "
          f"record(s) (repeat_region/lLTR/rLTR/lTSD/rTSD)", file=sys.stderr)
    if empty:
        print(f"[edtagff2bed] skipped {empty} zero-length interval(s)", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gff3", help="EDTA TEanno GFF3 (plain or .gz)")
    ap.add_argument("-o", "--output", default="-", help="output BED (default stdout)")
    args = ap.parse_args()

    with opener(args.gff3) as fh:
        if args.output == "-":
            convert(fh, sys.stdout)
        else:
            with open(args.output, "w") as out:
                convert(fh, out)


if __name__ == "__main__":
    main()
