#!/usr/bin/env python3
"""Convert an LTRDeNovo (NGSEP) GFF3 to the pipeline's 16-column BED.

LTRDeNovo predicts LTR retrotransposons only, and publishes coordinates as a
GFF3 with `source=NGSEP`.

FEATURE SELECTION. The GFF3 describes each insertion at two nested levels plus
optional sub-parts:

    repeat_region           the full element footprint, LTR to LTR
    transposable_element    the INTERNAL domain only, between the two LTRs
    five_prime_LTR /        the individual LTRs (structural calls only)
    three_prime_LTR
    target_site_duplication the TSD (structural calls only)

This script keeps `repeat_region`, one record per insertion. That choice is
load-bearing: for the 228 structurally-detected elements, `transposable_element`
is strictly *inside* `repeat_region` -- e.g. TE_12 spans 2500121-2504726 as a
repeat_region but only 2500361-2504487 as a transposable_element, with the
5' LTR at 2500121-2500360 and the 3' LTR at 2504488-2504726. Taking the
transposable_element interval would leave both LTRs unannotated, which for an
LTR-specialist tool is precisely the wrong half of the element to drop. For the
2,876 homology calls the two levels are coordinate-identical, so the choice only
matters for the structural ones -- but it always matters in the same direction.

Keeping repeat_region also avoids double-counting: emitting both levels would
report the same insertion twice and inflate coverage and tool support.

DIVERGENCE. LTRDeNovo writes no divergence, identity, or score for any feature
(GFF3 columns 6 and 8 are `.` throughout, and no attribute carries one). Every
record is therefore emitted with perc_div = NA, and the tool is registered
`rm_fields=no` in config/tools.tsv so it contributes to existence and
classification but not to the divergence signal.

CLASSIFICATION. The `classification` attribute uses Wicker et al. 2007 codes
under an `LTR/` prefix -- LTR/RLG (Gypsy), LTR/RLC (Copia), LTR/RLR
(Retrovirus) -- plus LTR/Unknown. These pass through unchanged; config/class_map.tsv
resolves them.

METHOD. `method=homology|structural` is preserved in the name field as a
suffix, so the two detection modes stay distinguishable in the browser without
adding a column outside the BED16 contract.

Usage:
    python scripts/ltrdenovogff2bed.py GCA_951799975.1_LTRDeNovo.gff.gz -o ltrdenovo.bed
"""

import argparse
import gzip
import sys

KEEP_TYPE = "repeat_region"


def attrs(field):
    out = {}
    for kv in field.strip().strip(";").split(";"):
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
        if f[2] != KEEP_TYPE:
            skipped += 1
            continue

        start, end = int(f[3]) - 1, int(f[4])   # GFF3 is 1-based inclusive
        if end <= start:
            empty += 1
            continue

        a = attrs(f[8])
        hit_id = a.get("ID", "")
        method = a.get("method", "")
        name = f"{hit_id}_{method}" if method else hit_id
        strand = f[6] if f[6] in ("+", "-") else "."

        outfile.write("\t".join([
            f[0], str(start), str(end),
            name, "0", strand,
            "NA",                                  # SW_score: not reported
            "NA", "NA", "NA", "NA",                # perc_div/del/ins, query_left
            a.get("classification", "Unknown"),
            "NA", "NA", "NA",                      # repeat_start/end/left
            hit_id,
        ]) + "\n")
        kept += 1

    print(f"[ltrdenovogff2bed] kept {kept} {KEEP_TYPE}(s); "
          f"skipped {skipped} other feature line(s)", file=sys.stderr)
    if empty:
        print(f"[ltrdenovogff2bed] skipped {empty} zero-length interval(s)",
              file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("gff", help="LTRDeNovo GFF3 (optionally .gz)")
    ap.add_argument("-o", "--out", required=True, help="output BED16")
    args = ap.parse_args()
    with opener(args.gff) as fh, open(args.out, "w") as out:
        convert(fh, out)


if __name__ == "__main__":
    main()