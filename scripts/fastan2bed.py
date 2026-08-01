#!/usr/bin/env python3
"""
Convert native FasTAN BED output to the pipeline's BED16 input format.

FasTAN (https://github.com/thegenemyers/FASTAN) is a tandem array finder. Its
native BED has five columns and no header:

    1 chrom       sequence name
    2 chromStart  0-based, half-open
    3 chromEnd
    4 period      estimated average unit size, in bp
    5 identity    average identity of the first-wave alignment, per mille
                  (1000 = units identical to each other)

Both column-4 and column-5 semantics are from the FasTAN README: with -m it
"produces a .1ano file of the intervals containining the detected tandem arrays
along with an estimate of the average unit size and the average identity of the
alignment".

WHAT MAPS WHERE, AND WHY
------------------------
score            = identity, verbatim. Already on BED's 0-1000 scale and
                   monotone in confidence, so no rescaling is invented.
strand           = '.'  A tandem array has no meaningful orientation.
perc_div         = (1000 - identity) / 10

    *** This is NOT the same quantity as RepeatMasker's perc_div. ***
    RepeatMasker reports divergence of a copy from its library consensus, which
    is an age proxy. FasTAN reports divergence of array units FROM EACH OTHER,
    which is an array-homogeneity proxy. They correlate loosely and mean
    different things. It is written to perc_div anyway because it is a real
    sequence-divergence percentage and the summary track's mean-divergence is
    explicitly a mean over tools that report one -- but tools.tsv records
    rm_fields=divergence_only for fastan so the distinction survives in the
    hub documentation. Use --no-divergence to emit NA instead.

SW_score, perc_del, perc_ins, query_left, repeat_start, repeat_end,
repeat_left  = NA. FasTAN performs no library alignment, so there is no
               consensus coordinate system and no score in those units.
               NA is correct here; 0 would be a lie the pipeline cannot detect.

repeat_class_family
    Default: the literal string `tandem`, which config/class_map.tsv maps to
    `repeat:tandem` at medium confidence.

    This is deliberate abstention. FasTAN detects tandem arrays; it does not
    classify them as satellite vs simple vs low-complexity. Emitting a
    subclass by default would manufacture a classification the tool never
    made, and the summary track would then report class agreement or conflict
    on evidence that does not exist.

    `--classify-period` opts in to the conventional size cut (period <= 6 ->
    Simple_repeat, period > 6 -> Satellite). That convention is real and
    widely used, but it is OUR inference, not FasTAN's call -- hence opt-in.

    *** Consequence of enabling it, measured. *** Segmentation cannot tell an
    inferred subclass from one a tool actually reported, so it scores ours as a
    real vote -- in BOTH directions. On the three-chromosome goby slice, over
    segments whose boundaries are identical in both builds, among bases where
    FasTAN votes:

        agreement deepened (synthetic consensus)   7.84 Mb
        agreement shallowed                        1.25 Mb
        new conflict, asserted by no tool          1.25 Mb
        unchanged                                  9.22 Mb

    Mean agreeDepth moves 3.135 -> 3.207 genome-wide. Turn the flag on only if
    you want the size split in the browser and accept that the reported
    agreement depth is then partly synthetic. See docs/SPECIFICATION.md §5.

name    `tandem_p<period>` so the period is visible on mouseover in every
        track, whether or not --classify-period was used.

hit_id  `TAN_<n>`, a stable within-file serial.

Usage:
    python fastan2bed.py fGobNig-tan.bed -o fastan.bed16.bed
    python fastan2bed.py fGobNig-tan.bed --classify-period | gzip > fastan.bed.gz
"""

import argparse
import gzip
import sys

BED16_HEADER = (
    "#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\t"
    "SW_score\tperc_div\tperc_del\tperc_ins\tquery_left\t"
    "repeat_class_family\trepeat_start\trepeat_end\trepeat_left\thit_id\n"
)

# Conventional tandem size classes. Only used with --classify-period.
MICROSAT_MAX_PERIOD = 6


def _open(path):
    if path == "-":
        return sys.stdin
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def convert(infile, outfile, classify_period=False, emit_divergence=True):
    outfile.write(BED16_HEADER)

    n_in = n_out = 0
    n_malformed = 0      # fewer than 5 fields
    n_empty = 0          # zero-length or inverted interval
    n_bad_number = 0     # non-integer where an integer is required
    periods_seen = set()

    for line in infile:
        line = line.rstrip("\n").rstrip("\r")
        if not line or line.startswith(("#", "track", "browser")):
            continue
        n_in += 1
        f = line.split("\t")
        if len(f) < 5:
            n_malformed += 1
            continue

        chrom = f[0]
        try:
            start, end = int(f[1]), int(f[2])
            period, identity = int(f[3]), int(f[4])
        except ValueError:
            n_bad_number += 1
            continue

        if end <= start:
            n_empty += 1
            continue

        periods_seen.add(period)

        if classify_period:
            label = "Simple_repeat" if period <= MICROSAT_MAX_PERIOD else "Satellite"
        else:
            label = "tandem"

        # Identity is per mille; BED score is 0-1000, so it transfers directly.
        score = max(0, min(1000, identity))
        div = "NA" if not emit_divergence else f"{(1000 - identity) / 10:.1f}"

        n_out += 1
        outfile.write(
            f"{chrom}\t{start}\t{end}\ttandem_p{period}\t{score}\t.\t"
            f"NA\t{div}\tNA\tNA\tNA\t"
            f"{label}\tNA\tNA\tNA\tTAN_{n_out}\n"
        )

    return {
        "records_in": n_in,
        "records_out": n_out,
        "malformed": n_malformed,
        "bad_number": n_bad_number,
        "empty_interval": n_empty,
        "distinct_periods": len(periods_seen),
        "period_min": min(periods_seen) if periods_seen else None,
        "period_max": max(periods_seen) if periods_seen else None,
    }


def main():
    p = argparse.ArgumentParser(
        description="Convert native FasTAN BED to pipeline BED16.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("input", help="FasTAN .bed (plain or .gz); '-' for stdin")
    p.add_argument("-o", "--output", default="-",
                   help="output BED16 (default: stdout)")
    p.add_argument("--classify-period", action="store_true",
                   help="infer Simple_repeat (period<=6) / Satellite (>6) from "
                        "period size. OFF by default, and think before turning "
                        "it on: FasTAN never makes this call, so the summary "
                        "track scores our inference as if it were a tool's own "
                        "vote. Measured on 3 goby chromosomes: 7.84 Mb of "
                        "synthetic class AGREEMENT plus 1.25 Mb of synthetic "
                        "CONFLICT, neither backed by evidence, and nothing "
                        "downstream can detect it. See docs/SPECIFICATION.md "
                        "section 5.")
    p.add_argument("--no-divergence", action="store_true",
                   help="emit NA for perc_div instead of (1000-identity)/10. "
                        "Use if you want the divergence track to carry only "
                        "library-consensus divergence.")
    a = p.parse_args()

    out = sys.stdout if a.output == "-" else open(a.output, "w")
    try:
        with _open(a.input) as fh:
            stats = convert(fh, out, a.classify_period, not a.no_divergence)
    finally:
        if out is not sys.stdout:
            out.close()

    for k, v in stats.items():
        print(f"{k}\t{v}", file=sys.stderr)
    dropped = stats["malformed"] + stats["bad_number"] + stats["empty_interval"]
    if dropped:
        print(f"WARNING: dropped {dropped} record(s); see counts above",
              file=sys.stderr)


if __name__ == "__main__":
    main()
