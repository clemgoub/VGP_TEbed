#!/usr/bin/env python3
"""Convert a REPET (TEannot) GFF3 to the pipeline's 16-column BED.

REPET's TEannot pipeline maps a de novo consensus library (TEdenovo) onto the
genome and publishes a GFF3 in which each COPY of a TE is a two-level feature
(confirmed with the REPET developers):

    match         one per copy: the start-stop footprint of the (possibly
                  interrupted) insertion
    match_part    the aligned fragment(s) of that copy; >1 when the copy is
                  fragmented by nested insertions or deletions

FEATURE SELECTION. This script emits one BED row per `match_part`, with the
parent match's copy ID in `hit_id` -- the documented purpose of BED16 col 16
("links fragments of one interrupted insertion"). The alternative, emitting
`match`, would annotate the gap interior between fragments as repeat: on
Gobius niger (Gnig_refTEs_redondant_features_merged.gff) the match level sums
310.6 Mb against 280.3 Mb of actually aligned fragments -- 30.8 Mb of gaps
that REPET itself does not claim. The match footprint is recoverable
downstream by grouping on hit_id; the gap bases are not recoverable in the
other direction. Verified on the same file: the union of match_part spans
equals the parent match span exactly for all 1,016,780 copies, so no element
extent is lost by this choice.

DIVERGENCE. Identity is taken from the `Identity` attribute of the
match_part (alignment identity of the fragment against the TEdenovo library
consensus), and emitted as perc_div = 100 - Identity. This is divergence from
a library consensus -- the same quantity RepeatMasker's perc_div estimates --
NOT a within-element figure, so it may feed the divergence track (EDTA
precedent). The match-level `AlignIdentity` is NOT used: it reads 0.00 on
196,332 copies whose fragments carry real identities, i.e. it is a
placeholder in merged features, not a measurement. Fragments without an
`Identity` attribute (11,408 / 1,128,771 observed) get perc_div = NA.
SW_score, perc_del, perc_ins and query_left are not reported by REPET: NA.

CLASSIFICATION. The Wicker code is parsed from `Wcode:<code>` at the start of
the parent match's `TargetDescription` and emitted verbatim in
repeat_class_family -- including `NA` (unclassified consensus, ~48% of copies
on the goby file; maps to bare `repeat` = existence-only abstention) and
compound codes like `RIX|DTX` (a consensus matching references from more than
one superfamily -- a chimeric or ambiguous consensus). Compound codes are
resolved by config/class_map.tsv to the deepest hierarchy level the codes
share (e.g. RIX|RSX -> repeat:TE:ClassI; RIX|DTX -> repeat:TE). Do not split
them here: choosing one code would assert a superfamily REPET itself did not.

SEQUENCE NAMES. The goby file's seqids are ENA FASTA headers with whitespace
stripped ("CATOHO010000001.1Gobiusnigergenome..."). When a seqid does not
look like a bare accession, the leading INSDC-style accession is extracted
(regex ^[A-Z]{2,6}\\d+\\.\\d+); seqids that are already clean pass through
unchanged. Rewrites are counted and reported.

COORDINATES. GFF3 1-based inclusive -> BED 0-based half-open (start-1, end).

hit_id is the short copy id (leading `ms<N>` of the match ID, REPET's own
numbering; verified unique genome-wide on the goby file). If a match ID does
not start with ms<N>_, the full ID is used. name = the target consensus name;
repeat_start/repeat_end = consensus coordinates of the fragment;
repeat_left = TargetLength - repeat_end (TargetLength from the parent match).

The file is read twice (pass 1: match records; pass 2: match_part records) so
no parent-before-child ordering is assumed. Round-trip accounting (matches
seen, parts seen, rows written, orphans, seqid rewrites) is printed to stderr;
orphan match_parts (parent never seen) are an error, not a silent drop.

Usage:
    python scripts/repetgff2bed.py Gnig_refTEs_redondant_features_merged.gff -o repet.bed
"""

import argparse
import gzip
import re
import sys

BED16_HEADER = (
    "#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\t"
    "SW_score\tperc_div\tperc_del\tperc_ins\tquery_left\t"
    "repeat_class_family\trepeat_start\trepeat_end\trepeat_left\thit_id\n"
)

ACC_RE = re.compile(r"^([A-Z]{2,6}\d+\.\d+)")
MS_RE = re.compile(r"^(ms\d+)_")
WCODE_RE = re.compile(r"Wcode:(\S+)")


def opener(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def attrs(field):
    out = {}
    for kv in field.rstrip(";").split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def clean_seqid(seqid, cache, n_rewritten):
    if seqid in cache:
        return cache[seqid], n_rewritten
    m = ACC_RE.match(seqid)
    if m and m.group(1) != seqid:
        cache[seqid] = m.group(1)
        return m.group(1), n_rewritten + 1
    cache[seqid] = seqid
    return seqid, n_rewritten


def parse_target(t):
    """'<name> <start> <end>' -> (name, start, end) with NA fallback."""
    parts = t.split()
    if len(parts) >= 3:
        return parts[0], parts[1], parts[2]
    return (parts[0] if parts else "NA"), "NA", "NA"


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("gff", help="REPET TEannot GFF3 (optionally .gz)")
    ap.add_argument("-o", "--out", default="/dev/stdout")
    args = ap.parse_args()

    # ---- pass 1: match records -> per-copy Wcode, TargetLength, short id
    matches = {}  # full match ID -> (wcode, target_length, short_id)
    n_match = 0
    with opener(args.gff) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9 or f[2] != "match":
                continue
            n_match += 1
            a = attrs(f[8])
            mid = a.get("ID", "")
            wm = WCODE_RE.search(a.get("TargetDescription", ""))
            wcode = wm.group(1) if wm else "NA"
            tlen = a.get("TargetLength", "NA")
            sm = MS_RE.match(mid)
            short = sm.group(1) if sm else mid
            matches[mid] = (wcode, tlen, short)

    # ---- pass 2: match_part records -> BED rows
    n_part = written = orphans = n_no_ident = n_rewritten = 0
    seq_cache = {}
    short_ids = set()
    out = sys.stdout if args.out == "/dev/stdout" else open(args.out, "w")
    with out, opener(args.gff) as fh:
        out.write(BED16_HEADER)
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9 or f[2] != "match_part":
                continue
            n_part += 1
            a = attrs(f[8])
            parent = a.get("Parent", "")
            if parent not in matches:
                orphans += 1
                continue
            wcode, tlen, short = matches[parent]
            short_ids.add(short)
            chrom, n_rewritten = clean_seqid(f[0], seq_cache, n_rewritten)
            start = int(f[3]) - 1
            end = int(f[4])
            tname, tstart, tend = parse_target(a.get("Target", ""))
            ident = a.get("Identity")
            if ident is None:
                n_no_ident += 1
                score, perc_div = "0", "NA"
            else:
                idv = float(ident)
                score = str(max(0, min(1000, round(idv * 10))))
                perc_div = f"{100.0 - idv:.1f}"
            if tlen != "NA" and tend != "NA":
                try:
                    rleft = str(int(tlen) - int(tend))
                except ValueError:
                    rleft = "NA"
            else:
                rleft = "NA"
            out.write(
                f"{chrom}\t{start}\t{end}\t{tname}\t{score}\t{f[6]}\t"
                f"NA\t{perc_div}\tNA\tNA\tNA\t"
                f"{wcode}\t{tstart}\t{tend}\t{rleft}\t{short}\n"
            )
            written += 1

    print(
        f"matches: {n_match}\n"
        f"match_parts: {n_part}\n"
        f"rows written: {written}\n"
        f"orphan parts (parent never seen, DROPPED): {orphans}\n"
        f"parts without Identity (perc_div=NA): {n_no_ident}\n"
        f"seqids rewritten to accession: {n_rewritten} "
        f"(of {len(seq_cache)} distinct)\n"
        f"distinct copy ids (hit_id): {len(short_ids)} "
        f"(unique per copy: {'OK' if len(short_ids) == n_match else 'CHECK'})",
        file=sys.stderr,
    )
    if orphans:
        sys.exit(f"ERROR: {orphans} orphan match_part records")


if __name__ == "__main__":
    main()
