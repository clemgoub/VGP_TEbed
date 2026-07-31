"""End-to-end driver: BED16 inputs -> validated UCSC/IGV hub directory.

    python -m vgptrack.cli build \
        --assembly GCA_951799975.1 \
        --sizes data/GCA_951799975.1.chrom.sizes \
        --alias data/GCA_951799975.1.chromAlias.txt \
        --bed rm2=/path/rm2.bed --bed edta=/path/edta.bed \
        --out hub

Every tool named with --bed must exist in config/tools.tsv. Tools listed in the
manifest but NOT given a --bed are treated as "did not run": they are excluded
from support denominators, and an EMPTY bigBed is still written for them so the
set of filenames is identical across every assembly directory (a UCSC contrib
requirement -- one trackDb serves all assemblies).
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

from . import bigfiles, hub, ingest, segment, summary, vocab


def _log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def build(args: argparse.Namespace) -> int:
    t_start = time.time()
    beds = dict(kv.split("=", 1) for kv in args.bed)
    outdir = os.path.join(args.out, args.assembly)
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(args.work, exist_ok=True)

    tools = vocab.ToolSet.load(args.tools)
    palette = vocab.Palette.load(args.palette)

    known = {t.tool_id for t in tools}
    unknown = set(beds) - known
    if unknown:
        sys.exit(f"error: --bed named tools absent from {args.tools}: {sorted(unknown)}\n"
                 f"       add a row for each, or correct the tool_id. Known: {sorted(known)}")

    hits, stats, unmapped, cmap, tools = ingest.ingest_all(
        beds, args.sizes, alias_path=args.alias,
        class_map_path=args.class_map, tools_path=args.tools)

    # The manifest may claim ran=yes; the BED files on THIS command line are the
    # authority for what actually ran in THIS build.
    for t in tools:
        t.ran = t.tool_id in beds
    ran = tools.subset([t.tool_id for t in tools if t.ran])
    if not ran:
        sys.exit("error: no --bed inputs given; nothing to build")
    _log(f"{len(ran)} tool(s) ran: {', '.join(ran.ids)}; "
         f"{len(tools) - len(ran)} placeholder(s)")
    _log(f"ingested {len(hits):,} hits")
    stats.to_csv(os.path.join(args.work, "ingest_stats.tsv"), sep="\t", index=False)
    if len(unmapped):
        unmapped.to_csv(os.path.join(args.work, "unmapped_labels.tsv"),
                        sep="\t", index=False)
        _log(f"WARNING: {len(unmapped)} label(s) unmapped -> "
             f"{args.work}/unmapped_labels.tsv (they fall back to bare 'repeat')")

    sizes = ingest.load_chrom_sizes(args.sizes)

    reg = segment.ClassRegistry.build(hits.canonical_path.unique())
    summary.set_registry(reg)

    seg = segment.segment_all(hits, ran, reg, sizes,
                              conflict_depth_threshold=args.conflict_threshold,
                              progress=_log)
    seg = segment.add_eligibility(seg, ran, reg)
    _log(f"{len(seg):,} segments over {int((seg.chromEnd - seg.chromStart).sum()):,} bp")
    seg.to_parquet(os.path.join(args.work, "segments.parquet"),
                   index=False, compression="zstd")

    elem = summary.merge_for_display(seg, sliver_bp=args.sliver)
    ts, te = summary.core_runs(seg, elem)
    bed = summary.build_summary_bed(elem, ran, reg, palette, thick=(ts, te))
    _log(f"{len(bed):,} display features")

    # --- summary track -----------------------------------------------------
    as_path = os.path.join(args.work, "repeatSummary.as")
    open(as_path, "w").write(summary.SUMMARY_AS)
    bed_path = os.path.join(args.work, "repeatSummary.bed")
    bed.sort_values(["chrom", "chromStart"], kind="stable").to_csv(
        bed_path, sep="\t", header=False, index=False)
    bigfiles.bed_to_bigbed(bed_path, args.sizes,
                           os.path.join(outdir, "repeatSummary.bb"),
                           as_file=as_path,
                           bed_type=bigfiles.bed_type_from_as(summary.SUMMARY_AS, 12),
                           extra_index=["name"])

    covered = int((seg.chromEnd - seg.chromStart).sum())
    info = bigfiles.bigfile_info(os.path.join(outdir, "repeatSummary.bb"))
    got = int(str(info.get("basesCovered", "0")).replace(",", ""))
    if got != covered:
        sys.exit(f"error: repeatSummary.bb covers {got:,} bp but segmentation "
                 f"produced {covered:,} bp -- display merge lost or invented bases")
    _log(f"repeatSummary.bb OK ({got:,} bp)")

    # --- signals -----------------------------------------------------------
    for name, df in summary.build_signals(seg).items():
        bg = os.path.join(args.work, f"{name}.bedGraph")
        df.to_csv(bg, sep="\t", header=False, index=False)
        bigfiles.bedgraph_to_bigwig(bg, args.sizes, os.path.join(outdir, f"{name}.bw"))
    _log("signal tracks written")

    # --- per-tool + discordance -------------------------------------------
    pt_as = os.path.join(args.work, "repeatByTool.as")
    open(pt_as, "w").write(hub.PERTOOL_AS)
    for t in tools:
        out_bb = os.path.join(outdir, f"repeat_{t.tool_id}.bb")
        if t.ran:
            df = hub.build_pertool_bed(hits[hits.tool_id == t.tool_id], t, palette)
        else:
            df = pd.DataFrame(columns=hub.pertool_columns())  # empty placeholder
        p = os.path.join(args.work, f"repeat_{t.tool_id}.bed")
        df.sort_values(["chrom", "chromStart"], kind="stable").to_csv(
            p, sep="\t", header=False, index=False)
        bigfiles.bed_to_bigbed(p, args.sizes, out_bb, as_file=pt_as,
                               bed_type=bigfiles.bed_type_from_as(hub.PERTOOL_AS, 9),
                               extra_index=["name"])
    _log("per-tool tracks written")

    u_as = os.path.join(args.work, "toolUnique.as")
    open(u_as, "w").write(hub.UNIQUE_AS)
    bu = hub.build_tool_unique_bed(seg, ran, palette)
    p = os.path.join(args.work, "toolUnique.bed")
    bu.sort_values(["chrom", "chromStart"], kind="stable").to_csv(
        p, sep="\t", header=False, index=False)
    bigfiles.bed_to_bigbed(p, args.sizes, os.path.join(outdir, "toolUnique.bb"),
                           as_file=u_as,
                           bed_type=bigfiles.bed_type_from_as(hub.UNIQUE_AS, 9),
                           extra_index=["name"])
    _log(f"toolUnique.bb: {len(bu):,} features")

    # --- hub configuration -------------------------------------------------
    hub.write_trackdb(ran, os.path.join(outdir, "trackDb.txt"), palette,
                      tools=list(tools), vocab_version=cmap.version)
    hub.write_docs(outdir, args.assembly, cmap, ran)

    hub.write_hub_files(args.out, args.assembly, args.twobit,
                        email=args.email, description=args.description)
    _log(f"hub written to {outdir}")

    if not args.no_check:
        rc, out = bigfiles.hub_check(os.path.join(args.out, "hub.txt"))
        if rc != 0:
            _log("hubCheck reported problems:\n" + out)
            return 1
        _log("hubCheck clean")

    _log(f"done in {time.time() - t_start:.0f}s")
    return 0



def cmd_session(args) -> int:
    """Write IGV session + genome descriptor for an already-built hub dir."""
    from . import igvsession

    hubdir = os.path.abspath(args.hub)
    if not os.path.isdir(hubdir):
        print(f"error: no such hub directory: {hubdir}", file=sys.stderr)
        return 1
    # The accession is the directory name -- that is the layout rule.
    accession = os.path.basename(hubdir.rstrip(os.sep))

    gpath = args.genome_out or os.path.join(
        os.path.dirname(os.path.abspath(args.out)) or ".",
        f"{accession}.genome.json")
    igvsession.write_genome_json(accession, gpath, name=args.name)
    igvsession.write_session(hubdir, args.out, accession, locus=args.locus)
    print(f"wrote {args.out} and {gpath}")
    print("  IGV: File > Genome > Load Genome from File...  then  File > Open Session...")
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="vgptrack", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("build", help="build a hub assembly directory")
    b.add_argument("--assembly", required=True, help="accession, e.g. GCA_951799975.1")
    b.add_argument("--sizes", required=True, help="chrom.sizes for the assembly")
    b.add_argument("--alias", help="chromAlias.txt (recommended: tools name sequences "
                                   "by different authorities)")
    b.add_argument("--bed", action="append", required=True, metavar="TOOL=PATH",
                   help="repeat to give one BED16 per tool")
    b.add_argument("--out", default="hub", help="hub root (default: hub)")
    b.add_argument("--work", default="work", help="intermediate directory")
    b.add_argument("--email", default="", help="contact address for hub.txt")
    b.add_argument("--description", default="", help="species/assembly description")
    b.add_argument("--twobit", help="path/URL to the assembly .2bit for UCSC")
    b.add_argument("--tools", default="config/tools.tsv")
    b.add_argument("--class-map", default="config/class_map.tsv")
    b.add_argument("--palette", default="config/palette.tsv")
    b.add_argument("--conflict-threshold", type=int, default=3,
                   help="taxonomic depth at or below which intra-tool self-overlap "
                        "counts as a genuine conflict (default 3 = order)")
    b.add_argument("--sliver", type=int, default=20,
                   help="segments this short are absorbed into neighbours (default 20)")
    b.add_argument("--no-check", action="store_true", help="skip hubCheck")
    b.set_defaults(func=build)

    s = sub.add_parser("session", help="write IGV session + genome files for a built hub")
    s.add_argument("--hub", required=True, help="assembly directory, e.g. hub/GCA_951799975.1")
    s.add_argument("--out", default="igv_session.xml", help="session XML path")
    s.add_argument("--genome-out", help="genome descriptor path "
                                        "(default: <accession>.genome.json beside --out)")
    s.add_argument("--locus", default="", help="opening locus, e.g. OX637595.1:3614000-3620400")
    s.add_argument("--name", default="", help="genome display name")
    s.set_defaults(func=cmd_session)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
