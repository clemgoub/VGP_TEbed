"""Ingest, validate and harmonize BED16 repeat annotations.

Reads the VGP_TEbed 16-column BED produced by each tool, validates it against
the assembly, resolves every raw class label through the curated vocabulary
(``config/class_map.tsv``), and classifies intra-tool overlap.

Intra-tool overlap policy (per the design brief): overlapping hits from the SAME
tool never add support -- support is counted per distinct tool via a bitmask, so
this is structural, not a rule that can be forgotten. But not all self-overlap
means the same thing, so it is resolved into three cases by geometry and class:

  nested        one hit strictly contains another and classes differ at or above
                order level -> legitimate biology (a SINE inside an LTR element).
                No penalty.
  redundant     high reciprocal overlap, compatible classes -> library redundancy
                (two family models matching the same locus). No penalty.
  selfconflict  substantial overlap with classes incompatible at class level or
                above -> the tool contradicts itself. The tool's classification
                vote at those bases becomes advisory: excluded from consensus but
                still reported, never silently dropped. The repeat CALL is
                unaffected -- the tool still saw something there.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .vocab import ClassMap, ToolSet, split_path

BED16_COLS = [
    "chrom", "chromStart", "chromEnd", "name", "score", "strand",
    "SW_score", "perc_div", "perc_del", "perc_ins", "query_left",
    "repeat_class_family", "repeat_start", "repeat_end", "repeat_left", "hit_id",
]
NUMERIC_MAYBE_NA = ["SW_score", "perc_div", "perc_del", "perc_ins",
                    "query_left", "repeat_start", "repeat_end", "repeat_left"]


def load_chrom_sizes(path) -> dict[str, int]:
    sizes = {}
    for line in open(path):
        if line.strip() and not line.startswith("#"):
            f = line.split()
            sizes[f[0]] = int(f[1])
    return sizes


def load_chrom_alias(path) -> dict[str, str]:
    """Map every known alias to the assembly's primary sequence name.

    The three tools do not agree on which naming authority to use (EDTA leads
    with CATOHO*, RM2/Pantera with OX*), so this is an explicit reconciliation
    step rather than an assumption.
    """
    alias = {}
    p = Path(path)
    if not p.exists():
        return alias
    for line in open(p):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        primary = f[0]
        for a in f:
            if a:
                alias[a] = primary
    return alias


@dataclass
class IngestResult:
    hits: pd.DataFrame
    stats: dict
    unmapped: pd.DataFrame


def read_bed16(path, tool_id: str, sizes: dict, alias: dict | None = None,
               chroms: set | None = None) -> tuple[pd.DataFrame, dict]:
    """Read and validate one tool's BED16. Returns (hits, stats)."""
    # Detect the header rather than assuming one: with a hardcoded header=0 a
    # file written without a header line silently loses its FIRST RECORD, which
    # no downstream check can catch. Column NAMES are irrelevant (positions are
    # authoritative and names are overridden) -- this only decides skip/no-skip.
    with (gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)) as fh:
        first = fh.readline()
    f0 = first.split("\t", 1)[0].strip().lstrip("#")
    has_header = first.startswith("#") or f0 in ("chrom", "chr", "chromosome")

    df = pd.read_csv(
        path, sep="\t", header=0 if has_header else None,
        names=BED16_COLS, comment=None,
        dtype={c: str for c in ["chrom", "name", "strand",
                                "repeat_class_family", "hit_id"]},
        na_values=["NA", "na", "."], keep_default_na=True,
        low_memory=False,
    )
    n_raw = len(df)
    for c in ["chromStart", "chromEnd"] + NUMERIC_MAYBE_NA:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    stats = {"tool_id": tool_id, "n_raw": n_raw}

    # --- validation
    bad_coord = df.chromStart.isna() | df.chromEnd.isna() | (df.chromEnd <= df.chromStart)
    stats["n_bad_coord"] = int(bad_coord.sum())
    df = df[~bad_coord].copy()
    df["chromStart"] = df.chromStart.astype(np.int64)
    df["chromEnd"] = df.chromEnd.astype(np.int64)

    # --- name reconciliation against the assembly
    if alias:
        df["chrom"] = df.chrom.map(lambda c: alias.get(c, c))
    known = df.chrom.isin(sizes)
    stats["n_unknown_seq"] = int((~known).sum())
    stats["unknown_seq_names"] = sorted(df.loc[~known, "chrom"].unique())[:10]
    df = df[known].copy()

    # --- clamp to sequence bounds (a hit may not run off the end)
    lens = df.chrom.map(sizes).astype(np.int64)
    over = df.chromEnd > lens
    stats["n_clamped"] = int(over.sum())
    df.loc[over, "chromEnd"] = lens[over]
    df = df[df.chromEnd > df.chromStart].copy()

    if chroms is not None:
        df = df[df.chrom.isin(chroms)].copy()

    stats["n_hits"] = len(df)
    stats["bp_naive"] = int((df.chromEnd - df.chromStart).sum())
    stats["n_seq"] = df.chrom.nunique()
    stats["n_labels"] = df.repeat_class_family.nunique()
    stats["n_families"] = df.name.nunique()
    stats["n_na_div"] = int(df.perc_div.isna().sum())
    stats["median_len"] = float((df.chromEnd - df.chromStart).median())
    df["tool_id"] = tool_id
    return df, stats


# Tools whose calls come from detecting sequence structure directly rather than
# from matching a consensus library. Labelling these "homology" on the mouseover
# asserts a library comparison that never happened.
_NON_HOMOLOGY_TOOLS = {"fastan": "tandem array structure"}


# EDTA encodes call provenance in its ID column: TE_homo_* are homology-based,
# TE_struc_* and LTRRT_* are structural. This is real information the other
# tools do not provide, so it is preserved rather than discarded.
def call_evidence(hit_id: pd.Series, tool_id: str) -> pd.Series:
    if tool_id in _NON_HOMOLOGY_TOOLS:
        return pd.Series([_NON_HOMOLOGY_TOOLS[tool_id]] * len(hit_id),
                         index=hit_id.index, dtype="category")
    if tool_id != "edta":
        return pd.Series(["homology"] * len(hit_id), index=hit_id.index, dtype="category")
    s = hit_id.fillna("")
    ev = np.where(s.str.startswith("TE_struc"), "structural",
         np.where(s.str.startswith("LTRRT"), "structural", "homology"))
    return pd.Series(ev, index=hit_id.index, dtype="category")


def harmonize(df: pd.DataFrame, cmap: ClassMap) -> pd.DataFrame:
    """Attach canonical class information. Raw labels are always preserved."""
    pairs = df[["repeat_class_family", "tool_id"]].drop_duplicates()
    recs = {}
    for lab, tool in pairs.itertuples(index=False):
        recs[(lab, tool)] = cmap.lookup(lab if isinstance(lab, str) else "", tool)
    key = list(zip(df.repeat_class_family, df.tool_id))
    df = df.copy()
    df["canonical_path"] = [recs[k]["canonical_path"] for k in key]
    df["is_te"] = pd.Categorical([recs[k]["is_te"] for k in key])
    df["map_confidence"] = pd.Categorical([recs[k]["confidence"] for k in key])
    df["uncertain"] = [bool(recs[k].get("uncertain", False)) for k in key]
    df["class_depth"] = [len(split_path(recs[k]["canonical_path"])) for k in key]
    return df


# --------------------------------------------------------------------------
# Intra-tool overlap classification
# --------------------------------------------------------------------------

def _compatible(pa: str, pb: str) -> bool:
    """True if two canonical paths are compatible (one is a prefix of the other)."""
    return pa.startswith(pb) or pb.startswith(pa)


def _divergence_depth(pa: str, pb: str) -> int:
    """Index of the first level at which two canonical paths differ.

    0 = they disagree about being a repeat at all, 1 = TE vs tandem,
    2 = ClassI vs ClassII, 3 = order, 4 = superfamily. Large = mild.
    Returns 99 when compatible.
    """
    a, b = split_path(pa), split_path(pb)
    for i in range(min(len(a), len(b))):
        if a[i] != b[i]:
            return i
    return 99


def classify_self_overlap(df: pd.DataFrame, recip_redundant=0.80,
                          min_conflict_overlap=0.30) -> pd.DataFrame:
    """Label each hit's relationship to overlapping hits from the same tool.

    Adds ``self_rel`` in {none, nested_inner, nested_host, redundant, selfconflict}.
    Only ``selfconflict`` carries a vote penalty downstream.
    """
    out = np.array(["none"] * len(df), dtype=object)
    # Shallowest divergence depth seen for each hit: 99 = no conflict, 4 =
    # superfamily-level (mild), <=3 = order or above (severe). The summariser
    # applies the penalty threshold, so the policy stays tunable rather than
    # baked into the geometry pass.
    cdepth = np.full(len(df), 99, dtype=np.int16)
    idx_all = df.index.to_numpy()
    pos = {ix: i for i, ix in enumerate(idx_all)}

    for (_tool, _chrom), g in df.groupby(["tool_id", "chrom"], observed=True, sort=False):
        if len(g) < 2:
            continue
        g = g.sort_values("chromStart")
        starts = g.chromStart.to_numpy()
        ends = g.chromEnd.to_numpy()
        paths = g.canonical_path.to_numpy()
        gidx = g.index.to_numpy()
        n = len(g)
        # Sweep: compare each hit against later hits that can still overlap.
        max_end = 0
        j_start = 0
        for i in range(n):
            si, ei, pi = starts[i], ends[i], paths[i]
            for j in range(i + 1, n):
                sj = starts[j]
                if sj >= ei:
                    break
                ej, pj = ends[j], paths[j]
                ov = min(ei, ej) - max(si, sj)
                if ov <= 0:
                    continue
                li, lj = ei - si, ej - sj
                contained = (si <= sj and ej <= ei) or (sj <= si and ei <= ej)
                recip = ov / max(1, min(li, lj))
                compat = _compatible(pi, pj)
                pa, pb = pos[gidx[i]], pos[gidx[j]]
                if contained and not compat and min(li, lj) / max(li, lj) < 0.5:
                    # A short element inside a long one with an incompatible
                    # class: this is nesting, which is real biology.
                    inner, host = (pb, pa) if li > lj else (pa, pb)
                    if out[inner] == "none":
                        out[inner] = "nested_inner"
                    if out[host] == "none":
                        out[host] = "nested_host"
                elif recip >= recip_redundant and compat:
                    for p in (pa, pb):
                        if out[p] == "none":
                            out[p] = "redundant"
                elif recip >= min_conflict_overlap and not compat:
                    # Substantial overlap, incompatible classes, neither clearly
                    # nested: the tool is contradicting itself here.
                    d = _divergence_depth(pi, pj)
                    out[pa] = "selfconflict"
                    out[pb] = "selfconflict"
                    cdepth[pa] = min(cdepth[pa], d)
                    cdepth[pb] = min(cdepth[pb], d)
    res = df.copy()
    res["self_rel"] = pd.Categorical(out)
    res["self_conflict_depth"] = cdepth
    return res


def ingest_all(bed_paths: dict[str, str], sizes_path, alias_path=None,
               class_map_path="config/class_map.tsv",
               tools_path="config/tools.tsv", chroms: set | None = None):
    """Full ingest of every tool. Returns (hits, stats_df, unmapped_df, cmap, tools)."""
    sizes = load_chrom_sizes(sizes_path)
    alias = load_chrom_alias(alias_path) if alias_path else {}
    cmap = ClassMap.load(class_map_path)
    tools = ToolSet.load(tools_path)

    frames, stats = [], []
    for tid, path in bed_paths.items():
        df, st = read_bed16(path, tid, sizes, alias, chroms)
        df = harmonize(df, cmap)
        df["evidence"] = call_evidence(df.hit_id, tid)
        df = classify_self_overlap(df)
        st["n_selfconflict"] = int((df.self_rel == "selfconflict").sum())
        st["n_nested_inner"] = int((df.self_rel == "nested_inner").sum())
        st["n_redundant"] = int((df.self_rel == "redundant").sum())
        st["n_uncertain"] = int(df.uncertain.sum())
        st["pct_superfamily"] = round(float((df.class_depth >= 5).mean() * 100), 1)
        st["pct_abstain"] = round(float((df.class_depth <= 1).mean() * 100), 1)
        frames.append(df)
        stats.append(st)
    hits = pd.concat(frames, ignore_index=True)
    return hits, pd.DataFrame(stats), cmap.unmapped_report(), cmap, tools
