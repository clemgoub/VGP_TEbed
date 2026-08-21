"""Build the repeatSummary BED12+ track and its companion bigWig signals.

The summary track answers, for any locus, the questions from the design brief:

  repeat or not          presence of a feature
  how many tools         `score` (support fraction), `n_support` field, and the
                         per-tool support string in the mouseover
  TE or other            consensus class path + colour
  classification agreed? `agree_depth` / `conflict_depth` fields, and the
                         itemRgb switches to the conflict colour when tools
                         disagree at or above class level
  relation to consensus  mean RepeatMasker divergence from the family model

DISPLAY MERGING. Segmentation emits a new segment whenever ANY per-tool state
changes, which yields 4.3M segments -- too fine to browse and mostly invisible
at any practical zoom. For display, adjacent segments are merged when they share
the display key (support mask, consensus class, agreement depth, conflict
depth). Merging is display-only: work/segments.parquet keeps the full
resolution, and the per-base bigWigs are built from the unmerged segments so no
signal is lost.

BED12 GEOMETRY. The previously-undefined BED12 fields are put to work:
  thickStart/thickEnd  the high-confidence CORE -- the contiguous stretch where
                       every eligible tool agrees. Renders as a solid block
                       inside a thinner outline, so a glance separates
                       "all tools agree" from "one tool's guess".
  itemRgb              consensus class colour (or the conflict colour)
  score                support fraction scaled 0-1000, driving greyscale density
"""

from __future__ import annotations

import numpy as np
import pandas as pd

# autoSql for the summary track. Field order must match the BED12+ output.
SUMMARY_AS = """table repeatSummary
"Multi-tool repeat annotation consensus"
    (
    string  chrom;          "Reference sequence chromosome or scaffold"
    uint    chromStart;     "Start position in chromosome"
    uint    chromEnd;       "End position in chromosome"
    string  name;           "Consensus class (short)"
    uint    score;          "Support fraction x1000"
    char[1] strand;         "+ or - or . for unknown"
    uint    thickStart;     "Start of high-confidence core"
    uint    thickEnd;       "End of high-confidence core"
    uint    reserved;       "RGB colour by consensus class"
    int     blockCount;     "Number of blocks"
    int[blockCount] blockSizes;  "Block sizes"
    int[blockCount] chromStarts; "Block start positions"
    string  consensusClass; "Consensus classification (full path)"
    uint    nSupport;       "Number of distinct tools calling a repeat here"
    uint    nEligible;      "Number of tools able to call this class"
    uint    nClassify;      "Number of tools asserting a class beyond bare repeat"
    lstring supportingTools; "Which tools support this call"
    string  agreement;      "Deepest level of classification agreement"
    string  conflict;       "Level at which tools first disagree, if any"
    lstring perToolClass;   "Classification from each tool"
    string  meanDivergence; "Mean divergence from family consensus (%)"
    string  evidence;       "Structural or homology evidence"
    string  flags;          "Quality flags"
    lstring mouseOver;      "Composed hover summary"
    )
"""

LEVEL_NAMES = {0: "none", 1: "repeat", 2: "TE vs tandem", 3: "class",
               4: "order", 5: "superfamily"}
CONFLICT_RGB = "0,0,0"


def merge_for_display(seg: pd.DataFrame, sliver_bp: int = 20) -> pd.DataFrame:
    """Merge segments into display ELEMENTS.

    Keying the merge on the support mask produces 4.1M features with a median
    length of 43 bp, because the dominant boundary cause (2.39M of 3.15M
    adjacent pairs) is tools disagreeing about where an element starts and
    stops. Those slivers are not separate biological features -- they are the
    boundary jitter itself, and the BED12 geometry is a better place to show it
    than the feature count.

    So an element is a maximal run of covered, genomically adjacent bases whose
    consensus classes remain COMPATIBLE (one path a prefix of the other). Tool
    boundary disagreement is absorbed into the element and re-emerges as the
    thick core (see ``_core_intervals``): outline = union of tools, solid core =
    where every eligible tool agrees. A genuine class change still splits the
    element, so distinct neighbouring repeats stay distinct.

    Full resolution is retained in work/segments.parquet and in the per-base
    bigWigs; nothing here is lost, only displayed differently.
    """
    seg = seg.sort_values(["chrom", "chromStart"], kind="stable").reset_index(drop=True)
    adj = np.zeros(len(seg), dtype=bool)
    adj[1:] = ((seg.chromStart.to_numpy()[1:] == seg.chromEnd.to_numpy()[:-1]) &
               (seg.chrom.to_numpy()[1:] == seg.chrom.to_numpy()[:-1]))
    cp = seg.consensus_path.to_numpy()
    compat = np.zeros(len(seg), dtype=bool)
    compat[1:] = np.fromiter(
        (a.startswith(b) or b.startswith(a) for a, b in zip(cp[:-1], cp[1:])),
        dtype=bool, count=len(seg) - 1)
    # Compatibility alone over-merges: a Jockey segment abuts a class-conflicted
    # segment whose consensus has been truncated to `repeat:TE`, which is a
    # prefix of the Jockey path and therefore "compatible" -- so two distinct
    # repeats fuse and the element inherits a conflict neither tool made at the
    # Jockey locus. Conflict state must therefore also match: a conflicted
    # stretch never merges into an unconflicted one.
    cf = seg.conflict_depth.to_numpy()
    same_conflict = np.zeros(len(seg), dtype=bool)
    same_conflict[1:] = (cf[1:] < 0) == (cf[:-1] < 0)
    # Sliver absorption. Where two tools' boundaries differ by a few bases, the
    # overlap edge produces a 1-20 bp segment carrying a spurious "conflict"
    # that exists only because one tool's element ends mid-way through another's.
    # Such slivers are boundary jitter, not biology: they are not allowed to
    # break an element. The threshold is deliberately small -- a real repeat
    # fragment below ~20 bp is not independently meaningful at browser
    # resolution, and the full-resolution segments retain the detail.
    seg_len = (seg.chromEnd - seg.chromStart).to_numpy()
    sliver = seg_len <= sliver_bp
    # A sliver never starts a new group, and never prevents its neighbours
    # from joining across it.
    same_conflict = same_conflict | sliver
    compat = compat | sliver
    new_group = ~(adj & compat & same_conflict)
    gid = np.cumsum(new_group) - 1

    lens = (seg.chromEnd - seg.chromStart).to_numpy()
    seg2 = seg.assign(
        _w=lens,
        _dw=np.where(np.isnan(seg.mean_div.to_numpy()), 0.0,
                     seg.mean_div.to_numpy() * lens),
        _dn=np.where(np.isnan(seg.mean_div.to_numpy()), 0, lens),
        _sw=seg.support_frac.to_numpy() * lens,
        # peak support marks the bases that will become the thick core
        _peak=(seg.support_frac.to_numpy() >= 0.999),
    )
    g = seg2.groupby(gid, sort=False)
    out = pd.DataFrame({
        "chrom": g.chrom.first().to_numpy(),
        "chromStart": g.chromStart.min().to_numpy(),
        "chromEnd": g.chromEnd.max().to_numpy(),
        # Union of every tool touching the element, and of every tool eligible.
        "mask": g["mask"].apply(lambda s: np.bitwise_or.reduce(s.to_numpy())).to_numpy(),
        "vote_mask": g.vote_mask.apply(lambda s: np.bitwise_or.reduce(s.to_numpy())).to_numpy(),
        "struct_mask": g.struct_mask.apply(lambda s: np.bitwise_or.reduce(s.to_numpy())).to_numpy(),
        "eligible_mask": g.eligible_mask.apply(lambda s: np.bitwise_or.reduce(s.to_numpy())).to_numpy(),
        "n_support_max": g.n_support.max().to_numpy(),
        "n_eligible": g.n_eligible.max().to_numpy(),
        # Deepest consensus reached anywhere in the element. Compatibility is
        # the merge condition, so the deepest path subsumes all the others.
        "agree_depth": g.agree_depth.max().to_numpy(),
        # SHALLOWEST conflict anywhere in the element is the honest summary: an
        # element containing one base of Class I / Class II disagreement is
        # class-conflicted. Using max() reported the mildest conflict present,
        # which understated disagreement.
        "conflict_depth": g.conflict_depth.apply(
            lambda s: (lambda v: v[v >= 0].min() if (v >= 0).any() else -1)(s.to_numpy())
        ).to_numpy(),
        "core_bp": g._peak.sum().to_numpy(),
        "_wsum": g._w.sum().to_numpy(),
    })
    idx_deep = g.agree_depth.idxmax().to_numpy()
    out["consensus_id"] = seg.consensus_id.to_numpy()[idx_deep]
    out["consensus_path"] = seg.consensus_path.to_numpy()[idx_deep]

    # An element can span a sub-interval where tools reach superfamily consensus
    # AND another where they diverge at class level. Reporting "superfamily
    # agreement" alongside "class conflict" is self-contradictory, so agreement
    # is capped at the first level of conflict and the consensus path is
    # truncated to match. The element then reads honestly: "three tools agree
    # this is a TE, they disagree about Class I vs Class II."
    cdv = out.conflict_depth.to_numpy()
    capped = np.where(cdv >= 0, np.minimum(out.agree_depth.to_numpy(), cdv),
                      out.agree_depth.to_numpy())
    need = capped < out.agree_depth.to_numpy()
    if need.any():
        cp_out = out.consensus_path.to_numpy().copy()
        trunc = [":".join(p.split(":")[:d]) if d > 0 else "repeat"
                 for p, d in zip(cp_out[need], capped[need])]
        cp_out[need] = trunc
        out["consensus_path"] = cp_out
        out["consensus_id"] = [registry_index_lookup(p) for p in cp_out]
    out["agree_depth"] = capped
    out["n_support"] = out.n_support_max
    dsum = g._dw.sum().to_numpy(); dn = g._dn.sum().to_numpy()
    out["mean_div"] = np.where(dn > 0, dsum / np.maximum(dn, 1), np.nan)
    # Length-weighted support fraction: how much of the element is agreed on.
    out["support_frac"] = g._sw.sum().to_numpy() / np.maximum(out._wsum.to_numpy(), 1)
    out["core_frac"] = out.core_bp / np.maximum(out._wsum, 1)
    # Per-tool class for the mouseover: take each tool's DEEPEST classification
    # anywhere in the element, not its class at the single deepest-consensus
    # segment. A tool that classifies most of an element but abstains on a few
    # bases should be reported by what it actually said.
    for t in [c[4:] for c in seg.columns if c.startswith("cls_")]:
        col = seg[f"cls_{t}"].to_numpy()
        dep = np.where(col > 0, registry_depth_lookup(col), -1)
        tmp = pd.DataFrame({"g": gid, "cls": col, "dep": dep})
        pick = tmp.groupby("g", sort=False).dep.idxmax().to_numpy()
        out[f"cls_{t}"] = col[pick]
    # CORROBORATION CAP, element level. Each segment's consensus already obeys
    # the quorum rule (segment.py _cascade_body), but the merge picks the
    # DEEPEST segment's consensus for the whole element -- and a lone-voter
    # segment legitimately deepens to its tool's full path. An element that
    # mixes edta-only stretches (consensus Helitron, sole voter) with
    # edta+pantera stretches (consensus ClassII, corroborated) would therefore
    # read "Helitron" again: the same overstatement the segment rule fixed,
    # re-introduced one level up (user-reported, same element class). Cap the
    # element's consensus at the deepest level its own per-tool classes
    # corroborate: the depth to which >=2 voting tools' element-level classes
    # follow the consensus path (a lone informative tool keeps full depth --
    # sole assertion is labelled as such downstream).
    reg = _REGISTRY
    if reg is not None:
        cons_ids = out.consensus_id.to_numpy()
        cons_depth = reg.depth[cons_ids]
        maxd = reg.prefix_at.shape[1]
        vm_out = out.vote_mask.to_numpy()
        cls_mat, vote_rows = [], []
        for t, b in _TOOL_BITS.items():
            if f"cls_{t}" in out:
                cls_mat.append(out[f"cls_{t}"].to_numpy())
                vote_rows.append(((vm_out >> b) & 1).astype(bool))
        if cls_mat:
            cls_mat = np.stack(cls_mat)            # (T, n)
            vote_rows = np.stack(vote_rows)        # (T, n)
            informative = vote_rows & (cls_mat > 0)
            n_inf_el = informative.sum(axis=0)
            # common-prefix depth of each tool's class with the consensus
            common = np.zeros(cls_mat.shape, dtype=np.int8)
            for d in range(maxd):
                ok = (informative
                      & (reg.depth[cls_mat] > d) & (cons_depth[None, :] > d)
                      & (reg.prefix_at[cls_mat, d] == reg.prefix_at[cons_ids, d][None, :]))
                common += ok & (common == d)
            # depth backed by >=2 voters: second-largest common depth
            part = np.partition(common, common.shape[0] - 2, axis=0)
            second = part[common.shape[0] - 2] if common.shape[0] >= 2 else common[0]
            corrob = np.where(n_inf_el >= 2, second, cons_depth).astype(np.int8)
            capped2 = np.minimum(out.agree_depth.to_numpy(), corrob)
            need2 = capped2 < out.agree_depth.to_numpy()
            if need2.any():
                cp2 = out.consensus_path.to_numpy().copy()
                cp2[need2] = [":".join(p.split(":")[:d]) if d > 0 else "repeat"
                              for p, d in zip(cp2[need2], capped2[need2])]
                out["consensus_path"] = cp2
                out["consensus_id"] = [registry_index_lookup(p)
                                       for p in out.consensus_path.to_numpy()]
                out["agree_depth"] = capped2

    # Per-tool coverage fraction of the merged element. A tool "supporting" an
    # element may cover 100% of its bases or 20% of them; the mask union that
    # builds `supportingTools` cannot tell those apart, and a user comparing
    # the summary against the per-tool tracks sees the difference immediately
    # (user-reported on a real element: EDTA end-to-end, WindowMasker four
    # slivers totalling 21%, both listed as "2/4 tools" support).
    lens_arr = lens.astype(np.float64)
    m_all = seg["mask"].to_numpy()
    tool_bits = {c[4:]: None for c in seg.columns if c.startswith("cls_")}
    for i_t, t in enumerate(tool_bits):
        covered = ((m_all >> _TOOL_BITS[t]) & 1) * lens_arr
        tmp2 = pd.DataFrame({"g": gid, "c": covered})
        out[f"cov_{t}"] = (tmp2.groupby("g", sort=False).c.sum().to_numpy()
                           / np.maximum(out._wsum.to_numpy(), 1))
    return out.drop(columns=["_wsum"])


# Set by merge_for_display's caller; a module-level hook keeps the registry out
# of the merge signature while still allowing depth-aware per-tool aggregation.
_REGISTRY = None


def registry_depth_lookup(cls_ids):
    if _REGISTRY is None:
        return np.zeros(len(cls_ids), dtype=np.int8)
    return _REGISTRY.depth[cls_ids]


def registry_index_lookup(path: str) -> int:
    if _REGISTRY is None:
        return 0
    return _REGISTRY.index.get(path, 0)


def set_registry(registry):
    global _REGISTRY
    _REGISTRY = registry


# Tool-id -> bit position, set alongside the registry; needed by
# merge_for_display to compute per-tool coverage fractions from the mask.
_TOOL_BITS: dict = {}


def set_tool_bits(bits: dict):
    global _TOOL_BITS
    _TOOL_BITS = dict(bits)


def core_runs(seg: pd.DataFrame, elements: pd.DataFrame, sliver_bp: int = 20) -> tuple[np.ndarray, np.ndarray]:
    """Longest contiguous full-support run inside each element.

    This is the thick core: the stretch where every eligible tool agrees a
    repeat is present. Where no such run exists the core is empty and kent
    renders an unfilled outline.
    """
    seg = seg.sort_values(["chrom", "chromStart"], kind="stable").reset_index(drop=True)
    adj = np.zeros(len(seg), dtype=bool)
    adj[1:] = ((seg.chromStart.to_numpy()[1:] == seg.chromEnd.to_numpy()[:-1]) &
               (seg.chrom.to_numpy()[1:] == seg.chrom.to_numpy()[:-1]))
    cp = seg.consensus_path.to_numpy()
    compat = np.zeros(len(seg), dtype=bool)
    compat[1:] = np.fromiter(
        (a.startswith(b) or b.startswith(a) for a, b in zip(cp[:-1], cp[1:])),
        dtype=bool, count=len(seg) - 1)
    cf = seg.conflict_depth.to_numpy()
    same_conflict = np.zeros(len(seg), dtype=bool)
    same_conflict[1:] = (cf[1:] < 0) == (cf[:-1] < 0)
    sliver = (seg.chromEnd - seg.chromStart).to_numpy() <= sliver_bp
    gid = np.cumsum(~(adj & (compat | sliver) & (same_conflict | sliver))) - 1
    full = seg.support_frac.to_numpy() >= 0.999
    s = seg.chromStart.to_numpy(); e = seg.chromEnd.to_numpy()
    n_el = len(elements)
    best_s = np.zeros(n_el, dtype=np.int64); best_e = np.zeros(n_el, dtype=np.int64)
    cur_s = cur_e = -1; cur_g = -1
    for i in range(len(seg)):
        gi = gid[i]
        if not full[i]:
            cur_s = cur_e = -1
            continue
        if cur_e == s[i] and cur_g == gi:
            cur_e = e[i]
        else:
            cur_s, cur_e, cur_g = s[i], e[i], gi
        if cur_e - cur_s > best_e[gi] - best_s[gi]:
            best_s[gi], best_e[gi] = cur_s, cur_e
    # Empty core -> thickStart == thickEnd == chromStart (unfilled outline)
    empty = best_e <= best_s
    best_s[empty] = elements.chromStart.to_numpy()[empty]
    best_e[empty] = elements.chromStart.to_numpy()[empty]
    return best_s, best_e


def build_summary_bed(seg: pd.DataFrame, tools, registry, palette,
                      min_len: int = 1, thick=None) -> pd.DataFrame:
    """Produce the BED12+11 summary table, ready for bedToBigBed.

    ``seg`` is the ELEMENT table from :func:`merge_for_display`; ``thick`` is the
    ``(start, end)`` pair from :func:`core_runs`.
    """
    keep = (seg.chromEnd - seg.chromStart) >= min_len
    if thick is not None:
        thick = (thick[0][keep.to_numpy()], thick[1][keep.to_numpy()])
    seg = seg[keep].copy()
    paths = np.array(registry.paths, dtype=object)
    tool_ids = [t.tool_id for t in tools]
    bits = {t.tool_id: t.bit for t in tools}

    cons = seg.consensus_path.to_numpy()
    # Colour: conflict at class level or above overrides the class colour, so a
    # Class I / Class II disagreement is visible without opening the mouseover.
    severe = (seg.conflict_depth.to_numpy() >= 0) & (seg.conflict_depth.to_numpy() <= 2)
    rgb = np.array([palette.color_for(p) for p in cons], dtype=object)
    rgb[severe] = CONFLICT_RGB

    short = np.array([palette.label_for(p) for p in cons], dtype=object)
    ad = seg.agree_depth.to_numpy()
    cd = seg.conflict_depth.to_numpy()
    # A feature where every tool classified but they disagree is NOT
    # "unclassified" -- that word means no tool ventured a class. Naming the
    # disputed case correctly is the difference between "nobody knows" and
    # "the tools contradict each other", which is the more interesting result.
    disputed = np.array([s.startswith("Repeat (unclassified)") for s in short]) & (cd >= 0)
    short = np.where(disputed, "Class disputed", short)

    def mask_names(m):
        return ",".join(t for t in tool_ids if m >> bits[t] & 1) or "none"

    umask, uinv = np.unique(seg["mask"].to_numpy(), return_inverse=True)
    mask_str = np.array([mask_names(m) for m in umask], dtype=object)[uinv]
    uelig, einv = np.unique(seg.eligible_mask.to_numpy(), return_inverse=True)
    elig_str = np.array([mask_names(m) for m in uelig], dtype=object)[einv]

    # Per-tool classification string for the mouseover: the whole point of the
    # summary track is that hovering tells you WHY the consensus is what it is.
    # Alongside it, count the tools whose classification BACKS the displayed
    # class: the name field reports {class} {n_backing}/{n_detecting}, so
    # "DNA 3/3" can no longer mean "one tool said DNA and two said Unknown" --
    # and "LINE 5/5" cannot count a tool that called the locus a satellite.
    # Compatibility is token-wise path prefix in either direction (a tool
    # asserting LINE:R2 backs a consensus of LINE; plain string startswith
    # would let ClassII match a ClassI consensus). When the consensus is bare
    # "repeat" (unclassified or disputed), the count falls back to every
    # non-advisory classifying tool, so "Class disputed 2/3" reads as
    # "two classifiers dispute, out of three detectors".
    def _backs(tool_path, cons_path):
        if not tool_path or tool_path == "repeat":
            return False
        if cons_path in ("repeat", ""):
            return True   # unclassified/disputed: count any real classification
        return (tool_path == cons_path
                or tool_path.startswith(cons_path + ":")
                or cons_path.startswith(tool_path + ":"))

    per_tool = []
    n_classify = np.zeros(len(seg), dtype=np.int8)
    classifiers = []
    cls_arrays = {t: seg[f"cls_{t}"].to_numpy() for t in tool_ids if f"cls_{t}" in seg}
    m_arr = seg["mask"].to_numpy()
    v_arr = seg.vote_mask.to_numpy()
    for i in range(len(seg)):
        parts = []
        cls_tools = []
        for t in tool_ids:
            if not (m_arr[i] >> bits[t] & 1):
                continue
            cid = cls_arrays.get(t, np.zeros(1, dtype=int))[i] if t in cls_arrays else 0
            lab = paths[cid] if cid else "no class"
            if not (v_arr[i] >> bits[t] & 1):
                lab += " (advisory)"
            elif cid and _backs(paths[cid], cons[i]):
                cls_tools.append(t)
            parts.append(f"{t}={lab}")
        per_tool.append("; ".join(parts) if parts else "none")
        n_classify[i] = len(cls_tools)
        classifiers.append(cls_tools)

    div = seg.mean_div.to_numpy()
    div_str = np.where(np.isnan(div), "not reported",
                       np.char.add(np.round(div, 1).astype(str), "%"))

    ev = np.where(seg.struct_mask.to_numpy() > 0, "structural + homology", "homology")

    # The denominator is "restricted" only relative to the tools that actually
    # RAN. fastLTR has not run, so comparing against len(tool_ids) flagged
    # essentially every feature and made the flag meaningless.
    n_ran = sum(1 for t in tools if getattr(t, "ran", True))
    ns = seg.n_support.to_numpy(); ne = seg.n_eligible.to_numpy()
    flags = []
    for i in range(len(seg)):
        f = []
        if ns[i] == 1:
            f.append("single-tool")
        if 0 <= cd[i] <= 2:
            f.append("class conflict")
        if ad[i] <= 1:
            f.append("unclassified")
        if ne[i] < n_ran:
            f.append(f"only {ne[i]} of {n_ran} tools can call this class")
        flags.append(",".join(f) if f else "none")

    if thick is None:
        ts = te = seg.chromStart.to_numpy()
    else:
        ts, te = thick
    n = len(seg)
    lens = seg.chromEnd.to_numpy() - seg.chromStart.to_numpy()
    bed = pd.DataFrame({
        "chrom": seg.chrom.to_numpy(),
        "chromStart": seg.chromStart.to_numpy(),
        "chromEnd": seg.chromEnd.to_numpy(),
        # {class} {n_classifying}/{n_detecting}: the numerator counts tools that
        # asserted a class (beyond bare "repeat"), the denominator tools that
        # detected the repeat at all. Existence support n/eligible stays in the
        # mouseover and the repeatSupport track. Previously the name showed
        # n_support/n_eligible, so "DNA 3/3" could mean one DNA vote plus two
        # Unknowns -- the natural reading of the label was false.
        "name": [f"{s} {int(a)}/{int(b)}" for s, a, b in
                 zip(short, n_classify, seg.n_support.to_numpy())],
        "score": np.clip((seg.support_frac.to_numpy() * 1000).round(), 0, 1000).astype(int),
        "strand": ".",
        "thickStart": ts,
        "thickEnd": te,
        "itemRgb": rgb,
        "blockCount": 1,
        "blockSizes": lens,
        "chromStarts": 0,
        "consensusClass": cons,
        "nSupport": seg.n_support.to_numpy(),
        "nEligible": seg.n_eligible.to_numpy(),
        "nClassify": n_classify,
        "supportingTools": mask_str,
        "agreement": [LEVEL_NAMES.get(int(x), str(x)) for x in ad],
        "conflict": ["none" if x < 0 else LEVEL_NAMES.get(int(x) + 1, str(x)) for x in cd],
        "perToolClass": per_tool,
        "meanDivergence": div_str,
        "evidence": ev,
        "flags": flags,
    })
    # ---- composed mouseOver -------------------------------------------------
    # UCSC shows ONE line on hover. The ordering answers the user's questions in
    # priority order: what is it, how many tools agree it is there, how deeply
    # they agree on what it is, and whether anything is wrong with the call.
    lab = [palette.label_for(c) for c in seg.consensus_path.to_numpy()]
    ns = seg.n_support.to_numpy(); ne = seg.n_eligible.to_numpy()
    agr = bed.agreement.to_numpy(); cfl = bed.conflict.to_numpy()
    st = bed.supportingTools.to_numpy(); dv = bed.meanDivergence.to_numpy()
    core = (bed.thickEnd.to_numpy() - bed.thickStart.to_numpy())
    span = (bed.chromEnd.to_numpy() - bed.chromStart.to_numpy())
    # Per-tool coverage of the element, so "supports" cannot read as
    # "spans": a tool covering 21% of the element in slivers and a tool
    # covering 100% both appear in supportingTools, and the difference is
    # exactly what a user comparing against the per-tool tracks sees.
    # Annotate each tool with its coverage %, omitting ">95%" as the
    # uninteresting common case: (edta, windowmasker 21%).
    cov_cols = {c[4:]: seg[c].to_numpy() for c in seg.columns if c.startswith("cov_")}
    def _tools_with_cov(i, names):
        out = []
        for t in names.split(","):
            f = cov_cols.get(t, None)
            v = float(f[i]) if f is not None else 1.0
            out.append(t if v > 0.95 else f"{t} {max(1, round(100 * v))}%")
        return ",".join(out)
    mouse = []
    for i in range(len(seg)):
        # "unclassified" means no tool said anything; when tools DID classify
        # but disagree, say so -- the two are very different states.
        name_i = lab[i]
        if name_i.startswith("Repeat (unclassified)") and cfl[i] != "none":
            name_i = "Class disputed"
        st_i = _tools_with_cov(i, st[i]) if st[i] != "none" else st[i]
        parts = [f"{name_i} | {ns[i]}/{ne[i]} tools ({st_i})"]
        # who actually asserted the class, when fewer than everyone did; and
        # never say "agree" on the word of a single classifier
        if 0 < n_classify[i] < ns[i]:
            parts.append(f"classified by {','.join(classifiers[i])} only")
        if n_classify[i] >= 2:
            parts.append(f"agree to {agr[i]}" if agr[i] != "none"
                         else "no classification agreement")
        elif n_classify[i] == 1:
            parts.append(f"sole assertion at {agr[i]}" if agr[i] != "none"
                         else "sole assertion")
        else:
            parts.append("unclassified by all")
        if cfl[i] != "none":
            parts.append(f"CONFLICT at {cfl[i]}")
        if core[i] > 0:
            parts.append(f"core {core[i]}/{span[i]} bp")
        else:
            parts.append("no full-support core")
        if dv[i] != "not reported":
            parts.append(f"div {dv[i]}")
        mouse.append(" | ".join(parts))
    bed["mouseOver"] = mouse

    return bed


def build_signals(seg: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """bedGraph tables for the companion bigWig signal tracks.

    Built from UNMERGED segments so per-base resolution is preserved.
    """
    base = seg[["chrom", "chromStart", "chromEnd"]]
    out = {}
    out["repeatSupport"] = base.assign(value=seg.n_support.to_numpy().astype(float))
    out["repeatSupportFrac"] = base.assign(value=seg.support_frac.to_numpy().round(4))
    d = seg.mean_div.to_numpy()
    ok = ~np.isnan(d)
    out["repeatDivergence"] = base[ok].assign(value=np.round(d[ok], 2))

    return out
