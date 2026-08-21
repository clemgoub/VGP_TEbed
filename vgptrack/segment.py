"""Per-base segmentation of multi-tool repeat annotations.

Rather than doing interval algebra (which becomes intractable with 4.9M
partially-overlapping intervals from tools that disagree about boundaries), this
paints per-base arrays for one sequence at a time and then run-length encodes
them. Every base carries the identity of the tools supporting it, so support is
counted per DISTINCT TOOL by construction: a tool with ten overlapping hits at a
locus sets its bit once and contributes exactly 1 to support, which is the
structural guarantee that intra-tool overlap can never inflate support.

Arrays painted per sequence (n = sequence length):
  mask        uint16 bitmask, one bit per tool -- who called a repeat here
  vote        uint16 bitmask of tools whose CLASS vote counts here
                     (excludes self-conflicting and uncertain calls)
  cls_id      int32  per-tool canonical class id, one array per tool
  div_sum     float32 / div_n uint8  -- for mean divergence per base
  te_mask     uint16 bitmask of tools calling this base a TE (not tandem/gene)

Memory: the largest goby sequence is 76.5 Mb, so a handful of uint16 arrays is
tens of MB. Sequences are processed one at a time and released.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .vocab import LEVELS, SCOPE_ELIGIBILITY, split_path


def _in_scope(tool, path: str) -> bool:
    """Does `path` fall inside `tool`'s declared detection scope?

    Same generous test `ToolSet.eligible` applies to the support denominator --
    at, above, or below any class in the scope -- but usable on a bare Tool, so
    `resolve_agreement` can gate classification votes without a ToolSet in hand.
    Keep the two in step: they answer the same question, one for the
    denominator and one for the class vote.
    """
    for allowed in SCOPE_ELIGIBILITY.get(tool.scope, ("repeat",)):
        if path.startswith(allowed) or allowed.startswith(path):
            return True
    return False

# Sentinel class id meaning "this tool asserted nothing beyond `repeat`".
ABSTAIN = 0


@dataclass
class ClassRegistry:
    """Interns canonical paths to small integer ids, and their level prefixes."""
    paths: list[str]
    index: dict[str, int]
    # For each id, the id of its prefix at each depth (for agreement by depth).
    prefix_at: np.ndarray  # shape (n_paths, max_depth)
    # Number of levels each path actually asserts. A tool whose path is shallower
    # than depth d has NOT dissented at d -- it simply did not resolve that far.
    depth: np.ndarray = None  # shape (n_paths,)

    @classmethod
    def build(cls, canonical_paths) -> "ClassRegistry":
        uniq = ["repeat"] + sorted(set(canonical_paths) - {"repeat", ""})
        index = {p: i for i, p in enumerate(uniq)}
        maxd = max(len(split_path(p)) for p in uniq)
        # Intern every prefix too, so prefix_at can reference real ids.
        allp = set(uniq)
        for p in uniq:
            parts = split_path(p)
            for d in range(1, len(parts) + 1):
                allp.add(":".join(parts[:d]))
        uniq = ["repeat"] + sorted(allp - {"repeat", ""})
        index = {p: i for i, p in enumerate(uniq)}
        prefix_at = np.zeros((len(uniq), maxd), dtype=np.int32)
        depth = np.zeros(len(uniq), dtype=np.int8)
        for p, i in index.items():
            parts = split_path(p)
            depth[i] = len(parts)
            for d in range(maxd):
                take = min(d + 1, len(parts))
                prefix_at[i, d] = index[":".join(parts[:take])] if take else 0
        out = cls(uniq, index, prefix_at)
        out.depth = depth
        return out

    def __len__(self):
        return len(self.paths)

    @property
    def max_depth(self) -> int:
        return self.prefix_at.shape[1]


def _paint(starts, ends, values, arr, mode="or"):
    """Paint values into an array over [start, end) intervals."""
    if mode == "or":
        for s, e, v in zip(starts, ends, values):
            arr[s:e] |= v
    elif mode == "set":
        for s, e, v in zip(starts, ends, values):
            arr[s:e] = v
    elif mode == "add":
        for s, e, v in zip(starts, ends, values):
            arr[s:e] += v
    elif mode == "max":
        for s, e, v in zip(starts, ends, values):
            np.maximum(arr[s:e], v, out=arr[s:e])
    return arr


def segment_sequence(hits_seq: pd.DataFrame, seq_len: int, tools, registry,
                     conflict_depth_threshold=3, drop_artefact=True):
    """Segment one sequence. Returns a DataFrame of run-length encoded segments.

    ``conflict_depth_threshold``: a self-conflict whose divergence depth is at or
    below this (i.e. the tool contradicts itself at order level or above) makes
    that tool's class vote advisory at those bases. Superfamily-level
    self-conflict (depth 4) is treated as mild and keeps its vote.
    """
    n_tools = len(tools)
    # Every per-base mask below is uint16, so bit 16 and above would be silently
    # truncated -- support counts would look plausible and be wrong. ToolSet
    # permits up to 64 tools, so this bound is the narrower one and must be
    # checked here rather than at load time. (Widened from uint8 2026-08-21
    # when TRF and WindowMasker became tools 8 and 9; the goby genome costs
    # ~76.5 MB extra per mask array at uint16, which is nothing.)
    max_bit = max((t.bit for t in tools), default=-1)
    if max_bit >= 16:
        raise ValueError(
            f"tool '{[t.tool_id for t in tools if t.bit == max_bit][0]}' has "
            f"bit {max_bit}, but the per-base masks are uint16 (max bit 15). "
            "Widen mask/vote_mask/te_mask/nonte_mask/struct_mask/eligible_mask "
            "(and the _POPCOUNT table) before adding a 17th tool.")
    mask = np.zeros(seq_len, dtype=np.uint16)
    vote_mask = np.zeros(seq_len, dtype=np.uint16)
    te_mask = np.zeros(seq_len, dtype=np.uint16)
    nonte_mask = np.zeros(seq_len, dtype=np.uint16)
    div_sum = np.zeros(seq_len, dtype=np.float32)
    div_n = np.zeros(seq_len, dtype=np.uint8)
    struct_mask = np.zeros(seq_len, dtype=np.uint16)
    # Per-tool class id at each base. Later hits overwrite earlier ones; to make
    # this deterministic, hits are painted in order of increasing class depth so
    # the most specific classification wins the base.
    # Rows are indexed by tool.bit, NOT by position in `tools`. Those differ
    # whenever the running subset is not a prefix of the manifest -- e.g. a
    # not-yet-run tool sitting above a running one, which is the normal state
    # once tools are appended over time. Sizing this by len(tools) silently
    # worked while the running tools happened to hold the lowest bits and then
    # raised IndexError the moment they did not.
    cls = np.zeros((max((t.bit for t in tools), default=-1) + 1, seq_len),
                   dtype=np.int32)

    for tool in tools:
        g = hits_seq[hits_seq.tool_id == tool.tool_id]
        if g.empty:
            continue
        if drop_artefact:
            g = g[~g.canonical_path.str.startswith("repeat:artefact")]
            if g.empty:
                continue
        bit = np.uint16(1 << tool.bit)
        s = g.chromStart.to_numpy(np.int64)
        e = g.chromEnd.to_numpy(np.int64)
        _paint(s, e, np.full(len(g), bit), mask, "or")

        # Class votes: exclude bases where this tool contradicts itself at order
        # level or above, and exclude '?'-flagged uncertain calls.
        ok = ~(((g.self_rel == "selfconflict") &
                (g.self_conflict_depth <= conflict_depth_threshold)) | g.uncertain)
        gv = g[ok.to_numpy()]
        if not gv.empty:
            _paint(gv.chromStart.to_numpy(np.int64), gv.chromEnd.to_numpy(np.int64),
                   np.full(len(gv), bit), vote_mask, "or")
            # Paint class ids, shallow first so deeper classifications win.
            gs = gv.sort_values("class_depth")
            ids = gs.canonical_path.map(registry.index).to_numpy(np.int32)
            _paint(gs.chromStart.to_numpy(np.int64), gs.chromEnd.to_numpy(np.int64),
                   ids, cls[tool.bit], "set")
            te = gv.is_te.astype(str).to_numpy()
            for sub, target in ((te == "yes", te_mask), (te == "no", nonte_mask)):
                if sub.any():
                    _paint(gv.chromStart.to_numpy(np.int64)[sub],
                           gv.chromEnd.to_numpy(np.int64)[sub],
                           np.full(int(sub.sum()), bit), target, "or")

        # Divergence: mean over tools reporting divergence FROM A CONSENSUS.
        # A tool whose perc_div measures something else (FasTAN: unit-to-unit
        # divergence within a tandem array) is excluded here so the genome-wide
        # track stays a single interpretable quantity. Its divergence is still
        # shown, labelled, on its own per-tool track.
        gd = g[g.perc_div.notna()] if tool.divergence_is_consensus else g.iloc[:0]
        if not gd.empty:
            _paint(gd.chromStart.to_numpy(np.int64), gd.chromEnd.to_numpy(np.int64),
                   gd.perc_div.to_numpy(np.float32), div_sum, "add")
            _paint(gd.chromStart.to_numpy(np.int64), gd.chromEnd.to_numpy(np.int64),
                   np.ones(len(gd), np.uint8), div_n, "add")
        gsx = g[g.evidence.astype(str) == "structural"]
        if not gsx.empty:
            _paint(gsx.chromStart.to_numpy(np.int64), gsx.chromEnd.to_numpy(np.int64),
                   np.full(len(gsx), bit), struct_mask, "or")

    return _encode_runs(mask, vote_mask, cls, te_mask, nonte_mask,
                        div_sum, div_n, struct_mask, tools, registry)


_POPCOUNT = np.array([bin(i).count("1") for i in range(65536)], dtype=np.uint8)


def _encode_runs(mask, vote_mask, cls, te_mask, nonte_mask, div_sum, div_n,
                 struct_mask, tools, registry):
    """Run-length encode the painted arrays into segments.

    A new segment starts whenever ANY of the state arrays changes, so a segment
    is a maximal run of bases with identical tool support, identical per-tool
    classification, and identical evidence type.
    """
    n = len(mask)
    if n == 0 or not mask.any():
        return pd.DataFrame()
    change = np.zeros(n, dtype=bool)
    change[0] = True
    for arr in (mask, vote_mask, te_mask, nonte_mask, struct_mask):
        change[1:] |= arr[1:] != arr[:-1]
    for row in cls:
        change[1:] |= row[1:] != row[:-1]
    starts = np.flatnonzero(change)
    ends = np.append(starts[1:], n)
    keep = mask[starts] != 0
    starts, ends = starts[keep], ends[keep]
    if len(starts) == 0:
        return pd.DataFrame()

    out = {
        "chromStart": starts.astype(np.int64),
        "chromEnd": ends.astype(np.int64),
        "mask": mask[starts],
        "vote_mask": vote_mask[starts],
        "te_mask": te_mask[starts],
        "nonte_mask": nonte_mask[starts],
        "struct_mask": struct_mask[starts],
        "n_support": _POPCOUNT[mask[starts]],
    }
    with np.errstate(invalid="ignore", divide="ignore"):
        dn = div_n[starts].astype(np.float32)
        out["mean_div"] = np.where(dn > 0, div_sum[starts] / np.maximum(dn, 1), np.nan)
    for t in tools:
        out[f"cls_{t.tool_id}"] = cls[t.bit][starts]
    return pd.DataFrame(out)


# --------------------------------------------------------------------------
# Agreement resolution, vectorized over segments
# --------------------------------------------------------------------------

def resolve_agreement(seg: pd.DataFrame, tools, registry) -> pd.DataFrame:
    """Compute consensus class, agreement depth and conflict depth per segment.

    Vectorized: for each depth d, take every voting tool's class prefix at depth
    d and test whether all informative voters share it. The deepest d at which
    they do is the consensus. Abstentions (class id 0 = `repeat`) never break
    agreement -- they are excluded from the comparison at every depth, which is
    what stops RepeatModeler's 39% Unknown from demoting well-supported loci.

    Scope is honoured the same way. A scope-restricted tool votes on
    classification only where the locus falls inside its scope; outside it the
    tool has no standing to classify, exactly as `add_eligibility` gives it no
    place in the support denominator there. Without this, FasTAN's `tandem` at a
    locus three TE tools call LTR is scored as a genuine TE-vs-tandem dispute
    and collapses the consensus to bare `repeat` -- 20.7 Mb of the goby genome
    when FasTAN was added. `tandem` is a real assertion, not an abstention, so
    the abstention rule above cannot absorb it; the scope test is what does.

    Resolution therefore runs twice: once over unrestricted tools to establish
    what kind of locus this is, then again admitting each restricted tool only
    where that provisional answer lies in its scope. Where the unrestricted
    tools say nothing, a restricted tool is admitted unconditionally -- it is
    then the only evidence available.
    """
    n = len(seg)
    if n == 0:
        return seg
    maxd = registry.max_depth
    tool_ids = [t.tool_id for t in tools]
    cls_cols = np.stack([seg[f"cls_{t}"].to_numpy(np.int32) for t in tool_ids])  # (T, n)
    votes = np.stack([(seg.vote_mask.to_numpy() >> t.bit & 1).astype(bool) for t in tools])
    informative = votes & (cls_cols != ABSTAIN)

    prefix = registry.prefix_at   # (n_paths, maxd)
    pdepth = registry.depth       # (n_paths,)

    def _cascade(inf, cons, ad, cd, still):
        _cascade_body(inf, cls_cols, prefix, pdepth, maxd, cons, ad, cd, still)

    def _resolve(inf):
        """Run the depth cascade over the given informative mask."""
        consensus = np.zeros(n, dtype=np.int32)
        agree_depth = np.zeros(n, dtype=np.int8)
        conflict_depth = np.full(n, -1, dtype=np.int8)
        still = inf.sum(axis=0) > 0
        _cascade(inf, consensus, agree_depth, conflict_depth, still)
        return consensus, agree_depth, conflict_depth

    # Pass 1 -- unrestricted tools only, to establish the kind of locus.
    unrestricted = np.array([t.scope == "general_homology" for t in tools])
    if unrestricted.any() and not unrestricted.all():
        prov_inf = informative & unrestricted[:, None]
        prov_consensus, _, _ = _resolve(prov_inf)
        # Where the unrestricted tools said nothing, anchor on the restricted
        # voter making the DEEPEST assertion. Admitting all of them instead
        # lets tools in disjoint scopes cross-veto: EDTA calling CACTA and
        # FasTAN calling `tandem` at the same locus scored as a dispute and
        # collapsed 5.6 Mb of specific TE classification to bare `repeat`,
        # even though a tandem array overlapping a TE is ordinary biology and
        # neither tool can adjudicate the other's question. Depth is the
        # tiebreak because identifying a TE superfamily is a stronger claim
        # than reporting that an array exists; a tool that asserts more
        # specifically anchors, and shallower out-of-scope votes are gated.
        no_prov = prov_inf.sum(axis=0) == 0
        if no_prov.any():
            restricted_inf = informative & ~unrestricted[:, None]
            depths = np.where(restricted_inf, pdepth[cls_cols], -1)
            deepest = np.argmax(depths, axis=0)
            fallback = cls_cols[deepest, np.arange(n)]
            prov_consensus = np.where(no_prov & (depths.max(axis=0) > 0),
                                      fallback, prov_consensus)
        prov_path = np.array(registry.paths, dtype=object)[prov_consensus]
        # Pass 2 -- admit each restricted tool only where its scope allows.
        uniq, inv = np.unique(prov_path.astype(str), return_inverse=True)
        for i, t in enumerate(tools):
            if unrestricted[i]:
                continue
            ok = np.array([_in_scope(t, p) for p in uniq], dtype=bool)[inv]
            informative[i] &= ok

    n_inf = informative.sum(axis=0)
    consensus = np.zeros(n, dtype=np.int32)          # id of consensus path
    agree_depth = np.zeros(n, dtype=np.int8)         # 0 = none beyond `repeat`
    conflict_depth = np.full(n, -1, dtype=np.int8)   # -1 = no conflict
    still = n_inf > 0
    _cascade(informative, consensus, agree_depth, conflict_depth, still)

    seg = seg.copy()
    seg["n_informative"] = n_inf.astype(np.int8)
    seg["consensus_id"] = consensus
    seg["agree_depth"] = agree_depth
    seg["conflict_depth"] = conflict_depth
    return seg


def _cascade_body(informative, cls_cols, prefix, pdepth, maxd,
                  consensus, agree_depth, conflict_depth, still):
    """Descend the class hierarchy, deepening consensus until voters diverge.

    Mutates `consensus`, `agree_depth` and `conflict_depth` in place.
    """
    for d in range(maxd):
        # A tool participates at depth d only if its own path actually asserts
        # that many levels. Pantera saying "ClassII" while RepeatModeler says
        # "ClassII:TIR:hAT" is NOT a conflict at order level -- Pantera simply
        # stopped at class. Treating it as dissent would penalise the coarser
        # vocabulary, and 62% of Pantera's hits stop at class level.
        reaches = informative & (pdepth[cls_cols] > d)
        pref_d = np.where(reaches, prefix[cls_cols, d], -1)
        first = np.max(pref_d, axis=0)
        n_reach = reaches.sum(axis=0)
        # Agreement at d requires at least one voter reaching d, and all who
        # reach d sharing the same prefix.
        same = np.all(np.where(reaches, pref_d == first, True), axis=0)
        deeper_ok = still & same & (n_reach > 0)
        # Conflict is recorded only where voters that BOTH reached d disagree.
        diverge = still & ~same & (conflict_depth < 0) & (n_reach > 1)
        conflict_depth[diverge] = d
        consensus[deeper_ok] = first[deeper_ok]
        agree_depth[deeper_ok] = d + 1
        # Keep descending while anyone still resolves deeper.
        still = deeper_ok


def segment_all(hits: "pd.DataFrame", tools, registry, sizes: dict,
                conflict_depth_threshold: int = 3, progress=None) -> "pd.DataFrame":
    """Segment every sequence in `sizes` that carries at least one hit.

    Sequences are processed largest-first so the memory high-water mark is hit
    early rather than at an unpredictable point; per-sequence frames are
    concatenated once at the end. `resolve_agreement` runs per sequence because
    it is O(segments) and benefits from the smaller working set.
    """
    import pandas as _pd
    by_chrom = dict(list(hits.groupby("chrom", observed=True)))
    tools = list(tools)
    out = []
    for i, (c, L) in enumerate(sorted(sizes.items(), key=lambda kv: -kv[1])):
        g = by_chrom.get(c)
        if g is None or len(g) == 0:
            continue
        s = segment_sequence(g, L, tools, registry,
                             conflict_depth_threshold=conflict_depth_threshold)
        if len(s) == 0:
            continue
        s = resolve_agreement(s, tools, registry)
        s.insert(0, "chrom", c)
        out.append(s)
        # `i` is a 0-based index, so report i+1 -- otherwise a single-sequence
        # run logs "segmented 0 sequences" while producing segments.
        if progress is not None and (i + 1) % 50 == 0:
            progress(f"  segmented {i + 1} sequences")
    if not out:
        raise ValueError("no sequence produced any segment -- check that hit "
                         "sequence names match the chrom.sizes file")
    return _pd.concat(out, ignore_index=True)


def add_eligibility(seg: pd.DataFrame, tools, registry) -> pd.DataFrame:
    """Attach the support DENOMINATOR: how many tools could have called this.

    A tool that cannot detect the consensus class (fastLTR at a LINE locus, EDTA
    at a satellite) is excluded rather than counted as a dissenting vote. This is
    what makes 2/2 at an LTR-only locus read as full agreement instead of 2/4.
    """
    paths = np.array(registry.paths, dtype=object)
    consensus_path = paths[seg.consensus_id.to_numpy()]
    n_elig = np.zeros(len(seg), dtype=np.int8)
    elig_mask = np.zeros(len(seg), dtype=np.uint16)
    # Only a handful of distinct consensus paths exist; compute per unique path.
    uniq, inv = np.unique(consensus_path.astype(str), return_inverse=True)
    ne = np.array([tools.n_eligible(p) for p in uniq], dtype=np.int8)
    em = np.array([tools.eligible_mask(p) for p in uniq], dtype=np.uint16)
    n_elig = ne[inv]
    elig_mask = em[inv]
    seg = seg.copy()
    seg["consensus_path"] = consensus_path
    # Observed evidence overrides declared scope. A tool that ACTUALLY called a
    # repeat here was evidently capable of calling it, whatever the manifest
    # says -- e.g. EDTA is declared structural_te and so nominally ineligible at
    # a satellite locus, but it does emit TE calls overlapping satellites. Union
    # the eligible mask with the observed mask so the denominator can never be
    # smaller than the numerator; without this, support_frac exceeds 1.0 on
    # ~0.34 Mb of tandem sequence.
    elig_mask = (elig_mask | seg["mask"].to_numpy()).astype(np.uint16)
    n_elig = _POPCOUNT[elig_mask]
    seg["n_eligible"] = n_elig
    seg["eligible_mask"] = elig_mask
    # Support fraction uses the eligible denominator, never the tool count.
    seg["support_frac"] = np.where(n_elig > 0,
                                   seg.n_support.to_numpy() / np.maximum(n_elig, 1), 0.0)
    return seg

