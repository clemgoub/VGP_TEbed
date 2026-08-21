#!/usr/bin/env python3
"""Auto-generated two-level repeat annotation report for a VGP_TEbed hub build.

Usage:
    python report/generate_report.py --assembly GCA_951799975.1 [--repo .]

Consumes:
    work/segments.parquet          (required -- the consensus segmentation)
    config/tools.tsv               (required -- tool manifest w/ scope + rm_fields)
    data/<ASSEMBLY>.chrom.sizes    (required -- genome denominator)
    report/data/families.parquet             (optional -- per-family evidence)
    report/data/dfam_shortlist_clusters.tsv  (optional -- cluster-level candidates)

Emits one self-contained HTML (figures embedded as base64 PNG):
    report/<ASSEMBLY>_report.html

Two levels:
  LEVEL 1 (technical): tool agreement/specificity, conflict structure, and the
  Dfam low-hanging-fruit shortlist with explicit per-gate pass/fail.
  LEVEL 2 (biological): tiered genome-wide repeat estimates (union upper bound,
  >=2-tool support, conflict-free consensus headline) with class breakdown and
  a confidence-aware divergence landscape.

Confidence semantics are read from config/tools.tsv (scope, rm_fields), so the
caveat text regenerates itself when the tool set changes.
"""
from __future__ import annotations

import argparse
import base64
import io
import os
import sys
from datetime import date

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ----------------------------------------------------------------------------- style
BASE, MID, SMALL = 9, 8, 7
plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 200, "font.size": BASE,
    "axes.titlesize": BASE, "axes.labelsize": BASE,
    "legend.fontsize": MID, "xtick.labelsize": SMALL, "ytick.labelsize": SMALL,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlelocation": "left", "font.family": "sans-serif",
})

# One hue per top-level repeat kind, threaded through every figure.
KIND_COLORS = {
    "TE": "#4477AA", "tandem": "#CCBB44", "unresolved": "#BBBBBB",
}
ORDER_COLORS = {  # within-TE orders sample the TE hue family + distinct hues
    "LINE": "#4477AA", "SINE": "#88CCEE", "LTR": "#AA3377", "DIRS": "#EE99AA",
    "PLE": "#44AA99", "TIR": "#117733", "Helitron": "#999933",
    "Maverick": "#DDCC77", "Crypton": "#332288", "ClassI*": "#6699CC",
    "ClassII*": "#66AA55", "TE*": "#7F7F7F",
}
TIER_COLOR = "#4477AA"


def fig_to_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def mb(bp) -> float:
    return bp / 1e6


# ----------------------------------------------------------------------------- inputs
def load_inputs(repo: str, assembly: str):
    seg = pd.read_parquet(
        os.path.join(repo, "work", "segments.parquet"),
        columns=["chrom", "chromStart", "chromEnd", "mask", "n_support",
                 "n_eligible", "agree_depth", "conflict_depth", "consensus_id",
                 "consensus_path", "support_frac", "mean_div"])
    seg["len"] = seg.chromEnd - seg.chromStart

    tools = pd.read_csv(os.path.join(repo, "config", "tools.tsv"),
                        sep="\t", comment="#")
    tools = tools[tools.tool_id.notna()]

    sizes = pd.read_csv(os.path.join(repo, "data", f"{assembly}.chrom.sizes"),
                        sep="\t", names=["chrom", "size"])
    genome_bp = int(sizes["size"].sum())

    ddir = os.path.join(repo, "report", "data")
    fam = shortlist = None
    fam_p = os.path.join(ddir, "families.parquet")
    sl_p = os.path.join(ddir, "dfam_shortlist_clusters.tsv")
    if os.path.exists(fam_p):
        fam = pd.read_parquet(fam_p)
    if os.path.exists(sl_p):
        shortlist = pd.read_csv(sl_p, sep="\t")
    return seg, tools, genome_bp, fam, shortlist


def path_kind(path: str) -> str:
    """repeat:TE:ClassI:LINE:Jockey -> TE ; repeat:tandem:simple -> tandem ; repeat -> unresolved"""
    parts = path.split(":")
    return parts[1] if len(parts) > 1 else "unresolved"


def path_order(path: str) -> str:
    """Order-level label for TE paths: LINE/SINE/LTR/TIR/...; ClassI*/ClassII* when order unresolved."""
    p = path.split(":")
    if len(p) >= 4:
        return p[3]
    if len(p) == 3:
        return p[2] + "*"      # ClassI* / ClassII* : class known, order not
    return "TE*"               # repeat:TE only


# ----------------------------------------------------------------------------- level 2 (biological)
def biological_level(seg, tools, genome_bp):
    out = {}
    rep = seg  # segments only exist where mask>0
    ln = rep["len"].to_numpy()

    tiers = [
        ("Union of all tools", np.ones(len(rep), bool),
         "any single tool called a repeat -- upper bound, includes singletons"),
        ("Supported by >=2 tools", (rep.n_support >= 2).to_numpy(),
         "at least two scope-eligible tools independently overlap"),
        ("Consensus (headline)", ((rep.n_support >= 2) & (rep.conflict_depth < 0)).to_numpy(),
         ">=2 tools and no classification conflict at any level"),
        ("Kind-resolved consensus", ((rep.n_support >= 2) & (rep.conflict_depth < 0)
                                     & (rep.agree_depth >= 2)).to_numpy(),
         "consensus loci where tools also agree the kind (TE vs tandem) or deeper"),
    ]
    out["tiers"] = [(name, mb(ln[m].sum()), 100 * ln[m].sum() / genome_bp, desc)
                    for name, m, desc in tiers]

    # class breakdown at headline tier
    head = rep[(rep.n_support >= 2) & (rep.conflict_depth < 0)].copy()
    head["kind"] = head.consensus_path.map(path_kind)
    out["kind_bp"] = head.groupby("kind")["len"].sum().sort_values(ascending=False)

    te = head[head["kind"] == "TE"].copy()
    te["order"] = te.consensus_path.map(path_order)
    out["order_bp"] = te.groupby("order")["len"].sum().sort_values(ascending=False)

    tan = head[head["kind"] == "tandem"].copy()
    tan["sub"] = tan.consensus_path.map(
        lambda p: p.split(":")[2] if len(p.split(":")) > 2 else "tandem*")
    out["tandem_bp"] = tan.groupby("sub")["len"].sum().sort_values(ascending=False)

    # divergence landscape: corroborated bp with true consensus divergence only
    land = head[head.mean_div.notna() & (head["kind"] == "TE")].copy()
    land["order"] = land.consensus_path.map(path_order)
    bins = np.arange(0, 41, 1)
    orders = out["order_bp"].index.tolist()
    mat = np.zeros((len(orders), len(bins) - 1))
    for i, o in enumerate(orders):
        sub = land[land["order"] == o]
        h, _ = np.histogram(sub.mean_div.clip(0, 39.99), bins=bins,
                            weights=sub["len"])
        mat[i] = h / 1e6
    out["landscape"] = (orders, bins, mat)
    out["landscape_bp"] = land["len"].sum()
    out["headline_bp"] = head["len"].sum()
    out["genome_bp"] = genome_bp

    # auto-caveats from tools.tsv semantics
    cav = []
    div_only = tools[tools.rm_fields == "divergence_only"].tool_id.tolist()
    if div_only:
        cav.append(f"<b>{', '.join(div_only)}</b>: <code>perc_div</code> measures "
                   "unit-to-unit homogeneity within an array, not divergence from a "
                   "library consensus. Excluded from the divergence landscape.")
    no_div = tools[tools.rm_fields == "no"].tool_id.tolist()
    if no_div:
        cav.append(f"<b>{', '.join(no_div)}</b>: no divergence of any kind; these "
                   "contribute to existence/classification but never to age estimates.")
    for scope, msg in [("ltr_only", "only detects LTR retrotransposons"),
                       ("tandem", "only detects tandem repeats"),
                       ("structural_te", "detects TEs only, not tandem repeats")]:
        tl = tools[tools.scope == scope].tool_id.tolist()
        if tl:
            cav.append(f"<b>{', '.join(tl)}</b>: {msg}; excluded from the support "
                       "denominator outside that scope (absence of its call is not a dissent).")
    notrun = tools[tools.ran == "no"].tool_id.tolist()
    if notrun:
        cav.append(f"<b>{', '.join(notrun)}</b>: not run on this assembly (empty track).")
    cav.append("The union tier counts every singleton call, including "
               "over-represented k-mer masking; it is an upper bound, not an estimate.")
    cav.append("The headline tier requires two independent tools and zero class "
               "conflict; repeats genuinely detectable by only one method "
               "(e.g. young LTR insertions found structurally) are pushed below it.")
    out["caveats"] = cav
    return out


def fig_tiers(bio) -> str:
    names = [t[0] for t in bio["tiers"]][::-1]
    vals = [t[1] for t in bio["tiers"]][::-1]
    pcts = [t[2] for t in bio["tiers"]][::-1]
    fig, ax = plt.subplots(figsize=(6.4, 2.2))
    shades = ["#B7C8DC", "#8FA9C8", "#6690B4", TIER_COLOR]
    ax.barh(names, vals, color=shades, height=0.62)
    for i, (v, p) in enumerate(zip(vals, pcts)):
        ax.text(v + 6, i, f"{v:,.0f} Mb ({p:.1f}%)", va="center", fontsize=MID)
    ax.set_xlabel("repeat-covered sequence (Mb)")
    ax.set_title("Genome repeat content falls from %.0f%% to %.0f%% as confidence tightens"
                 % (bio["tiers"][0][2], bio["tiers"][2][2]))
    ax.set_xlim(0, max(vals) * 1.28)
    ax.margins(y=0.08)
    return fig_to_b64(fig)


def fig_class_breakdown(bio) -> str:
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 2.9),
                             gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axes[0]
    ob = bio["order_bp"] / 1e6
    colors = [ORDER_COLORS.get(o, "#7F7F7F") for o in ob.index]
    ax.barh(ob.index[::-1], ob.values[::-1], color=colors[::-1], height=0.7)
    for i, v in enumerate(ob.values[::-1]):
        ax.text(v + ob.max() * 0.02, i, f"{v:,.1f}", va="center", fontsize=SMALL)
    ax.set_xlabel("consensus TE sequence (Mb)")
    ax.set_title("TE orders at consensus confidence")
    ax.set_xlim(0, ob.max() * 1.18)

    ax = axes[1]
    kb = bio["kind_bp"] / 1e6
    tb = bio["tandem_bp"] / 1e6
    rows = [("TE", kb.get("TE", 0), KIND_COLORS["TE"])]
    rows += [(f"tandem: {s}", v, KIND_COLORS["tandem"]) for s, v in tb.items()]
    for k in kb.index:
        if k in ("TE", "tandem"):
            continue
        label = "kind unresolved" if k == "unresolved" else k
        rows.append((label, kb[k], KIND_COLORS.get(k, "#7F7F7F")))
    labs = [r[0] for r in rows][::-1]
    vals = [r[1] for r in rows][::-1]
    cols = [r[2] for r in rows][::-1]
    ax.barh(labs, vals, color=cols, height=0.7)
    for i, v in enumerate(vals):
        ax.text(v + max(vals) * 0.02, i, f"{v:,.1f}", va="center", fontsize=SMALL)
    ax.set_xlabel("consensus sequence (Mb)")
    ax.set_title("TE vs tandem partition")
    ax.set_xlim(0, max(vals) * 1.2)
    fig.tight_layout(w_pad=2.0)
    return fig_to_b64(fig)


def fig_landscape(bio) -> str:
    orders, bins, mat = bio["landscape"]
    fig, ax = plt.subplots(figsize=(7.2, 3.0))
    bottom = np.zeros(mat.shape[1])
    x = bins[:-1]
    for i, o in enumerate(orders):
        ax.bar(x, mat[i], bottom=bottom, width=1.0, align="edge",
               color=ORDER_COLORS.get(o, "#7F7F7F"), label=o, linewidth=0)
        bottom += mat[i]
    ax.set_xlabel("divergence from library consensus (%)")
    ax.set_ylabel("Mb per 1% bin")
    ax.set_title("Divergence landscape of consensus-confidence TE sequence "
                 "(%.0f Mb with consensus divergence)" % mb(bio["landscape_bp"]))
    ax.legend(ncol=2, frameon=False, fontsize=SMALL)
    ax.margins(x=0.01)
    return fig_to_b64(fig)


# ----------------------------------------------------------------------------- level 1 (technical)
def technical_level(seg, tools):
    out = {}
    order = tools.tool_id.tolist()
    bit = {t: 1 << i for i, t in enumerate(order)}
    ran = [t for t in order if
           tools.set_index("tool_id").loc[t, "ran"] == "yes"]
    m = seg["mask"].to_numpy()
    ln = seg["len"].to_numpy()
    cov = {t: (m & bit[t]) > 0 for t in ran}

    J = pd.DataFrame(np.eye(len(ran)), index=ran, columns=ran)
    for i, ti in enumerate(ran):
        for tj in ran[i + 1:]:
            inter = ln[cov[ti] & cov[tj]].sum()
            union = ln[cov[ti] | cov[tj]].sum()
            J.loc[ti, tj] = J.loc[tj, ti] = inter / union if union else np.nan
    out["jaccard"] = J

    nsup = seg.n_support.to_numpy()
    confl = (seg.conflict_depth >= 0).to_numpy()
    rows = []
    for t in ran:
        c = cov[t]
        tot = ln[c].sum()
        rows.append(dict(
            tool=t, covered_Mb=mb(tot),
            solo_frac=ln[c & (nsup == 1)].sum() / tot if tot else np.nan,
            conflict_frac=ln[c & confl].sum() / tot if tot else np.nan,
            scope=tools.set_index("tool_id").loc[t, "scope"]))
    out["spec"] = pd.DataFrame(rows)

    # agreement-depth spectrum (bp by agree_depth), and conflict depth
    out["agree_bp"] = seg.groupby("agree_depth")["len"].sum()
    out["conflict_bp"] = seg[seg.conflict_depth >= 0].groupby(
        "conflict_depth")["len"].sum()
    out["support_bp"] = seg.groupby("n_support")["len"].sum()
    return out


def fig_jaccard(tech) -> str:
    J = tech["jaccard"]
    n = len(J)
    fig, ax = plt.subplots(figsize=(0.62 * n + 1.6, 0.62 * n + 1.2))
    im = ax.imshow(J.to_numpy(), cmap="Blues", vmin=0, vmax=1)
    ax.set_xticks(range(n), J.columns, rotation=45, ha="right")
    ax.set_yticks(range(n), J.index)
    for i in range(n):
        for j in range(n):
            v = J.iloc[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=SMALL,
                    color="white" if v > 0.55 else "#333333")
    ax.set_title("Pairwise overlap (length-weighted Jaccard)")
    fig.colorbar(im, ax=ax, shrink=0.75, label="Jaccard")
    return fig_to_b64(fig)


def fig_specificity(tech) -> str:
    sp = tech["spec"].sort_values("covered_Mb", ascending=True)
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 0.32 * len(sp) + 1.4),
                             sharey=True)
    ax = axes[0]
    corro = sp.covered_Mb * (1 - sp.solo_frac)
    solo = sp.covered_Mb * sp.solo_frac
    ax.barh(sp.tool, corro, color="#4477AA", label="corroborated (>=2 tools)")
    ax.barh(sp.tool, solo, left=corro, color="#CCCCCC", label="this tool only")
    ax.set_xlabel("covered sequence (Mb)")
    ax.set_title("Coverage and how much of it is corroborated")
    ax.legend(frameon=False, fontsize=SMALL, loc="lower right")
    ax = axes[1]
    ax.barh(sp.tool, 100 * sp.conflict_frac, color="#AA3377")
    for i, v in enumerate(100 * sp.conflict_frac):
        ax.text(v + 0.4, i, f"{v:.1f}%", va="center", fontsize=SMALL)
    ax.set_xlabel("% of covered bp on class-conflicted segments")
    ax.set_title("Involvement in classification conflicts")
    ax.set_xlim(0, max(100 * sp.conflict_frac) * 1.25)
    fig.tight_layout(w_pad=2.0)
    return fig_to_b64(fig)


def fig_agreement_depth(tech) -> str:
    # agree_depth = number of hierarchy levels agreed:
    # 1 = "repeat" only, 2 = kind (TE vs tandem), 3 = class, 4 = order, 5 = superfamily
    labels = {0: "no agreed level\n(conflict or no vote)", 1: "repeat only",
              2: "agree to kind\n(TE vs tandem)", 3: "agree to class",
              4: "agree to order", 5: "agree to\nsuperfamily"}
    ab = tech["agree_bp"]
    fig, ax = plt.subplots(figsize=(7.0, 2.6))
    idx = sorted(ab.index)
    vals = [mb(ab[i]) for i in idx]
    labs = [labels.get(i, str(i)) for i in idx]
    colors = ["#BBBBBB" if i < 1 else "#4477AA" for i in idx]
    ax.bar(range(len(idx)), vals, color=colors, width=0.7)
    for i, v in enumerate(vals):
        ax.text(i, v + max(vals) * 0.02, f"{v:,.0f}", ha="center", fontsize=SMALL)
    ax.set_xticks(range(len(idx)), labs, fontsize=SMALL)
    ax.set_ylabel("Mb")
    ax.set_title("How deep does multi-tool class agreement go")
    ax.set_ylim(0, max(vals) * 1.15)
    return fig_to_b64(fig)


# ----------------------------------------------------------------------------- Dfam shortlist section
GATES = [
    ("G1_copies", ">=5 near-full-length copies",
     "copies reaching both ends of the tool consensus (5% tolerance) and >=80% of its span"),
    ("G2_msa_depth", "MSA depth >=3 over >=99% of consensus",
     "per-position copy depth from consensus coordinates -- Dfam's seed-alignment requirement"),
    ("G3_class_agree", ">=60% of bp class-corroborated",
     "family footprint lies on segments where independent tools agree to class "
     "level or deeper (agree_depth >= 3)"),
    ("G4_not_tandem", "<=30% tandem-tool overlap",
     "footprint not dominated by FasTAN/TRF/Satellome calls (satellite masquerading as TE)"),
    ("G5_alignable", "median divergence <=25%",
     "copies close enough to the consensus for a clean MSA"),
    ("G6_te_locus", ">=50% of bp on TE-consensus loci",
     "the multi-tool consensus resolves these loci as TE"),
]


def dfam_section(fam, shortlist):
    out = {}
    if fam is None:
        return None
    evaluated = fam[~fam.is_simple & fam.cons_len.notna()]
    funnel = [("TE families with consensus coordinates", len(evaluated))]
    for g, label, _ in GATES:
        funnel.append((label, int((evaluated[g] & ~evaluated.is_simple).sum())))
    funnel.append(("pass ALL gates", int(evaluated.dfam_pass.sum())))
    out["funnel"] = funnel
    out["by_tool"] = evaluated[evaluated.dfam_pass].groupby("tool").size()
    out["rescued"] = fam[fam.rescued_class.notna()]
    out["completeness"] = evaluated[evaluated.dfam_pass].cons_completeness.value_counts()
    out["shortlist"] = shortlist
    return out


def fig_funnel(df) -> str:
    names = [f[0] for f in df["funnel"]]
    vals = [f[1] for f in df["funnel"]]
    fig, ax = plt.subplots(figsize=(6.8, 0.34 * len(names) + 1.0))
    colors = ["#8FA9C8"] * (len(names) - 1) + ["#117733"]
    ax.barh(names[::-1], vals[::-1], color=colors[::-1], height=0.64)
    for i, v in enumerate(vals[::-1]):
        ax.text(v + max(vals) * 0.01, i, f"{v:,}", va="center", fontsize=SMALL)
    ax.set_xlabel("families")
    ax.set_title("Dfam candidate gates (each gate counted independently)")
    ax.set_xlim(0, max(vals) * 1.12)
    return fig_to_b64(fig)


# ----------------------------------------------------------------------------- html
CSS = """
body{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;max-width:1060px;
     margin:2em auto;padding:0 1.2em;color:#1a1a1a;line-height:1.45}
h1{font-size:1.5em;border-bottom:2px solid #4477AA;padding-bottom:.25em}
h2{font-size:1.2em;margin-top:2em;color:#274b6d}
h3{font-size:1.02em;margin-top:1.4em}
table{border-collapse:collapse;font-size:.82em;margin:.8em 0}
th,td{border:1px solid #d8d8d8;padding:.28em .55em;text-align:right}
th{background:#eef2f6}
td:first-child,th:first-child{text-align:left}
img{max-width:100%;height:auto;display:block;margin:.6em 0}
.caveat{background:#fbf7ec;border-left:4px solid #CCBB44;padding:.6em .9em;
        font-size:.86em;margin:.8em 0}
.note{color:#555;font-size:.86em}
.pass{color:#117733;font-weight:600}.fail{color:#AA3377}
code{background:#f2f2f2;padding:0 .25em;border-radius:3px;font-size:.9em}
.tocbox{background:#f4f7fa;padding:.7em 1em;border-radius:6px;font-size:.9em}
"""


def html_table(df, max_rows=25) -> str:
    return df.head(max_rows).to_html(index=False, border=0, escape=True,
                                     float_format=lambda x: f"{x:,.2f}")


def build_html(assembly, bio, tech, dfam, figs, tools) -> str:
    ran = tools[tools.ran == "yes"]
    H = []
    H.append(f"<html><head><meta charset='utf-8'><title>{assembly} repeat "
             f"annotation report</title><style>{CSS}</style></head><body>")
    H.append(f"<h1>Multi-tool repeat annotation report &mdash; {assembly}</h1>")
    H.append(f"<p class='note'>Generated {date.today().isoformat()} from the "
             f"VGP_TEbed consensus segmentation ({len(ran)} tools run). "
             "All Mb figures are length-weighted; support counts respect each "
             "tool's declared detection scope.</p>")
    H.append("<div class='tocbox'><b>Contents</b> &mdash; "
             "<a href='#bio'>1. Biological summary</a> &middot; "
             "<a href='#tech'>2. Technical: tool agreement</a> &middot; "
             "<a href='#dfam'>3. Dfam low-hanging fruit</a> &middot; "
             "<a href='#methods'>Methods &amp; caveats</a></div>")

    # ---- level 2
    H.append("<h2 id='bio'>1. Biological summary &mdash; confidence-tiered</h2>")
    hl = bio["tiers"][2]
    H.append(f"<p>The headline consensus estimate: <b>{hl[1]:,.0f} Mb "
             f"({hl[2]:.1f}% of the {bio['genome_bp']/1e6:,.0f} Mb assembly)</b> is "
             "covered by repeats supported by at least two independent tools with no "
             "classification conflict. The single-tool union reaches "
             f"{bio['tiers'][0][2]:.1f}%; the difference is annotation that exists "
             "on one tool's authority only and is quarantined, not discarded.</p>")
    H.append(f"<img src='data:image/png;base64,{figs['tiers']}'>")
    H.append("<table><tr><th>tier</th><th>Mb</th><th>% genome</th><th>definition</th></tr>")
    for name, v, p, desc in bio["tiers"]:
        H.append(f"<tr><td>{name}</td><td>{v:,.1f}</td><td>{p:.1f}</td>"
                 f"<td style='text-align:left'>{desc}</td></tr>")
    H.append("</table>")
    H.append(f"<img src='data:image/png;base64,{figs['classes']}'>")
    H.append(f"<img src='data:image/png;base64,{figs['landscape']}'>")
    H.append("<div class='caveat'><b>Confidence caveats (auto-generated from the "
             "tool manifest):</b><ul>")
    for c in bio["caveats"]:
        H.append(f"<li>{c}</li>")
    H.append("</ul></div>")

    # ---- level 1
    H.append("<h2 id='tech'>2. Technical &mdash; tool agreement and specificity</h2>")
    H.append(f"<img src='data:image/png;base64,{figs['jaccard']}'>")
    H.append(f"<img src='data:image/png;base64,{figs['spec']}'>")
    H.append(f"<img src='data:image/png;base64,{figs['depth']}'>")

    # ---- dfam
    H.append("<h2 id='dfam'>3. Dfam low-hanging fruit</h2>")
    if dfam is None:
        H.append("<p class='note'>Family-level evidence tables not built for this "
                 "assembly; run the family pipeline to populate this section.</p>")
    else:
        H.append("<p>A family is a deposition candidate when it passes all six "
                 "gates below. <b>Near-full-length is measured against the tool's "
                 "own consensus, which may itself be incomplete</b> &mdash; the "
                 "completeness column flags families whose consensus shows "
                 "inconsistent implied length across hits (<code>len_inconsistent</code>) "
                 "or thin edge support (<code>weak_5p/weak_3p</code>). Candidates are "
                 "grouped into cross-tool clusters (reciprocal &gt;=50% footprint "
                 "overlap); the cluster, not the single best program, is the "
                 "deposition unit &mdash; the seed alignment should be rebuilt from "
                 "the pooled copies of all members (phase 2).</p>")
        H.append("<table><tr><th>gate</th><th>criterion</th></tr>")
        for _, label, desc in GATES:
            H.append(f"<tr><td>{label}</td><td style='text-align:left'>{desc}</td></tr>")
        H.append("</table>")
        H.append(f"<img src='data:image/png;base64,{figs['funnel']}'>")
        bt = dfam["by_tool"]
        H.append("<p>Passing families by source program: " +
                 ", ".join(f"<b>{t}</b>: {n:,}" for t, n in bt.items()) + ".</p>")
        comp = dfam["completeness"]
        n_ok = int(comp.get("ok", 0))
        H.append(f"<p>Consensus completeness among passers: <b>{n_ok:,}</b> clean; " +
                 ", ".join(f"{k}: {v:,}" for k, v in comp.items() if k != "ok") +
                 ". A <code>len_inconsistent</code> flag usually means the library "
                 "entry is a truncated or chimeric consensus &mdash; exactly the "
                 "families where an MSA-rebuilt consensus will improve on the "
                 "library (phase 2 target list).</p>")
        resc = dfam["rescued"]
        H.append(f"<p><b>Unknown rescue:</b> {len(resc):,} families labelled "
                 "Unknown/bare-DNA by their own program sit on loci where the "
                 "cross-tool consensus confidently resolves a TE class. These are "
                 "free classification improvements for the library.</p>")
        if dfam["shortlist"] is not None:
            sl = dfam["shortlist"].copy()
            keep = ["cluster_id", "n_tools", "tools", "pooled_copies",
                    "pooled_full_len", "best_member", "cons_len_range",
                    "best_completeness", "class_consensus", "median_div",
                    "any_rescued"]
            H.append("<h3>Top 25 cluster-level candidates</h3>")
            H.append(html_table(sl[keep], 25))
            H.append(f"<p class='note'>Full list: report/data/"
                     f"dfam_shortlist_clusters.tsv ({len(sl):,} clusters).</p>")

    # ---- methods
    H.append("<h2 id='methods'>Methods &amp; caveats</h2>")
    H.append("<p class='note'>Segmentation: the genome is cut wherever the set of "
             "supporting tools changes; each segment carries a bitmask of "
             "supporting tools, per-tool class votes mapped to a canonical "
             "hierarchy (repeat : kind : class : order : superfamily), an "
             "eligibility-scoped support count, and the deepest agreed consensus "
             "path. Tools vote only within their declared scope; a tandem-only "
             "tool is not a dissenting vote at a TE locus. Divergence is averaged "
             "only over tools whose <code>perc_div</code> is true divergence from "
             "a library consensus. Family evidence: per-copy consensus "
             "coordinates give a per-position depth profile per family; "
             "near-full-length counts use a 5% edge tolerance. Class "
             "corroboration comes from intersecting each family's footprint with "
             "the segmentation. Tool manifest: <code>config/tools.tsv</code>; "
             "class map: <code>config/class_map.tsv</code>.</p>")
    H.append("<table><tr><th>tool</th><th>version</th><th>scope</th>"
             "<th>divergence semantics</th><th>ran</th></tr>")
    for _, r in tools.iterrows():
        H.append(f"<tr><td>{r.tool_id}</td><td>{r.version}</td><td>{r.scope}</td>"
                 f"<td>{r.rm_fields}</td><td>{r.ran}</td></tr>")
    H.append("</table>")
    H.append("</body></html>")
    return "\n".join(H)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--assembly", required=True)
    ap.add_argument("--repo", default=".")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    seg, tools, genome_bp, fam, shortlist = load_inputs(args.repo, args.assembly)
    print(f"[report] {len(seg):,} segments, genome {genome_bp/1e6:,.0f} Mb, "
          f"families: {'yes' if fam is not None else 'no'}", file=sys.stderr)

    bio = biological_level(seg, tools, genome_bp)
    tech = technical_level(seg, tools)
    dfam = dfam_section(fam, shortlist)

    figs = {
        "tiers": fig_tiers(bio),
        "classes": fig_class_breakdown(bio),
        "landscape": fig_landscape(bio),
        "jaccard": fig_jaccard(tech),
        "spec": fig_specificity(tech),
        "depth": fig_agreement_depth(tech),
    }
    if dfam is not None:
        figs["funnel"] = fig_funnel(dfam)

    html = build_html(args.assembly, bio, tech, dfam, figs, tools)
    out = args.out or os.path.join(args.repo, "report",
                                   f"{args.assembly}_report.html")
    with open(out, "w") as fh:
        fh.write(html)
    print(f"[report] wrote {out} ({os.path.getsize(out)/1e6:.1f} MB)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
