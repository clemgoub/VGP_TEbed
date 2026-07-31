"""Synthetic BED16 generator on real assembly coordinates.

Generates spec-compliant BED16 (VGP_TEbed input format) for a set of mock tools
so the whole pipeline can be built and validated before real EDTA / RM2 /
Pantera / fastLTR outputs arrive. Swapping in real files is a path change.

Design: a ground-truth set of insertions is generated FIRST, then each tool is
simulated *observing* that truth with its own error modes (boundary jitter,
misclassification, dropout, fragmentation, library redundancy). Ground truth is
therefore independent of anything the pipeline computes, which is what makes it
usable as a validation target rather than a circular check.

Injected scenarios, each tagged in truth.tsv by `scenario`:
  concordant        all eligible tools agree on extent and class
  jitter            all agree it exists; boundaries differ by 10-500 bp
  unique            exactly one tool calls it
  classconflict     tools agree on extent, disagree at superfamily or above
  unknownlabel      one tool reports Unknown; others classify
  nested            a younger element inserted into an older one
  fragmented        one insertion broken into fragments sharing an ID
  redundant         one tool emits near-identical overlapping hits (library redundancy)
  selfconflict      one tool emits overlapping hits with incompatible classes
  ltronly           an LTR element; only here can fastLTR contribute
  tandem            satellite/simple repeat; only homology tools are eligible
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
import pandas as pd

BED16_COLS = [
    "chrom", "chromStart", "chromEnd", "name", "score", "strand",
    "SW_score", "perc_div", "perc_del", "perc_ins", "query_left",
    "repeat_class_family", "repeat_start", "repeat_end", "repeat_left", "ID",
]

# Element archetypes: (canonical family stem, RM-style class, typical full length,
# length sd, relative abundance, mean divergence, divergence sd)
ARCHETYPES = [
    ("L1-Gn",       "LINE/L1",              6000, 1800, 0.16,  12.0, 7.0),
    ("CR1-Gn",      "LINE/CR1",             4200, 1500, 0.13,  18.0, 8.0),
    ("RTE-Gn",      "LINE/RTE-BovB",        3200, 1100, 0.06,  20.0, 7.0),
    ("Gypsy-Gn",    "LTR/Gypsy",            5600, 1600, 0.12,   9.0, 6.0),
    ("Copia-Gn",    "LTR/Copia",            5100, 1400, 0.05,  11.0, 6.0),
    ("ERV1-Gn",     "LTR/ERV1",             7200, 2000, 0.04,  14.0, 7.0),
    ("tRNA-SINE-Gn","SINE/tRNA",             280,   70, 0.15,  16.0, 8.0),
    ("MIR-Gn",      "SINE/MIR",              230,   60, 0.06,  25.0, 6.0),
    ("hAT-Gn",      "DNA/hAT-Charlie",      2400,  900, 0.08,  15.0, 7.0),
    ("TcMar-Gn",    "DNA/TcMar-Tc1",        1700,  700, 0.07,  17.0, 7.0),
    ("Helitron-Gn", "RC/Helitron",          3800, 1500, 0.03,  19.0, 8.0),
    ("Sat-Gn",      "Satellite",             340,  120, 0.03,   6.0, 4.0),
    ("Simple-Gn",   "Simple_repeat",          60,   30, 0.02,   0.0, 0.0),
]

# Classification confusions a tool may make, keyed by true class. These are the
# realistic ones -- confusions within an order, plus the Gypsy/Copia mix-up that
# actually happens when only a partial pol domain is visible.
CONFUSIONS = {
    "LTR/Gypsy":       ["LTR/Copia", "LTR/Unknown"],
    "LTR/Copia":       ["LTR/Gypsy", "LTR/Unknown"],
    "LTR/ERV1":        ["LTR/ERVK", "LTR/Unknown"],
    "LINE/L1":         ["LINE/Unknown", "LINE/RTE-BovB"],
    "LINE/CR1":        ["LINE/L2", "LINE/Unknown"],
    "LINE/RTE-BovB":   ["LINE/RTE-X", "LINE/Unknown"],
    "SINE/tRNA":       ["SINE/Unknown", "SINE/MIR"],
    "SINE/MIR":        ["SINE/tRNA", "SINE/Unknown"],
    "DNA/hAT-Charlie": ["DNA/hAT-Ac", "DNA/Unknown"],
    "DNA/TcMar-Tc1":   ["DNA/TcMar-Tigger", "DNA/Unknown"],
    "RC/Helitron":     ["DNA/Unknown", "Unknown"],
    "Satellite":       ["Simple_repeat", "Unknown"],
    "Simple_repeat":   ["Low_complexity"],
}

# How each mock tool translates a RepeatMasker-style class into its own dialect.
EDTA_CODES = {
    "LINE/L1": "RIL", "LINE/CR1": "RIJ", "LINE/RTE-BovB": "RIT",
    "LINE/L2": "RIJ", "LINE/RTE-X": "RIT", "LINE/Unknown": "RIX",
    "LTR/Gypsy": "RLG", "LTR/Copia": "RLC", "LTR/ERV1": "RLR",
    "LTR/ERVK": "RLR", "LTR/Unknown": "RLX",
    "SINE/tRNA": "RST", "SINE/MIR": "RST", "SINE/Unknown": "RSX",
    "DNA/hAT-Charlie": "DTA", "DNA/hAT-Ac": "DTA", "DNA/TcMar-Tc1": "DTT",
    "DNA/TcMar-Tigger": "DTT", "DNA/Unknown": "DTX", "RC/Helitron": "DHH",
    "Satellite": "XXX", "Simple_repeat": "XXX", "Low_complexity": "XXX",
    "Unknown": "XXX",
}

FASTLTR_TERMS = {"LTR/Gypsy": "Gypsy", "LTR/Copia": "Copia",
                 "LTR/ERV1": "unknown_LTR", "LTR/ERVK": "unknown_LTR",
                 "LTR/Unknown": "unknown_LTR"}


@dataclass
class ToolProfile:
    """Per-tool observation model."""
    tool_id: str
    detect: float            # P(detect an eligible element)
    jitter_sd: float         # bp, boundary noise sd
    misclass: float          # P(report a confused class)
    unknown_rate: float      # P(report Unknown instead of a class)
    frag_rate: float         # P(split an element into fragments)
    redundant_rate: float    # P(emit an extra near-identical hit)
    selfconflict_rate: float # P(emit an overlapping hit with an incompatible class)
    scope: str
    rm_fields: bool
    dialect: str             # rm | edta | fastltr


PROFILES = {
    "rm2": ToolProfile("rm2", 0.92, 60, 0.10, 0.22, 0.14, 0.06, 0.02,
                       "general_homology", True, "rm"),
    "edta": ToolProfile("edta", 0.74, 180, 0.14, 0.10, 0.05, 0.03, 0.05,
                        "structural_te", False, "edta"),
    "pantera": ToolProfile("pantera", 0.88, 45, 0.07, 0.09, 0.10, 0.08, 0.01,
                           "general_homology", True, "rm"),
    "fastltr": ToolProfile("fastltr", 0.80, 25, 0.06, 0.15, 0.02, 0.01, 0.00,
                           "ltr_only", False, "fastltr"),
}

_SCOPE_OK = {
    "general_homology": lambda cls: True,
    "structural_te": lambda cls: not cls.startswith(("Satellite", "Simple_repeat",
                                                     "Low_complexity")),
    "ltr_only": lambda cls: cls.startswith("LTR/"),
    "tandem": lambda cls: cls.startswith(("Satellite", "Simple_repeat",
                                          "Low_complexity")),
}


def load_chrom_sizes(path, top_n=None, names=None) -> dict[str, int]:
    sizes = {}
    with open(path) as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                f = line.split()
                sizes[f[0]] = int(f[1])
    if names:
        return {k: sizes[k] for k in names}
    if top_n:
        top = sorted(sizes.items(), key=lambda kv: -kv[1])[:top_n]
        return dict(top)
    return sizes


def generate_truth(chrom_sizes: dict[str, int], density=0.28, rng=None,
                   region_frac=1.0) -> pd.DataFrame:
    """Generate the ground-truth insertion set.

    ``density`` is the target fraction of simulated sequence covered by repeats;
    ``region_frac`` restricts generation to the first fraction of each chromosome
    so the prototype stays small while coordinates remain real.
    """
    rng = rng or np.random.default_rng(20260731)
    stems = [a[0] for a in ARCHETYPES]
    classes = {a[0]: a[1] for a in ARCHETYPES}
    lens = {a[0]: (a[2], a[3]) for a in ARCHETYPES}
    divs = {a[0]: (a[5], a[6]) for a in ARCHETYPES}
    probs = np.array([a[4] for a in ARCHETYPES], float)
    probs /= probs.sum()

    rows, eid = [], 0
    for chrom, full_len in chrom_sizes.items():
        span = int(full_len * region_frac)
        target_bp = int(span * density)
        placed_bp, guard = 0, 0
        occupied: list[tuple[int, int]] = []
        while placed_bp < target_bp and guard < 200000:
            guard += 1
            stem = stems[rng.choice(len(stems), p=probs)]
            mu, sd = lens[stem]
            # Truncated: most copies are 5'-truncated fragments, a minority full length.
            full = rng.random() < 0.30
            length = int(max(60, rng.normal(mu, sd) if full else rng.uniform(0.15, 0.85) * mu))
            start = int(rng.integers(0, max(1, span - length)))
            end = start + length
            if any(start < oe and end > os_ for os_, oe in occupied[-60:]):
                continue
            occupied.append((start, end))
            dmu, dsd = divs[stem]
            div = float(np.clip(rng.normal(dmu, dsd), 0.0, 45.0))
            cons_len = int(mu)
            if full:
                cons_start, cons_end = 1, cons_len
            else:  # 5' truncation: match starts partway into the consensus
                cons_start = int(rng.integers(1, max(2, cons_len - length)))
                cons_end = min(cons_len, cons_start + length)
            rows.append(dict(
                elem_id=f"E{eid:06d}", chrom=chrom, start=start, end=end,
                family=f"{stem}_{rng.integers(1, 40)}", true_class=classes[stem],
                strand="+" if rng.random() < 0.5 else "-",
                perc_div=round(div, 1), full_length=full,
                cons_len=cons_len, cons_start=cons_start, cons_end=cons_end,
                scenario="concordant",
            ))
            eid += 1
            placed_bp += length
    truth = pd.DataFrame(rows).sort_values(["chrom", "start"]).reset_index(drop=True)
    return _assign_scenarios(truth, rng)


def _assign_scenarios(truth: pd.DataFrame, rng) -> pd.DataFrame:
    """Tag a subset of elements with the specific scenarios to be exercised."""
    n = len(truth)
    idx = rng.permutation(n)
    scen = np.array(["concordant"] * n, dtype=object)

    def take(pool, k):
        chunk, rest = pool[:k], pool[k:]
        return chunk, rest

    pool = list(idx)
    for name, frac in [("jitter", 0.14), ("unique", 0.10), ("classconflict", 0.09),
                       ("unknownlabel", 0.08), ("fragmented", 0.07),
                       ("redundant", 0.05), ("selfconflict", 0.04)]:
        chunk, pool = take(pool, int(n * frac))
        scen[chunk] = name
    truth["scenario"] = scen
    # Structural scenarios are determined by class, not by the draw.
    truth.loc[truth.true_class.str.startswith("LTR/"), "scenario"] = np.where(
        truth.loc[truth.true_class.str.startswith("LTR/"), "scenario"] == "concordant",
        "ltronly", truth.loc[truth.true_class.str.startswith("LTR/"), "scenario"])
    tand = truth.true_class.str.startswith(("Satellite", "Simple_repeat"))
    truth.loc[tand, "scenario"] = "tandem"
    return truth


def _nest_elements(truth: pd.DataFrame, rng, frac=0.05) -> pd.DataFrame:
    """Insert short elements inside longer ones -- the `nested` scenario."""
    long_e = truth[(truth.end - truth.start) > 2500]
    if long_e.empty:
        return truth
    k = max(1, int(len(long_e) * frac))
    hosts = long_e.sample(min(k, len(long_e)), random_state=int(rng.integers(1e6)))
    rows = []
    for i, (_, h) in enumerate(hosts.iterrows()):
        length = int(rng.integers(180, 420))
        off = int(rng.integers(400, max(401, h.end - h.start - length - 200)))
        s = int(h.start + off)
        rows.append(dict(
            elem_id=f"N{i:05d}", chrom=h.chrom, start=s, end=s + length,
            family=f"tRNA-SINE-Gn_{rng.integers(1, 40)}", true_class="SINE/tRNA",
            strand="+" if rng.random() < 0.5 else "-",
            perc_div=round(float(np.clip(rng.normal(14, 6), 0, 45)), 1),
            full_length=True, cons_len=280, cons_start=1, cons_end=length,
            scenario="nested", host_elem=h.elem_id,
        ))
    nested = pd.DataFrame(rows)
    truth = truth.copy()
    truth.loc[truth.elem_id.isin(hosts.elem_id), "scenario"] = "nested_host"
    out = pd.concat([truth, nested], ignore_index=True)
    return out.sort_values(["chrom", "start"]).reset_index(drop=True)


def _translate(cls: str, dialect: str) -> str:
    if dialect == "edta":
        return EDTA_CODES.get(cls, "XXX")
    if dialect == "fastltr":
        return FASTLTR_TERMS.get(cls, "unknown_LTR")
    return cls


def _incompatible(cls: str, rng) -> str:
    """A class from a different Wicker class -- for the selfconflict scenario."""
    if cls.startswith(("LINE", "SINE", "LTR")):
        return rng.choice(["DNA/hAT-Charlie", "DNA/TcMar-Tc1", "RC/Helitron"])
    return rng.choice(["LINE/L1", "LTR/Gypsy", "SINE/tRNA"])


def observe(truth: pd.DataFrame, profile: ToolProfile, rng) -> pd.DataFrame:
    """Simulate one tool observing the ground truth."""
    scope_ok = _SCOPE_OK[profile.scope]
    out, hit_id = [], 0

    def emit(e, s, t, cls, ident, div=None, frag=None):
        nonlocal hit_id
        s, t = int(max(0, s)), int(t)
        if t <= s:
            return
        div = e.perc_div if div is None else div
        length = t - s
        rec = {
            "chrom": e.chrom, "chromStart": s, "chromEnd": t,
            "name": e.family if profile.dialect != "edta" else f"TE_{abs(hash(e.family)) % 99999:05d}",
            "score": min(1000, int(length * 0.6 + (40 - min(div, 40)) * 8)),
            "strand": e.strand,
            "repeat_class_family": _translate(cls, profile.dialect),
            "ID": ident,
        }
        if profile.rm_fields:
            span = max(1, e.cons_end - e.cons_start)
            cs = int(e.cons_start + (s - e.start) * span / max(1, e.end - e.start))
            ce = int(cs + length * span / max(1, e.end - e.start))
            cs, ce = max(1, cs), max(2, ce)
            rec.update({
                "SW_score": int(length * 2.2 + (40 - min(div, 40)) * 12),
                "perc_div": round(float(div), 1),
                "perc_del": round(float(max(0, rng.normal(1.8, 1.2))), 1),
                "perc_ins": round(float(max(0, rng.normal(1.2, 0.9))), 1),
                "query_left": 0,
                "repeat_start": cs, "repeat_end": ce,
                "repeat_left": max(0, int(e.cons_len) - ce),
            })
        else:
            rec.update({k: "NA" for k in ("SW_score", "perc_div", "perc_del",
                                          "perc_ins", "query_left", "repeat_start",
                                          "repeat_end", "repeat_left")})
        out.append(rec)
        hit_id += 1

    for e in truth.itertuples():
        if not scope_ok(e.true_class):
            continue
        scen = e.scenario
        # Detection
        p = profile.detect
        if scen == "unique":
            # Exactly one tool sees it; pick deterministically from the element id.
            owner = PROFILE_ORDER[int(e.elem_id[1:], 36) % len(PROFILE_ORDER)]
            if profile.tool_id != owner:
                continue
            p = 1.0
        if rng.random() > p:
            continue

        # Class reported
        cls = e.true_class
        if scen == "unknownlabel" and profile.tool_id in ("rm2", "edta"):
            cls = "Unknown"
        elif scen == "classconflict" and profile.tool_id in ("edta", "fastltr"):
            cls = rng.choice(CONFUSIONS.get(e.true_class, ["Unknown"]))
        elif rng.random() < profile.unknown_rate:
            cls = "Unknown"
        elif rng.random() < profile.misclass:
            cls = rng.choice(CONFUSIONS.get(e.true_class, ["Unknown"]))

        # Boundaries
        jit = profile.jitter_sd * (4.0 if scen == "jitter" else 1.0)
        s = e.start + int(rng.normal(0, jit))
        t = e.end + int(rng.normal(0, jit))
        if t - s < 30:
            s, t = e.start, e.end
        div = float(np.clip(rng.normal(e.perc_div, 1.5), 0, 45))
        ident = f"{profile.tool_id}_{hit_id}"

        # Fragmentation: one insertion reported as pieces sharing an ID
        if scen == "fragmented" or rng.random() < profile.frag_rate:
            if t - s > 600:
                nfrag = int(rng.integers(2, 4))
                cuts = sorted(rng.choice(range(s + 150, t - 150), nfrag - 1, replace=False)) \
                    if t - s > 700 else [(s + t) // 2]
                bounds = [s] + list(cuts) + [t]
                shared = f"{profile.tool_id}_{hit_id}"
                for i in range(len(bounds) - 1):
                    gap = int(rng.integers(20, 120))
                    emit(e, bounds[i], bounds[i + 1] - gap, cls, shared, div)
                continue
        emit(e, s, t, cls, ident, div)

        # Library redundancy: a second near-identical hit, same class
        if scen == "redundant" or rng.random() < profile.redundant_rate:
            emit(e, s + int(rng.normal(0, 40)), t + int(rng.normal(0, 40)),
                 cls, f"{profile.tool_id}_{hit_id}", div)

        # Self-conflict: overlapping hit with an incompatible class
        if scen == "selfconflict" or rng.random() < profile.selfconflict_rate:
            ov_s = s + int((t - s) * 0.25)
            ov_t = s + int((t - s) * 0.85)
            emit(e, ov_s, ov_t, _incompatible(e.true_class, rng),
                 f"{profile.tool_id}_{hit_id}", div)

    df = pd.DataFrame(out, columns=BED16_COLS)
    return df.sort_values(["chrom", "chromStart"]).reset_index(drop=True)


PROFILE_ORDER = ["rm2", "edta", "pantera", "fastltr"]


def write_bed16(df: pd.DataFrame, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("#" + "\t".join(BED16_COLS) + "\n")
        df.to_csv(fh, sep="\t", header=False, index=False)
    return path


def simulate_all(chrom_sizes_path, outdir="data/synthetic", top_n=3,
                 region_frac=0.02, density=0.30, seed=20260731):
    """Full simulation: truth set + one BED16 per tool."""
    rng = np.random.default_rng(seed)
    sizes = load_chrom_sizes(chrom_sizes_path, top_n=top_n)
    truth = generate_truth(sizes, density=density, rng=rng, region_frac=region_frac)
    truth = _nest_elements(truth, rng)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    truth.to_csv(outdir / "truth.tsv", sep="\t", index=False)

    beds = {}
    for tid in PROFILE_ORDER:
        df = observe(truth, PROFILES[tid], rng)
        beds[tid] = write_bed16(df, outdir / f"{tid}.bed")
    # Region actually simulated, for honest genome-fraction denominators.
    regions = pd.DataFrame([{"chrom": c, "start": 0, "end": int(l * region_frac)}
                            for c, l in sizes.items()])
    regions.to_csv(outdir / "simulated_regions.bed", sep="\t",
                   header=False, index=False)
    return truth, beds
