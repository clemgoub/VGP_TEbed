"""Tool manifest, class vocabulary, palette, and the agreement calculator.

The canonical class vocabulary is treated as replaceable INPUT. Nothing in this
module hardcodes a class name: every mapping decision comes from
``config/class_map.tsv``, whose ``#vocabulary_version`` header propagates into
the hub documentation so a rendered track always declares which vocabulary
produced it. When the unified Dfam/Repbase classification is released, swapping
the file and re-running ``harmonize`` + ``summarize`` is sufficient -- the
per-base segmentation is class-agnostic and does not need recomputing.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

# Ordered names of the levels in a canonical path. Index = depth.
LEVELS = ["repeat", "kind", "class", "order", "superfamily"]

# Labels that assert "something is here" but abstain from classification.
# These are UNINFORMATIVE, never CONFLICTING -- a tool reporting Unknown must not
# demote a locus that three other tools agree is LTR/Gypsy.
UNINFORMATIVE_PATHS = {"repeat", ""}


def _read_tsv(path) -> tuple[pd.DataFrame, dict]:
    """Read a #-commented TSV, returning the frame and any #key: value headers."""
    meta = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            body = line[1:].strip()
            if ":" in body and not body.startswith(" "):
                k, _, v = body.partition(":")
                if k and " " not in k.strip():
                    meta[k.strip()] = v.strip()
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str).fillna("")
    df.columns = [c.strip() for c in df.columns]
    return df, meta


@dataclass
class Tool:
    tool_id: str
    short_label: str
    long_label: str
    version: str
    scope: str
    rm_fields: bool
    library: str
    color: str
    priority: int
    ran: bool
    notes: str = ""
    bit: int = 0  # assigned by ToolSet
    rm_fields_raw: str = ""

    @property
    def divergence_is_consensus(self) -> bool:
        """Is this tool's perc_div divergence FROM A LIBRARY CONSENSUS?

        The genome-wide repeatDivergence track averages this across tools, so
        every contributor must be measuring the same thing. RepeatMasker-style
        divergence (copy vs consensus) is an age proxy. FasTAN's is unit-to-unit
        divergence within a tandem array -- an array-homogeneity measure that is
        not comparable and would make the mean uninterpretable.

        `mixed` counts as yes: EDTA's homology calls carry genuine consensus
        divergence and its structural calls carry NA, so nothing incomparable
        enters the mean. Tools that are not consensus-based still show their
        divergence on their own per-tool track, where it is labelled.
        """
        return self.rm_fields_raw.strip().lower() in ("yes", "mixed")


# Which canonical classes each scope is capable of calling. A tool only enters
# the support denominator for a locus whose consensus class it could have found.
SCOPE_ELIGIBILITY = {
    "general_homology": ("repeat",),          # anything
    "ltr_only": ("repeat:TE:ClassI:LTR", "repeat:TE:ClassI:DIRS"),
    "structural_te": ("repeat:TE",),
    "tandem": ("repeat:tandem",),
}


@dataclass
class ToolSet:
    tools: dict[str, Tool] = field(default_factory=dict)

    @classmethod
    def load(cls, path="config/tools.tsv") -> "ToolSet":
        df, _ = _read_tsv(path)
        ts = cls()
        for bit, (_, r) in enumerate(df.iterrows()):
            if bit >= 64:
                raise ValueError("more than 64 tools; widen the bitmask dtype")
            ts.tools[r["tool_id"]] = Tool(
                tool_id=r["tool_id"], short_label=r["short_label"],
                long_label=r["long_label"], version=r["version"], scope=r["scope"],
                rm_fields=r["rm_fields"].lower() == "yes",
                rm_fields_raw=r["rm_fields"], library=r["library"],
                color=r["color"], priority=int(r["priority"]),
                ran=r["ran"].lower() == "yes", notes=r.get("notes", ""), bit=bit,
            )
        return ts

    def __getitem__(self, k):
        return self.tools[k]

    def __iter__(self):
        return iter(self.tools.values())

    def __len__(self):
        return len(self.tools)

    @property
    def ids(self) -> list[str]:
        return list(self.tools)

    def subset(self, tool_ids) -> "ToolSet":
        """A ToolSet restricted to *tool_ids*, PRESERVING each tool's bit.

        Bit positions must not be renumbered: masks computed against the full
        manifest elsewhere in the pipeline would silently point at the wrong
        tool. The subset therefore has stable bits with gaps.
        """
        keep = set(tool_ids)
        return ToolSet({k: v for k, v in self.tools.items() if k in keep})

    def mask_to_ids(self, mask: int) -> list[str]:
        return [t.tool_id for t in self if mask >> t.bit & 1]

    def eligible(self, tool: Tool, path: str) -> bool:
        """Could *tool* have called a locus of class *path*?

        Eligibility is deliberately generous: a tool counts as eligible if the
        consensus class is at, above, or below any class in its scope. Above,
        because a locus classified only as `repeat` might be anything; below,
        because an ltr_only tool is eligible for a locus called LTR:Gypsy.
        """
        if not tool.ran:
            return False
        for allowed in SCOPE_ELIGIBILITY.get(tool.scope, ("repeat",)):
            if path.startswith(allowed) or allowed.startswith(path):
                return True
        return False

    def n_eligible(self, path: str) -> int:
        return sum(1 for t in self if self.eligible(t, path))

    def eligible_mask(self, path: str) -> int:
        m = 0
        for t in self:
            if self.eligible(t, path):
                m |= 1 << t.bit
        return m


@dataclass
class ClassMap:
    exact: dict[tuple[str, str], dict] = field(default_factory=dict)
    prefix: list[tuple[str, str, dict]] = field(default_factory=list)
    version: str = "unknown"
    meta: dict = field(default_factory=dict)
    unmapped: dict = field(default_factory=dict)

    @classmethod
    def load(cls, path="config/class_map.tsv") -> "ClassMap":
        df, meta = _read_tsv(path)
        cm = cls(version=meta.get("vocabulary_version", "unknown"), meta=meta)
        for _, r in df.iterrows():
            rec = {
                "canonical_path": r["canonical_path"].strip(),
                "canonical_v2": r.get("canonical_v2", "").strip(),
                "is_te": r["is_te"].strip(),
                "confidence": r["confidence"].strip(),
                "raw_rule": r["raw_label"],
            }
            label, tool = r["raw_label"].strip(), r["tool_id"].strip()
            if label.endswith("*"):
                cm.prefix.append((label[:-1].lower(), tool, rec))
            else:
                cm.exact[(label.lower(), tool)] = rec
        # Longest prefix first, so a specific rule always beats a general one.
        cm.prefix.sort(key=lambda x: -len(x[0]))
        return cm

    @staticmethod
    def _normalize(label: str) -> tuple[str, bool]:
        """Strip decoration that carries no class information.

        Returns ``(clean_label, uncertain)``. A trailing ``?`` is RepeatMasker's
        uncertainty marker (e.g. ``SINE?``, ``DNA/hAT-Ac?``): the class is a best
        guess, not a confident call. That is semantically distinct from both a
        confident call and ``Unknown``, so we keep the class -- discarding it
        would throw away real information -- but flag it. Downstream the flag
        demotes the mapping confidence, and callers may treat an uncertain vote
        as advisory in the same way a self-conflicting tool's vote is.
        """
        lab = (label or "").strip()
        if lab.startswith("#"):
            lab = lab[1:]
        uncertain = False
        while lab.endswith("?"):
            lab, uncertain = lab[:-1], True
        # A '?' can also sit on the superfamily only: "DNA/hAT?-Ac"
        if "?" in lab:
            lab, uncertain = lab.replace("?", ""), True
        return lab, uncertain

    def lookup(self, label: str, tool_id: str) -> dict:
        """Resolve a raw label to a canonical record. Never raises.

        Unmapped labels return an uninformative record and are tallied in
        ``self.unmapped`` for the curation report -- they are never silently
        binned into a class.
        """
        lab, uncertain = self._normalize(label)
        cands = [lab]
        if "/" in lab:  # EDTA-style Class/CODE -- also try the bare code
            cands.append(lab.split("/", 1)[1])

        # Collect every matching rule across all candidate spellings, then take
        # the most specific one. Specificity is (depth of canonical path, length
        # of the matched rule) -- NOT the order candidates were tried. This is
        # what makes "MITE/DTT" resolve to Tc1Mariner via the DTT rule rather
        # than stopping at bare TIR via the MITE rule; EDTA labels MITEs that
        # way routinely, and the superfamily is the informative part.
        best, best_key = None, None
        for ci, cand in enumerate(cands):
            low = cand.lower()
            hits = []
            for tid in (tool_id, "*"):
                if (low, tid) in self.exact:
                    hits.append((self.exact[(low, tid)], len(low), tid, True))
            for pref, tid, rec in self.prefix:
                if tid in (tool_id, "*") and low.startswith(pref):
                    hits.append((rec, len(pref), tid, False))
            for rec, rule_len, tid, is_exact in hits:
                key = (
                    # An exact rule matching the FULL label (ci == 0) outranks
                    # any prefix rule, regardless of path depth: it is by
                    # definition the most specific statement the map can make
                    # about that label. Without this, REPET's compound Wicker
                    # codes break: the exact rule `RLX|RYX -> repeat:TE:ClassI`
                    # (the deepest level the two codes share) would lose to the
                    # prefix rule `RLX* -> ...:LTR`, asserting an order the
                    # tool itself left ambiguous. Restricted to ci == 0 because
                    # an exact hit on the bare-code FALLBACK candidate (the
                    # `Class/CODE` split) must not outrank a rule written for
                    # the full label: `SINE/tRNA` splits to `tRNA`, whose exact
                    # rule (repeat:multigene:tRNA) would otherwise beat the
                    # intended SINE/tRNA* prefix rule.
                    is_exact and ci == 0,
                    len(split_path(rec["canonical_path"])),  # deeper path wins
                    rule_len,                                # longer rule wins
                    tid != "*",                              # tool-specific wins
                )
                if best_key is None or key > best_key:
                    best, best_key = rec, key
        if best is not None:
            if uncertain:
                # Keep the class, demote the confidence, and mark it so the
                # summariser can treat the vote as advisory.
                rec = dict(best)
                rec["confidence"] = "low" if rec["confidence"] == "high" else "low"
                rec["uncertain"] = True
                return rec
            return best

        key = (label, tool_id)
        self.unmapped[key] = self.unmapped.get(key, 0) + 1
        return {"canonical_path": "repeat", "canonical_v2": "", "is_te": "unknown",
                "confidence": "none", "raw_rule": "<unmapped>", "uncertain": uncertain}

    def unmapped_report(self) -> pd.DataFrame:
        rows = [{"raw_label": lab, "tool_id": tid, "n_hits": n}
                for (lab, tid), n in self.unmapped.items()]
        return (pd.DataFrame(rows, columns=["raw_label", "tool_id", "n_hits"])
                .sort_values("n_hits", ascending=False))


@dataclass
class Palette:
    rules: list[tuple[str, str, str]] = field(default_factory=list)  # path, rgb, label
    conflict_color: str = "0,0,0"

    @classmethod
    def load(cls, path="config/palette.tsv") -> "Palette":
        df, _ = _read_tsv(path)
        p = cls()
        for _, r in df.iterrows():
            if r["match_path"] == "CONFLICT":
                p.conflict_color = r["color"]
                continue
            p.rules.append((r["match_path"], r["color"], r["label"]))
        p.rules.sort(key=lambda x: -len(x[0]))  # longest prefix wins
        return p

    def color_for(self, canonical_path: str) -> str:
        for pref, rgb, _ in self.rules:
            if canonical_path.startswith(pref):
                return rgb
        return "128,128,128"

    def label_for(self, canonical_path: str) -> str:
        for pref, _, lab in self.rules:
            if canonical_path.startswith(pref):
                return lab
        return "Repeat"


# --------------------------------------------------------------------------
# Agreement calculator
# --------------------------------------------------------------------------

def split_path(p: str) -> list[str]:
    return [x for x in (p or "").split(":") if x]


@dataclass
class Agreement:
    consensus_path: str      # deepest path all informative voters support
    agree_depth: str         # level name of that depth, or "none"
    conflict_depth: str      # level name where voters first diverge, or "" if none
    conflict_detail: str     # human-readable "rm2:Gypsy|edta:Copia", or ""
    n_informative: int       # voters that asserted anything beyond "repeat"
    n_abstain: int           # voters that were Unknown/NA


def agreement(votes: dict[str, str], weights: dict[str, float] | None = None) -> Agreement:
    """Deepest level at which all informative voters agree.

    ``votes`` maps tool_id -> canonical path. Paths equal to "repeat" (or empty)
    are treated as ABSTENTIONS: they neither deepen nor break consensus. This is
    the policy that keeps a RepeatModeler "Unknown" from demoting a locus that
    three tools call LTR/Gypsy.

    ``weights`` optionally downweights a tool's vote (used for tools that
    contradict themselves at the locus -- see ``selfConflict`` in overlap.py).
    A weight of 0 makes the vote advisory only: it is excluded from consensus
    but still reported in ``conflict_detail`` so nothing is hidden.
    """
    weights = weights or {}
    informative = {t: split_path(p) for t, p in votes.items()
                   if p and p not in UNINFORMATIVE_PATHS}
    abstain = len(votes) - len(informative)

    voting = {t: p for t, p in informative.items() if weights.get(t, 1.0) > 0}
    advisory = {t: p for t, p in informative.items() if weights.get(t, 1.0) <= 0}

    if not voting:
        # Everyone abstained or was downweighted: repeat is asserted, class is not.
        detail = ";".join(f"{t}:{':'.join(p)}(downweighted)" for t, p in advisory.items())
        return Agreement("repeat", "none", "", detail, 0, abstain)

    depth = 0
    max_depth = min(len(p) for p in voting.values())
    while depth < max_depth and len({p[depth] for p in voting.values()}) == 1:
        depth += 1

    consensus = ":".join(next(iter(voting.values()))[:depth]) if depth else "repeat"

    conflict_depth, detail = "", ""
    if depth < max(len(p) for p in voting.values()):
        divergent = {p[depth] for p in voting.values() if len(p) > depth}
        if len(divergent) > 1:
            conflict_depth = LEVELS[depth] if depth < len(LEVELS) else f"depth{depth}"
            detail = "|".join(
                f"{t}:{p[depth]}" for t, p in sorted(voting.items()) if len(p) > depth
            )
    if advisory:
        adv = "|".join(f"{t}:{':'.join(p)}(downweighted)" for t, p in sorted(advisory.items()))
        detail = f"{detail}||{adv}" if detail else adv

    agree_depth = LEVELS[depth - 1] if 0 < depth <= len(LEVELS) else (
        "none" if depth == 0 else f"depth{depth}")
    return Agreement(consensus, agree_depth, conflict_depth, detail,
                     len(informative), abstain)
