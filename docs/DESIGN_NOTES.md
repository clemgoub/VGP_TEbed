# Design notes

Original brainstorming and format-design notes for the integrated repeat track,
preserved from the pre-implementation README. The sections that were open
questions at the time are now answered in
[SPECIFICATION.md](SPECIFICATION.md); where that is the case it is marked
inline. The **input format requirement** below is still current and normative —
that is the format every tool must provide.

> These are working notes, not a specification. They record how we got here,
> including reasoning that was later revised.

---

## Concept / Brainstorming

We are looking at annotating each base of a genome. For each base we ask: is it a TE? If yes, what type of TE?
Hierarchical structure of annotation. We have multiple tools, and each tools can have overlapping annotations.
1. The smallest unit is a hit made by one tool, represented by an interval.
2. The same tool can have overlapping hits
3. Multiple tools can have overlapping hits

For a given base of the genome, we can ask how confident we are that this is a repeat. We can do that by asking how many tool agree that this is a repeat. But what if one tools has 2 or more overlaping hits for that base? Is is further evidence, or less or is it irrelevant (I think it is irrelevant and instead it says that there is redundancy in the library).

For a given base, we have multiple questions:
- is it a repeat?
- what type of repeat? (classification)
  - TE / non-TE
  - class I / class II
  - order
  - finer grain?
- what is the consensus in the library?
- what is the divergence to consensus?
- what is the orientation?
- what is the position in the consensus

---

## Deliverable Specification

> **Status: superseded.** The delivered format is BED12+11 in a bigBed, not a
> plain `.bed` — UCSC hubs require binary formats, and BED12 was needed for the
> thick/thin geometry that carries boundary disagreement. The field list below
> was the starting sketch; the delivered schema is in
> [SPECIFICATION.md §6](SPECIFICATION.md#6-track-reference). The BED12 field
> descriptions here remain an accurate reference for the standard columns.

describe here the format specification of the integrated track (it will be `.bed`)

### `.bed` track

UCSC eventually convert any track to bed for display, so Hiram recommend to focus on delivering a bed track.

The first 12 fields are constrained by the bed standard:

1. `chr`: chromosome
2. `start`: feature start  (0-based)
3. `end`: feature end (1-based)
4. `name`: unique instance name
5. `score`: *undefined*
6. `strand`: `+`, `-` or `.` if conflict
7. `thickStart`: The starting position at which the feature is drawn thickly (for example, the start codon in gene displays). When there is no thick part, thickStart and thickEnd are usually set to the chromStart position.
8. `thickEnd`: The ending position at which the feature is drawn thickly (for example the stop codon in gene displays).
9. `itemRgb`: An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
10. `blockCount`: The number of blocks (exons) in the BED line.
11. `blockSizes`: A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
12. `blockStarts`: A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount.

The rest is up for grabs. The current order is irrelevant, the goal is to list the information we need to display.

13. `support`: number of methods supporting a repeat in this region `<int>`
14. `methods`: list of methods supporting a repeat in this regions `<string>` (eg: `RM2;EDTA`)
15.

---

## Input format requirement

> **Status: promoted.** This section is now maintained as the normative
> specification in [INPUT_FORMAT.md](INPUT_FORMAT.md), which documents the
> validation rules the pipeline actually applies. The text below is the
> original design note, kept for context.

The goal is to be able to link each hit to a consensus in the library produced by each tool. This is straightforward using Repeatmasker output (`.out`), but EDTA outputs are more complex as it blends structural and homology calls. Tools using standalone `RepeatMasker` for annotation (e.g. `RepeatModeler2`, `Pantera`) will use a conversion script `.out --> .bed` while `EDTA` outputs will be converted to `.bed` with a custom convertor in order to account for structure-based annotation.

As standard input for the integrated track, each method must provide a custom `.bed` with the following specification:

| # | Column | Description |
|---|--------|-------------|
| 1 | `chrom` | Query sequence name (scaffold/chromosome). |
| 2 | `chromStart` | Match start, 0-based (RepeatMasker query begin − 1). |
| 3 | `chromEnd` | Match end, half-open (RepeatMasker query end). |
| 4 | `name` | Repeat element name (RepeatMasker "matching repeat"). |
| 5 | `score` | Smith–Waterman score, capped to the BED-legal range 0–1000. Use `SW_score` (col 7) for the true value. |
| 6 | `strand` | `+` for `+`, `-` for RepeatMasker's `C` (complement). |
| 7 | `SW_score`* | Raw Smith–Waterman alignment score (uncapped). |
| 8 | `perc_div`* | Percent substitutions vs. the consensus. |
| 9 | `perc_del`* | Percent bases deleted vs. the consensus. |
| 10 | `perc_ins`* | Percent bases inserted vs. the consensus. |
| 11 | `query_left`* | Bases left in the query after the match (RepeatMasker parenthesized value, unwrapped). |
| 12 | `repeat_class_family` | Repeat class/family (e.g. `LINE/L1`, `LTR/ERV1`). |
| 13 | `repeat_start`* | Match start position in the repeat consensus. |
| 14 | `repeat_end`* | Match end position in the repeat consensus. |
| 15 | `repeat_left`* | Bases left in the repeat consensus after the match (parenthesized value, unwrapped). |
| 16 | `ID` | copy ID, linking fragments of a single interrupted insertion. |

\*`RepeatMasker`-specific fields (structure-based annotations do not have such information), `NA` otherwise.

### Input format conversion scripts

Conversion scripts can be found in `scripts/`

#### `rmout2bed.py`

Convert `RepeatMasker` `.out` file into `.bed`. Tab-delimited, BED6+10. The first line is a #-prefixed header (skipped by bedtools). Coordinates are converted from RepeatMasker's 1-based, fully-closed intervals to BED's 0-based, half-open convention (chromStart = query_begin − 1, chromEnd = query_end). Records that would produce a zero-length or negative interval (chromEnd ≤ chromStart) are dropped, and malformed lines (fewer than 15 fields) are skipped; both are reported as counts to stderr.

| # | Column | Description |
|---|--------|-------------|
| 1 | `chrom` | Query sequence name (scaffold/chromosome). |
| 2 | `chromStart` | Match start, 0-based (RepeatMasker query begin − 1). |
| 3 | `chromEnd` | Match end, half-open (RepeatMasker query end). |
| 4 | `name` | Repeat element name (RepeatMasker "matching repeat"). |
| 5 | `score`* | Smith–Waterman score, capped to the BED-legal range 0–1000. Use `SW_score` (col 7) for the true value. |
| 6 | `strand` | `+` for `+`, `-` for RepeatMasker's `C` (complement). |
| 7 | `SW_score` | Raw Smith–Waterman alignment score (uncapped). |
| 8 | `perc_div` | Percent substitutions vs. the consensus. |
| 9 | `perc_del` | Percent bases deleted vs. the consensus. |
| 10 | `perc_ins` | Percent bases inserted vs. the consensus. |
| 11 | `query_left` | Bases left in the query after the match (RepeatMasker parenthesized value, unwrapped). |
| 12 | `repeat_class_family` | Repeat class/family (e.g. `LINE/L1`, `LTR/ERV1`). |
| 13 | `repeat_start` | Match start position in the repeat consensus. |
| 14 | `repeat_end` | Match end position in the repeat consensus. |
| 15 | `repeat_left` | Bases left in the repeat consensus after the match (parenthesized value, unwrapped). |
| 16 | `RM_ID` | RepeatMasker element ID, linking fragments of a single interrupted insertion. |

> Repeat-consensus coordinate order in the source `.out` depends on strand; the script reorders columns 13–15 to repeat_start, repeat_end, repeat_left according to strand.

> Parenthesized RepeatMasker values ((n), meaning "remaining") are stripped to plain integers in columns 11 and 15.

> Only columns 1–6 are standard BED6; colums 7–16 are extra RepeatMasker fields carried through and are ignored by standard BED tools.

#### `XXX` EDTA to `.bed` (to add)

Convert `EDTA` outputs (gff3?) into `.bed`  -- please add description of EDTA specific fields here.

| # | Column | Description |
|---|--------|-------------|
| 1 | `chrom` | Query sequence name (scaffold/chromosome). |
| 2 | `chromStart` | Match start, 0-based (RepeatMasker query begin − 1). |
| 3 | `chromEnd` | Match end, half-open (RepeatMasker query end). |
| 4 | `name` | Repeat element name (RepeatMasker "matching repeat"). |
| 5 | `score`* | Smith–Waterman score, capped to the BED-legal range 0–1000. Use `SW_score` (col 7) for the true value. |
| 6 | `strand` | `+` for `+`, `-` for RepeatMasker's `C` (complement). |
| 7 | `SW_score`* | Raw Smith–Waterman alignment score (uncapped). |
| 8 | `perc_div`* | Percent substitutions vs. the consensus. |
| 9 | `perc_del`* | Percent bases deleted vs. the consensus. |
| 10 | `perc_ins`* | Percent bases inserted vs. the consensus. |
| 11 | `query_left`* | Bases left in the query after the match (RepeatMasker parenthesized value, unwrapped). |
| 12 | `repeat_class_family` | Repeat class/family (e.g. `LINE/L1`, `LTR/ERV1`). |
| 13 | `repeat_start`* | Match start position in the repeat consensus. |
| 14 | `repeat_end`* | Match end position in the repeat consensus. |
| 15 | `repeat_left`* | Bases left in the repeat consensus after the match (parenthesized value, unwrapped). |
| 16 | `TE_ID` | `<PLEASE UPDATE WITH EDTA DEFINITION OF TE_ID>` |

\*`RepeatMasker`-specific fields (structure-based annotations do not have such information), `NA` otherwise.

---

## Integration tool

> **Status: implemented.** The pipeline is the `vgptrack` package; see the
> [README quickstart](../README.md#quickstart) and
> [SPECIFICATION.md](SPECIFICATION.md).

---

## Scope-gated classification votes

A tool only votes on *classification* where the locus falls inside its declared
`scope`. This is separate from the support denominator (`add_eligibility`),
which had honoured scope from the start; agreement resolution did not, and the
gap only became visible when the first narrow-scope tool was added.

The failure mode: FasTAN calls `repeat:tandem`, which is a real assertion, not
an abstention — so the existing abstention rule (class id 0 never breaks
agreement) could not absorb it. At a locus three TE tools called
`ClassI:LTR`, FasTAN's `tandem` registered as a genuine dispute at depth 1 and
collapsed the consensus to bare `repeat`. Measured against a three-tool
baseline on identical segment boundaries, this destroyed **20.7 Mb** of
agreement and manufactured **14.4 Mb** of conflict, against 2.7 Mb gained.

Resolution now runs in two passes:

1. Resolve over unrestricted (`general_homology`) tools to establish what kind
   of locus this is.
2. Admit each scope-restricted tool's class vote only where that provisional
   answer lies inside its scope.

Where the unrestricted tools say nothing, the restricted voter making the
**deepest** assertion anchors. Admitting all of them instead lets tools in
disjoint scopes cross-veto — EDTA calling CACTA and FasTAN calling `tandem` at
the same locus scored as a dispute, collapsing a further 5.6 Mb, though a
tandem array overlapping a TE is ordinary biology and neither tool can
adjudicate the other's question. Depth breaks the tie because naming a TE
superfamily is a stronger claim than reporting that an array exists.

After both changes, adding FasTAN is strictly additive: 6.85 Mb of agreement
gained, none lost, no new conflicts. Genome-wide, sequence left unresolved at
bare `repeat` falls 96.9 → 47.2 Mb, TE-classified sequence rises 313.5 → 357.3
Mb, and tandem-classified sequence rises 52.9 → 57.6 Mb. The gate applies to EDTA too (`structural_te`), which was not the
original intent but is the same principle: on the bundled test slice it
deepens 9,087 bp of satellite sequence that a TE-only tool had been
overriding, and shallows nothing.

**Consequence for new tools.** A tool declaring a narrow `scope` in
`config/tools.tsv` contributes existence evidence everywhere but classification
evidence only within scope. Declaring `general_homology` for a specialist tool
re-opens this bug. `_in_scope` in `segment.py` and `ToolSet.eligible` in
`vocab.py` answer the same question for the class vote and the denominator
respectively — keep them in step.

---

## Open questions

Carried forward from the brainstorming above; these are the most useful places
for outside input.

- **Is multi-tool support the right confidence measure?** Tools that share a
  library or an algorithmic lineage are not independent evidence — RM2 and EDTA
  share RepeatMasker machinery and agree at J = 0.70, against J ≈ 0.53 for
  either versus Pantera. Support currently counts tools equally.
- **Should tools be weighted?** Currently they are not, deliberately: weighting
  requires a benchmark we do not have.
- **Divergence is not comparable across tools** — each reports against its own
  consensus. The current signal track mixes them. Should it be per-tool only?
- **Position in consensus and orientation** are carried per tool but not yet
  summarized across tools.
- **Relationship to the consensus library** (`1.5` in the original notes) is
  carried as `consensusRange` per tool; linking to the actual library FASTA is
  not implemented.
- **How should the vocabulary handle tandem repeats and satellites?** They are
  mapped as non-TE repeats, but the TE-centric hierarchy fits them poorly.
