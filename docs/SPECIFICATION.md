# VGP Multi-Tool Repeat Track Hub — Format Specification

**Status:** prototype, validated end-to-end on `GCA_951799975.1` (*Gobius niger*, black goby).
**Vocabulary version:** `dfam-rm-2026.1-draft` (provisional).
**Validation:** `hubCheck` clean, full track validation.

---

## 1. What this hub is

Repeat annotation from several independent tools, published as **evidence rather
than as a verdict**. Nothing is filtered on the basis of agreement; a locus called
by one tool alone is displayed as prominently as one called by all. Where tools
disagree — about presence, about extent, or about classification — the
disagreement is encoded in the display rather than resolved away.

This is deliberately *not* a merged annotation. A merge requires deciding whose
boundary and whose classification wins, and that decision is exactly what the
official pipeline is being built to make. Until then, the honest product is the
evidence.

---

## 2. Directory layout

Per the UCSC contributed-tracks layout:

```
hub.txt
genomes.txt
groups.txt
documentation.html
GCA_951799975.1/
    trackDb.txt
    documentation.html
    repeatSummary.bb          bigBed 12+11
    repeatSupport.bw          bigWig
    repeatSupportFrac.bw      bigWig
    repeatDivergence.bw       bigWig
    toolUnique.bb             bigBed 9+4
    repeat_rm2.bb             bigBed 9+7
    repeat_edta.bb            bigBed 9+7
    repeat_pantera.bb         bigBed 9+7
    repeat_fastltr.bb         bigBed 9+7  (empty — did not run)
```

**Two rules govern this layout and both have consequences in the code:**

1. **Filenames are identical in every assembly directory.** One `trackDb.txt`
   serves all assemblies, so it cannot reference an assembly and cannot reference
   a file that might be missing. A tool that did not run therefore gets an
   **empty but schema-valid bigBed**, not a missing file. `repeat_fastltr.bb`
   above is such a file: 0 features, full 18-field schema, opens correctly.
2. **Binary formats only.** No BED, no GFF, no gzip. Everything is bigBed or
   bigWig, built with the UCSC utilities and indexed for name search.

`genomes.txt` is **appended** to, never overwritten — building a second assembly
must not erase the first.

---

## 3. The unit of display: the element

The pipeline computes per-base state, then run-length encodes it into
**segments** — maximal runs of bases with identical tool support, identical
consensus class, and identical conflict state. On the black goby this yields
4.31 M segments over 431.3 Mb.

Segments are too fine to browse: the median is 42 bp, because adjacent segments
usually differ only in *which* tools cover them. That fragmentation is
**boundary disagreement**, which is signal, not noise.

So the display feature is the **element**: a maximal run of adjacent covered
bases whose consensus classification stays compatible and whose conflict state
does not change. 1.92 M elements, median 116 bp. Boundary disagreement is
absorbed into the element and re-expressed as BED12 geometry:

| Geometry | Meaning |
|---|---|
| **Thick block** | longest contiguous run where *every eligible tool* called it |
| **Thin extent** | union of all tool calls at this locus |
| **Entirely thin** | no position in the element has full eligible support |

Two guards keep elements honest:

- **Sliver absorption (≤20 bp).** A segment shorter than this neither starts a
  new element nor blocks its neighbours from joining across it. Without this,
  single-base boundary jitter splits one repeat into three features.
- **Conflict state must match to merge.** Class compatibility alone over-fuses:
  a truncated consensus path is a prefix of a longer one, so a conflicted stretch
  would absorb a clean neighbouring element and poison its classification.

**No bases are lost.** The element table's total covered bases equals the
segmentation's exactly (431,250,929 bp), and this is asserted at build time —
the CLI aborts if the summary bigBed's `basesCovered` disagrees with the
segmentation.

Full per-base resolution remains available in the three bigWig signals, which
are **not** display-merged.

---

## 4. Support: what "how many tools" means

**Support counts distinct tools, never calls.** Support is carried as a bitmask
with one bit per tool, so two overlapping calls from one tool set the same bit
and contribute 1. This is structural, not a rule applied afterwards — it cannot
be bypassed by a tool that emits redundant or nested annotations.

**The denominator is the eligible tool count, not the tool count.** A structural
LTR finder that cannot detect SINEs is not counted as dissenting at a SINE. For
each locus, `nEligible` is the number of tools whose declared detection scope
covers the consensus class.

Two corrections were needed here and both are in the code:

- *Observed evidence overrides declared scope.* If a tool declared out-of-scope
  for a class nonetheless emits calls there, it joins the denominator — otherwise
  support could exceed 100%. An assertion enforces `nSupport ≤ nEligible`.
- *Tools that did not run are excluded entirely* from both numerator and
  denominator, and from the "restricted denominator" flag.

---

## 5. Classification agreement: a depth, not a boolean

Tools classify at different resolutions. Measured on the real data:

| Tool | resolves to superfamily | stops at class | abstains |
|---|---|---|---|
| EDTA | 78% | — | — |
| RepeatModeler2 | 40% | — | 39% (`Unknown`) |
| Pantera | 8.6% | 62% | — |

A Pantera call of `DNA` alongside an RM2 call of `DNA/hAT-Ac` is **agreement at
class level**, not disagreement. Scoring it as disagreement would make Pantera
look like it contradicts everything.

Agreement is therefore reported as a **depth** on the hierarchy:

```
1 repeat  2 TE/non-TE  3 class  4 order  5 superfamily
```

- `agreeDepth` — deepest level at which all informative votes coincide.
- `conflictDepth` — shallowest level at which any two votes diverge; `-1` if none.
- **Abstention** (`Unknown`, `unclassified`, `NA`) neither deepens agreement nor
  creates conflict. A tool that declines to classify is not a dissenting vote.
- **Agreement is capped at the first conflict.** A feature cannot claim
  superfamily consensus and class-level conflict simultaneously; the consensus
  path is truncated to match.

`Class disputed` and `Repeat (unclassified)` are **different labels**: the first
means every tool classified and they contradict each other; the second means
nobody ventured a class.

### Self-conflict

A single tool overlapping itself with incompatible classifications is graded by
the depth of the disagreement, not treated as a binary. Overlap is first
classified as *nested*, *library-redundant*, or *self-contradictory*; only
self-contradiction at or below the configured depth (default 3, order level)
downweights that tool's **classification vote**. It never affects that tool's
**repeat call** — the tool still saw a repeat there.

---

## 6. Track reference

### repeatSummary.bb — `bigBed 12+11`

Standard BED12 plus:

| Field | Type | Meaning |
|---|---|---|
| `consensusClass` | string | harmonized class at the agreed depth |
| `nSupport` | int | distinct tools calling a repeat here |
| `nEligible` | int | tools capable of calling this class |
| `supportingTools` | string | which ones |
| `agreement` | string | depth name (`superfamily`, `class`, …) |
| `conflict` | string | depth name or `none` |
| `perToolClass` | string | each tool's deepest statement in this element |
| `coreBp` | int | length of the full-support core |
| `meanDivergence` | float | length-weighted, over tools reporting one |
| `flags` | string | plain-language caveats |
| `mouseOver` | string | composed hover line |

Indexed on `name` (`searchIndex`, `extraIndex`), `mouseOverField mouseOver`,
`itemRgb on`.

The hover answers, in order: **what it is → how many tools out of how many could
→ how deeply they agree → what conflicts → core size → divergence.**

### Signals — `bigWig`, full per-base resolution

`repeatSupport` (1–3, mean 2.23) · `repeatSupportFrac` (0.33–1.0) ·
`repeatDivergence` (0–59.2%, mean 12.2%).

### toolUnique.bb — `bigBed 9+4`

Intervals called by **exactly one** tool: 1.52 M features, 91.1 Mb. Built from
**segments**, not elements — an element carries the maximum support seen anywhere
within it, so building from elements lost 43% of genuinely single-tool bases.

The mouseover names the tools that were **capable** of calling that class and did
not, which is a stronger statement than a tool that could never have seen it.

### repeat_<tool>.bb — `bigBed 9+7`

Full unmodified per-tool output. **`rawClass` preserves the tool's own label
verbatim**, beside the harmonized `canonicalClass`. Also carries `hitId`,
`divergence`, `consensusRange`, `evidence` (homology vs structural), and
`overlapStatus` in plain language.

---

## 7. Class vocabulary

`config/class_map.tsv` is **data, not code**. It carries a version header, a
provisional-status marker, per-rule confidence, and a reserved `canonical_v2`
column for the forthcoming unified Dfam/Repbase classification. Re-running the
pipeline against a new map version rebuilds every track; nothing in the code
hardcodes a class name.

Coverage on the real data: **100% of 135 distinct observed labels**.

Design points worth knowing:

- **Uncertainty is first class.** A trailing `?` (RepeatMasker convention) means
  *best guess*, distinct from both a confident call and `Unknown`. The class is
  kept, confidence is demoted, and the vote is treated as advisory.
- **Lookup prefers rule specificity over candidate order.** `MITE/DTT` resolves
  through `DTT` → `Tc1Mariner`, not through `MITE*` → bare `TIR`. Getting this
  wrong discards superfamily information across a large fraction of EDTA's
  Class II calls.
- Colours are capped at 8 (UCSC guidance): colour carries **class**, geometry
  and hover carry everything finer.

---

## 8. Reproducing

```bash
python -m vgptrack.cli build \
    --assembly GCA_951799975.1 \
    --sizes   GCA_951799975.1.chrom.sizes \
    --alias   GCA_951799975.1.chromAlias.txt \
    --bed rm2=rm2.bed --bed edta=edta.bed --bed pantera=pantera.bed \
    --out hub --description "Gobius niger (black goby)" \
    --twobit https://hgdownload.soe.ucsc.edu/hubs/GCA/951/799/975/GCA_951799975.1/GCA_951799975.1.2bit
```

Swapping in real output from another tool is a `--bed` argument plus a row in
`config/tools.tsv`. Tools in the manifest without a `--bed` become empty
placeholders automatically. The build asserts the bp invariant and runs
`hubCheck` before exiting.

**Chromosome naming.** Tools lead with different naming authorities (INSDC vs
UCSC vs assembly-internal), so `--alias` is effectively required; without it,
sequence names that do not match `chrom.sizes` are dropped at ingest and reported
in the QC table.

---

## 9. Known limitations

- **Support measures agreement, not correctness.** Tools sharing a library or an
  algorithmic lineage are not independent evidence. RM2 and EDTA share
  RepeatMasker machinery, which inflates their pairwise agreement (J = 0.70)
  relative to either against Pantera (J ≈ 0.53).
- **Divergence is not comparable between tools** — each reports its own estimate
  against its own consensus.
- **Discovery has not saturated.** One tool finds 320 Mb on average, two find
  401 Mb, three find 431 Mb. The curve is still climbing; a fourth tool would
  add more.
- **Boundaries rarely agree**: only 32–39% of calls have a start within 10 bp of
  the nearest start in another tool (median offset 50–64 bp). This is why the
  thick/thin geometry exists.
- **The vocabulary is provisional** and will be revised on the unified
  Dfam/Repbase release.
- **fastLTR has not been run** on this assembly; its track is an empty
  placeholder and it is absent from all denominators.
