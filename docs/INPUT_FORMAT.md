# BED16 input format

**This is the normative input specification.** Every tool feeding the
integrated track must provide its annotation in this format. If you are writing
a converter for a new tool, this document is the contract.

> Status: stable in structure, provisional in vocabulary. The columns below are
> not expected to change. The permitted values of `repeat_class_family` will be
> revised when the unified Dfam/Repbase classification is released.

---

## 1. File shape

Tab-separated, one row per annotated interval, 16 columns, in the order below.

- **Header line: optional.** A first line beginning with `#`, or whose first
  field is `chrom`, is treated as a header and skipped. Otherwise the first line
  is read as data. Column *names* in a header are ignored — **position is
  authoritative**.
- **Compression: optional.** `.gz` is read transparently.
- **Sorting: not required.** The pipeline sorts internally.
- **Missing values:** `NA`, `na`, or `.` in any column marked *nullable*.
- Only columns 1–6 are standard BED6, so the file is readable by `bedtools` and
  loadable in a genome browser as-is. Columns 7–16 are carried through and
  ignored by standard BED tools.

## 2. Columns

| # | Column | Type | Null? | Description |
|---|---|---|---|---|
| 1 | `chrom` | string | no | Query sequence name. Must match `chrom.sizes` **or** resolve through the assembly's chromAlias — see §4. |
| 2 | `chromStart` | int | no | Match start, **0-based**. |
| 3 | `chromEnd` | int | no | Match end, **half-open**. Must be > `chromStart`. |
| 4 | `name` | string | no | Repeat element / family name as the tool calls it. |
| 5 | `score` | int 0–1000 | no | BED-legal score. Capped; use col 7 for the true value. |
| 6 | `strand` | `+` / `-` | no | RepeatMasker `C` must be converted to `-`. |
| 7 | `SW_score` | int | yes | Raw Smith–Waterman score, uncapped. |
| 8 | `perc_div` | float | yes | Percent substitution vs. the consensus. **Feeds the divergence track.** |
| 9 | `perc_del` | float | yes | Percent deleted vs. consensus. |
| 10 | `perc_ins` | float | yes | Percent inserted vs. consensus. |
| 11 | `query_left` | int | yes | Bases remaining in the query after the match. Parentheses stripped. |
| 12 | `repeat_class_family` | string | no | Classification. **Feeds harmonization** — see §3. |
| 13 | `repeat_start` | int | yes | Match start in the repeat consensus. |
| 14 | `repeat_end` | int | yes | Match end in the repeat consensus. |
| 15 | `repeat_left` | int | yes | Bases remaining in the consensus after the match. Parentheses stripped. |
| 16 | `hit_id` | string | no | Copy ID linking fragments of one interrupted insertion. |

Columns 7–11 and 13–15 are RepeatMasker-specific. **Structural callers should
write `NA`, not zero.** `NA` means "this tool cannot report this"; `0` is a
measurement. A structural LTR call with `perc_div = 0` claims a perfect match to
a consensus that was never compared.

### Two columns do real work

**`repeat_class_family` (col 12)** drives classification agreement. It is mapped
onto a shared hierarchy through `config/class_map.tsv`. Emit the tool's own
label — do **not** pre-translate it into someone else's vocabulary. Unmapped
labels are reported in `qc/unmapped_labels.tsv` rather than silently dropped;
add a mapping row for anything that appears there.

Use `Unknown` (or `NA`) when the tool genuinely does not classify. This is
treated as **abstention** — it neither deepens agreement nor creates conflict.

**Abstain wherever your tool abstains. Never fill in a plausible class.** A
guessed label cuts both ways, and the second case is the dangerous one:

- Guess wrong and you manufacture **false disagreement** — the locus reports a
  conflict that no tool actually raised.
- Guess *right* and you manufacture **false agreement** — the locus reports
  consensus at a depth only one tool truly asserted, and `agreeDepth` reads
  deeper than the evidence supports.

False agreement is worse because it looks like a result. By the time a label
reaches segmentation it is indistinguishable from a tool's own call, so nothing
downstream can flag it. If a coarse label is all your tool supports, emit the
coarse label: `repeat:tandem` with three tools abstaining below it is an honest
statement about what is known.

Where a size- or length-based convention exists that users may still want,
expose it as an **opt-in flag**, not a default, and say in the flag's help what
it does to agreement depth. `scripts/fastan2bed.py --classify-period` is the
worked example.

**`hit_id` (col 16)** links fragments of one interrupted insertion. It must be
unique per tool per assembly. EDTA additionally encodes provenance here
(`TE_homo_*` = homology, `TE_struc_*` / `LTRRT_*` = structural), which the
pipeline preserves as the `evidence` field.

## 3. Validation, and what gets dropped

Applied per tool at ingest; all counts land in `qc/ingest_stats.tsv`.

| Check | Action | Reported as |
|---|---|---|
| `chromEnd` ≤ `chromStart`, or non-numeric coordinates | row dropped | `n_bad_coord` |
| `chrom` not in `chrom.sizes` after alias resolution | row dropped | `n_unknown_seq` (+ first 10 names) |
| `chromEnd` beyond sequence length | **clamped** to length | `n_clamped` |
| `repeat_class_family` not in the class map | row kept, label reported | `qc/unmapped_labels.tsv` |

**Always check `n_unknown_seq` after adding a tool.** A converter that emits the
wrong naming authority produces a valid-looking file where every row is silently
discarded. A large `n_unknown_seq` with a plausible-looking name list is the
signature.

Overlapping calls from a single tool are **not** an error and are not
deduplicated. Support counts distinct tools via a bitmask, so one tool's
redundant calls cannot inflate it; self-overlap is separately classified as
nested, library-redundant or self-contradictory.

## 4. Sequence naming

Tools disagree about naming authority — in the black goby run, EDTA led with
`CATOHO*` while RepeatModeler2 and Pantera used INSDC `OX*` accessions for the
same sequences. Emit whatever your tool produces and pass `--alias` at build
time; every alias in that file resolves to the assembly's primary name.

Without `--alias`, non-matching names are dropped at ingest, so it is
effectively required on real assemblies.

## 5. Minimal valid example

```
#chrom	chromStart	chromEnd	name	score	strand	SW_score	perc_div	perc_del	perc_ins	query_left	repeat_class_family	repeat_start	repeat_end	repeat_left	hit_id
OX637595.1	15848	16090	(ACTACT)n	171	+	171	9.4	0.0	0.0	76485642	Simple_repeat	1	242	0	1
OX637595.1	21044	21395	L2-3_DR	421	-	421	22.1	1.4	0.7	76480337	LINE/L2	2891	3247	1204	2
```

A structural caller with no consensus alignment:

```
OX637595.1	48210	53887	LTR_retro_1	NA	+	NA	NA	NA	NA	NA	LTR/Gypsy	NA	NA	NA	TE_struc_1
```

Working examples: `tests/data/*.bed.gz`.

## 6. Converters

In `scripts/`. See [DESIGN_NOTES.md](DESIGN_NOTES.md#input-format-conversion-scripts)
for per-converter detail.

| Tool | Source format supplied | Script | Production expectation |
|---|---|---|---|
| RepeatModeler2 | RepeatMasker `.out` | `rmout2bed.py` | `.out` — no change |
| Pantera | BED16 (already conformant) | none needed | `.out` from the RepeatMasker run |
| fastLTR | RepeatMasker `.out` | `rmout2bed.py` | `.out` — no change |
| EDTA | GFF3 (`*.TEanno.gff3`) | `edtagff2bed.py` | converter until EDTA emits `.out`/BED16 |
| FasTAN | native BED | `fastan2bed.py` | converter until FasTAN emits BED16 |
| LTRDeNovo | native GFF3 (NGSEP) | `ltrdenovogff2bed.py` | converter until it emits `.out`/BED16 |
| Satellome | native BED5 | `satellome2bed.py` | converter until it emits BED16 |

**Three of seven tools already supply a standard input.** RepeatModeler2 and
fastLTR provide RepeatMasker `.out`; Pantera's was supplied as a conformant
BED16 directly. Those need no per-tool code — `rmout2bed.py` is the generic
`.out` reader, not a fastLTR-specific script.

The remaining four are the real format gap: EDTA, FasTAN, LTRDeNovo and
Satellome emit native formats, and each converter encodes a decision the tool
should be making itself (which GFF3 feature level represents the element;
whether an `identity` field is consensus divergence or a within-element
measure; whether a length filter applied upstream is part of the annotation).
Each `docs/<TOOL>.md` records that decision so it can be raised upstream.

### Prototype-only accommodations

Distinct from format conversion, and expected to disappear:

- **`infer_rm_rename.py`** — the pilot fastLTR `.out` had query names shortened
  to `_J0000000` placeholders by an upstream `.fa.mod` step, with no mapping
  published; names were recovered by matching sequence lengths (120/126 unique;
  7 hits on size-colliding unplaced scaffolds dropped as unresolvable). Future
  fastLTR output carries the assembly's real contig IDs, so this step retires
  and fastLTR ingests through plain `rmout2bed.py`.

If you are a tool author reading this: emitting RepeatMasker `.out`, or the
BED16 of §2, removes your tool from the converter table entirely.

Converting from RepeatMasker `.out`, remember: coordinates are 1-based
fully-closed, so `chromStart = query_begin - 1`; strand `C` becomes `-`;
parenthesized "remaining" values are stripped to plain integers; and columns
13–15 are ordered by strand in the source, so they must be reordered to
`repeat_start`, `repeat_end`, `repeat_left`.

### Checking a new converter

```bash
python -m vgptrack.cli build --assembly TEST \
    --sizes your.chrom.sizes --alias your.chromAlias.txt \
    --bed yourtool=yourtool.bed --out /tmp/hubtest
```

Then read `qc/ingest_stats.tsv`: `n_hits` should be close to `n_raw`,
`n_unknown_seq` should be 0, and any label in `qc/unmapped_labels.tsv` needs a
row in `config/class_map.tsv`.
