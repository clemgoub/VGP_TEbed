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

## 6. Converters — the supported ingest path for known tools

**Policy: for a supported tool, the native output + its converter is the
preferred input.** The converter is written against the tool's own format and
documents exactly what it keeps, transforms and drops — so nothing is lost to
a hand conversion, and the load-bearing decisions (feature level, identity
semantics, coordinate convention) are made once, reviewed, and versioned in
this repo rather than re-made ad hoc by every data submitter. The BED16 of §2
is the **interchange format for tools without a converter**: any tool can join
the build today by supplying it, and a converter can be added later when the
tool's native format warrants one.

In `scripts/`. See [DESIGN_NOTES.md](DESIGN_NOTES.md#input-format-conversion-scripts)
for per-converter detail.

| Tool | Native input | Script | Native fields beyond BED16, and where they go |
|---|---|---|---|
| RepeatModeler2 | RepeatMasker `.out` | `rmout2bed.py` | none carried beyond BED16; the overlap `*` flag is dropped (see residual losses) |
| fastLTR | RepeatMasker `.out` | `rmout2bed.py` | same — overlap `*` flag dropped |
| Pantera | RepeatMasker `.out` (BED16 supplied for the pilot) | `rmout2bed.py` | same — overlap `*` flag dropped |
| EDTA | GFF3 (`*.TEanno.gff3`) | `edtagff2bed.py` | `method` (structural/homology) → `hit_id` prefix; `identity` → `perc_div` (structural LTR identity deliberately NOT emitted as divergence — within-element, not vs consensus) |
| FasTAN | native BED (period, identity) | `fastan2bed.py` | `period` → `name` (`tandem_p<period>`); unit identity → `score` + `perc_div` (flagged `divergence_only` in the manifest) |
| LTRDeNovo | native GFF3 (NGSEP) | `ltrdenovogff2bed.py` | `method` → `name` suffix; **LTR/TSD sub-features currently dropped** — see below |
| Satellome | native BED5 | `satellome2bed.py` | none — col 5 is redundant (verified) |
| TRF | GenArk `simpleRepeat.bb` (bigBed 4+12) | `simplerepeat2bed.py` | period+copyNum → `name` (`p<N>_x<N>_<unit>`); TRF score → `SW_score`; 100−perMatch → `perc_div` (`divergence_only`); base comp/entropy/unit seq not carried (rederivable from the public track) |
| WindowMasker | GenArk `windowMasker.bb` (bigBed 3) | `windowmasker2bed.py` | nothing to carry — bare intervals; emits `Unknown` (pure existence, abstains from class) |
| REPET | TEannot GFF3 (`match`/`match_part`) | `repetgff2bed.py` | one row per `match_part` (fragment), copy id → `hit_id`; fragment `Identity` → `perc_div` (consensus divergence; match-level `AlignIdentity` is a placeholder and NOT used); `Wcode` → `repeat_class_family` verbatim incl. compound codes; **`OtherTargets` and `TargetDescription` evidence detail (coding/struct annotations) not carried** — see docs/REPET.md |

Where a native field has no BED16 column, the converters route it into `name`
or `hit_id` so it survives into the per-tool track mouseover. Known residual
losses, kept honest here rather than hidden:

- **LTRDeNovo sub-features**: `five_prime_LTR` / `three_prime_LTR` spans
  (228 structural elements) and `target_site_duplication` spans + sequences
  (118) are not representable in a flat BED16 row. The element footprint and
  detection method survive; the internal architecture does not. If wanted,
  the natural home is BED12 blocks on the per-tool track (thick = internal
  domain, blocks = LTRs), which is a display change, not an ingest change.
- **RepeatMasker `.out` overlap flag** (final `*` column) is not carried.

### Prototype-only accommodations

Distinct from format conversion, and expected to disappear:

- **`infer_rm_rename.py`** — the pilot fastLTR `.out` had query names shortened
  to `_J0000000` placeholders by an upstream `.fa.mod` step, with no mapping
  published; names were recovered by matching sequence lengths (120/126 unique;
  7 hits on size-colliding unplaced scaffolds dropped as unresolvable). Future
  fastLTR output carries the assembly's real contig IDs, so this step retires
  and fastLTR ingests through plain `rmout2bed.py`.

If you are a tool author whose tool is not in the table: the fastest route in
is the BED16 of §2 (or RepeatMasker `.out`). If your native format carries
information BED16 cannot express, open an issue — that is the case for writing
a converter.

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
