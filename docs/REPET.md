# REPET (TEannot) — tool 10, bit 9

Input: `Gnig_refTEs_redondant_features_merged.gff` (840 MB GFF3, source column
`Gnig_TEannotGr2_REPET_TEs`). REPET's TEannot pipeline maps the TEdenovo
de novo consensus library onto the assembly by homology. Converter:
`scripts/repetgff2bed.py`.

## Native format

Each TE **copy** is a two-level GFF3 feature (structure confirmed with the
REPET developers):

| type | count | role |
|---|---|---|
| `match` | 1,016,780 | one per copy: start–stop footprint of the (possibly interrupted) insertion |
| `match_part` | 1,128,771 | the aligned fragment(s); >1 when the copy is fragmented |

91% of copies are single-fragment. Verified on the full file: the union of
`match_part` spans equals the parent `match` span exactly for **all** copies,
and the short copy id (`ms<N>`) is unique genome-wide.

## Feature selection: match_part, with the copy id in hit_id

The `match` level sums 310.6 Mb; the fragments sum 280.3 Mb. The 30.8 Mb
difference is **gap interior** in fragmented copies — bases between fragments
that REPET does not itself claim as aligned repeat (gaps up to >10 kb;
38,191 sibling fragment pairs overlap, which merged coverage absorbs).
Emitting `match` would paint those gaps as repeat and manufacture support at
bases with no evidence. So the converter emits one BED row per `match_part`
and links fragments of a copy through `hit_id` = the parent's `ms<N>` id —
exactly what BED16 col 16 is for. The copy footprint is recoverable by
grouping on `hit_id`; the reverse (removing gap bases from a match-level
ingest) is not.

## Divergence: fragment Identity, not match AlignIdentity

`match_part.Identity` is the alignment identity of the fragment against the
TEdenovo library consensus → `perc_div = 100 − Identity` is genuine
consensus divergence (same quantity as RepeatMasker's column; EDTA
precedent), so `divergence_is_consensus` is true and REPET feeds the
divergence track. The match-level `AlignIdentity` is **not** used: it reads
`0.00` on 196,332 copies whose fragments carry real identities — a
placeholder in merged features, not a measurement. 11,408 fragments
(1.0%) carry no `Identity` attribute and get `perc_div = NA`.
`SW_score`, `perc_del`, `perc_ins`, `query_left` are never reported: NA
(hence `rm_fields=mixed`).

## Classification: Wicker codes, including compounds

Class comes from `Wcode:<code>` at the head of the parent match's
`TargetDescription`, emitted verbatim. 36 labels observed:

- `NA` — 485,905 copies (48%): unclassified consensus → bare `repeat`
  (existence support, abstains from class votes). Like RM2's 39% Unknown,
  the abstention policy is load-bearing here, not a corner case.
- 13 plain codes (`RIX`, `DTX`, `DXX-MITE`, `RSX`, `RXX-LARD`, `RLX`, …).
  Five already existed in `class_map.tsv` scoped to `edta` and were widened
  to `*` after verifying identical semantics; `RPX`/`RYX`/`DHX`/`DYX` are new
  `*` rules; `DXX-MITE`/`RXX-LARD`/`RXX-TRIM` are repet-scoped (LARD/TRIM →
  `repeat:TE:ClassI:LTR`: non-autonomous LTR-derivatives, order asserted,
  superfamily not; DXX-MITE → TIR, same reasoning as the edta MITE rule).
- **21 compound codes** (`RIX|DTX`, `RLX|RYX`, … ~61k copies): the consensus
  matched references from more than one superfamily — chimeric or ambiguous.
  Each maps **exactly** to the deepest hierarchy level its component codes
  share (`RIX|RSX` → `repeat:TE:ClassI`; `RIX|DTX` → `repeat:TE`). Mapping to
  the first code would assert a superfamily REPET declined to choose.

Compound codes required a precedence fix in `ClassMap.lookup`: an exact rule
matching the full label now outranks any prefix rule (previously the deeper
`RLX*` prefix rule would have beaten the exact `RLX|RYX` rule and asserted
LTR). The fix is restricted to the full label — an exact hit on the bare-code
fallback candidate (the `Class/CODE` split) does not get the boost, otherwise
`SINE/tRNA` → split candidate `tRNA` → `repeat:multigene:tRNA` would have
beaten the intended `SINE/tRNA*` rule. Verified: zero lookup changes across
all 185 (label, tool) pairs observed in the current 9-tool inputs.

## Scope: structural_te

This GFF is TEannot's **TE annotation only** — REPET annotates SSR/tandem
repeats in a separate output not ingested here. Declaring `general_homology`
would count REPET as a dissenting-capable tool at tandem loci it never looks
at (the exact bug documented in DESIGN_NOTES.md). If the SSR output is added
later it should be a separate manifest row, like the FasTAN/TRF/Satellome
family. `structural_te` is used in its eligibility sense — "TE classes only" —
even though the evidence type is homology to a de novo library
(evidence reads "homology" in the per-tool track, which is correct).

## Sequence names

Seqids are ENA FASTA headers with the whitespace stripped
(`CATOHO010000001.1Gobiusnigergenomeassembly...`). The converter extracts the
leading INSDC accession (`^[A-Z]{2,6}\d+\.\d+`) when the seqid is not already
clean; 275/275 distinct seqids mapped to `chrom.sizes`, zero unknown. The 23
missing scaffolds (298 in the assembly) simply have no REPET hits.

## Prototype accommodations vs. production invariants

Real invariants (belong in production): match_part-level ingest with copy id
in `hit_id`; fragment `Identity` (not `AlignIdentity`) as the divergence
source; verbatim Wcode emission incl. `NA` and compounds; compound → LCA
mapping via exact rules; `structural_te` scope; round-trip accounting.

Prototype accommodations (flag as temporary): the seqid accession extraction
exists because this file was built from ENA headers with spaces stripped —
production REPET runs on properly named assemblies should not need it. The
`OtherTargets` attribute (11,247 matches) is not carried.

## Round-trip (observed, full file)

```
matches:       1,016,780
match_parts:   1,128,771
rows written:  1,128,771   (zero dropped, zero orphans)
BED span:      280,287,568 bp  == sum of match_part spans, byte-identical
seqids:        275/275 mapped after accession extraction
copy ids:      1,016,780 distinct == number of matches
```
