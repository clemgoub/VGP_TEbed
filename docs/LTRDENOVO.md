# LTRDeNovo (NGSEP)

`scope=ltr_only` · `rm_fields=no` · `ran=yes` (GCA_951799975.1)

## Source

Supplied locally as a native GFF3 (`GCA_951799975.1_LTRDeNovo.gff.gz`, 81 KB),
not yet on GenomeArk. `source=NGSEP` on every feature line.

## Feature selection — the load-bearing decision

The GFF3 describes each insertion at two nested levels, plus sub-parts for
structurally-detected elements:

| type | n | what it spans |
|---|---|---|
| `repeat_region` | 3,104 | the full element, LTR to LTR |
| `transposable_element` | 3,104 | the **internal domain only**, between the LTRs |
| `five_prime_LTR` | 228 | the 5' LTR (structural calls only) |
| `three_prime_LTR` | 228 | the 3' LTR (structural calls only) |
| `target_site_duplication` | 118 | the TSD (structural calls only) |

`scripts/ltrdenovogff2bed.py` keeps `repeat_region`. For the 2,876 homology
calls the two levels are coordinate-identical, so the choice is invisible; for
the 228 structural calls `transposable_element` is strictly *inside*
`repeat_region`. TE_12 is representative:

```
repeat_region          2500121-2504726   (4,606 bp)
  five_prime_LTR       2500121-2500360
  transposable_element 2500361-2504487   (4,127 bp)
  three_prime_LTR      2504488-2504726
```

Taking `transposable_element` would leave both LTRs unannotated — for an
LTR-specialist tool, the wrong half of the element to drop. Keeping only one
level also avoids double-counting the same insertion, which would inflate both
coverage and tool support.

## No divergence

Verified across all 6,782 feature lines: GFF3 column 6 (score) and column 8
(phase) are `.` throughout, and the only attributes present are `ID`, `Parent`,
`Ontology_term`, `classification`, `method` and `tsd`. There is no identity,
score, or divergence field of any kind.

Every record is therefore emitted with `perc_div=NA` and the tool is registered
`rm_fields=no`, so it contributes to existence and classification calls but
never to the divergence signal. This is the same treatment fastLTR would have
received had GenomeArk not published a RepeatMasker re-annotation for it.

## Classification and the class-map change

`classification` uses Wicker et al. 2007 codes under an `LTR/` prefix:

| label | n | canonical path |
|---|---|---|
| `LTR/RLG` | 1,133 | `repeat:TE:ClassI:LTR:Gypsy` |
| `LTR/RLC` | 1,109 | `repeat:TE:ClassI:LTR:Copia` |
| `LTR/Unknown` | 861 | `repeat:TE:ClassI:LTR` |
| `LTR/RLR` | 1 | `repeat:TE:ClassI:LTR:ERV` |

`config/class_map.tsv` already carried Wicker rules, but scoped `tool=edta`.
Left alone, all 2,243 `RLG`/`RLC`/`RLR` calls would have resolved to bare
`repeat`, silently discarding the superfamily on 72% of this tool's output.

All six LTR codes — `RLC`, `RLG`, `RLB`, `RLR`, `RLE`, `RLX` — were widened to
`tool=*`. Only `RLG`/`RLC`/`RLR` were actually observed in this GFF3; the other
three were widened with them because `RL*` is one closed, published series in
Wicker et al. 2007, and splitting it by tool would mean an unobserved sibling
code silently degrades to bare `repeat` the first time some tool emits it —
exactly the failure this change fixes.

The non-LTR series (`DT*`, `RI*`, `RS*`) stay `edta`-scoped. They are larger and
their per-code semantics vary more between tools, so widening them would assert
cross-vocabulary agreement that has not been checked against real output. EDTA's
own resolution of `DTA`, `DTT`, `RIL`, `RIR`, `RSX`, `RST` and `RPP` was
re-tested after the change and is unaffected.

## Method attribute

`method=homology` (2,876) or `method=structural` (228) is preserved as a suffix
on the BED `name` field — `TE_12_structural` — so the two detection modes stay
distinguishable in the browser without adding a column outside the BED16
contract. The bare `ID` goes in `hit_id`.

## Observed content

3,104 elements, 7.07 Mb, 66 of 298 contigs, lengths 210 bp – 19,595 bp, real
`OX*` contig names, zero coordinate faults against `chrom.sizes`. 123 adjacent
`repeat_region` pairs overlap, which the segmenter handles as it does for any
tool.

## Contribution, and a caveat worth reading

After merging its 1.416 Mb of self-overlap (nested element calls), LTRDeNovo
covers 5.654 Mb, of which only 0.145 Mb is called by no other tool. 97.4% of its
territory is corroborated by at least one other tool.

Against the other LTR specialist it is strikingly non-redundant: fastLTR and
LTRDeNovo intersect over just 1.506 Mb, a Jaccard of 0.102, and 73.4% of
LTRDeNovo's bp are not called by fastLTR at all. Two tools nominally solving the
same problem overlap less than any other pair in the manifest.

That non-redundancy comes with a caveat. Only **37.3%** of LTRDeNovo's territory
is classified `LTR*` by the summary's independently derived consensus, against
**98.1%** for fastLTR. Splitting by detection method shows the effect is not
uniform:

| method | n | bp overlapped | summary says LTR |
|---|---|---|---|
| `structural` | 228 | 1.605 Mb | **85.4%** |
| `homology` | 2,876 | 5.465 Mb | **13.7%** |

The 228 structural calls behave like a well-specified LTR detector. The 2,876
homology calls largely land on territory the other five tools classify as hAT
(0.78 Mb), Jockey LINEs (0.56 Mb), TIR elements (0.54 Mb) and tRNA SINEs
(0.31 Mb) — i.e. not LTR retrotransposons. Structural calls are also 3.7×
longer (median 6,120 bp vs 1,269 bp), consistent with full-length elements
versus short fragmentary matches.

Two readings are possible and this pipeline cannot distinguish them: either the
homology stage is over-calling LTR on non-LTR repeats, or the other five tools
are misclassifying genuinely LTR-derived sequence. What the track does is make
the disagreement visible per-locus rather than averaging it away. Since
`method` is preserved in the BED `name` field, the two modes can be told apart
in the browser — worth doing before treating a homology call as an LTR
annotation.

## Conversion

```bash
python scripts/ltrdenovogff2bed.py GCA_951799975.1_LTRDeNovo.gff.gz -o inputs/ltrdenovo.bed
```

Verified round-trip: 3,104 records out, all 3,104 source IDs present, zero
coordinate mismatches against the GFF3 (1-based inclusive → 0-based half-open),
and an identical 7.070 Mb total span.