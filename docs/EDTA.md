# EDTA

`scope=structural_te` · `rm_fields=mixed` · `ran=yes` (GCA_951799975.1)

## Source

GenomeArk publishes EDTA coordinates only as GFF3 — there is no RepeatMasker
`.out` in the bucket:

```
s3://genomeark/downstream_analyses/repeats/systematic_annotations/EDTA-v2.3.2/GCA_951799975.1/GCA_951799975.1.fa.mod.EDTA.TEanno.gff3
```

`scripts/edtagff2bed.py` converts it to the 16-column BED. Despite the
`.fa.mod` in the filename, this file carries real assembly sequence names
(`CATOHO010000001.1`), so no renaming is needed — unlike the fastLTR
RepeatMasker output.

## Feature selection: one record per element

The GFF3 mixes whole elements with their sub-parts. A structurally-detected
LTR insertion appears as six records — a `repeat_region` container holding
`lTSD`/`rTSD` target-site duplications, an `LTRRT` element, and its
`lLTR`/`rLTR` terminal repeats. Ingesting all of them would count one locus up
to three times and inflate both coverage and apparent tool support.

The converter keeps `TE_homo_*` (1,389,367), `TE_struc_*` (3,994) and
`LTRRT_*` (1,173), and drops 5,863 sub-feature records.

## Divergence

EDTA reports `identity` (0–1 against the library consensus) rather than
RepeatMasker's percent divergence. Homology calls convert as
`perc_div = 100 × (1 − identity)`, the same quantity.

**Structural calls are emitted as `NA`.** EDTA does write an `identity` for
them, but it is an LTR-to-LTR identity — a within-element measure of how much
the two terminal repeats have drifted apart, not divergence from a consensus.
Reading the two as the same number would corrupt the mean-divergence signal.
This reproduces the NA pattern of the earlier RepeatMasker-derived input
exactly: 3,994 + 1,173 = **5,167 NA records**, matching to the record.

## Agreement with the previous input

The GFF3-derived BED reproduces the RepeatMasker-derived one closely:

| | GFF3 (v2.3.2) | previous (v2.3.1) |
|---|---|---|
| records | 1,394,534 | 1,394,714 |
| coverage | 372.35 Mb | 372.3 Mb |
| NA divergence | 5,167 | 5,172 |
| labels | 38 | 37 |

The differences are the EDTA version, not the conversion.

## `Simple_repeat` vs `scope=structural_te`

The 38th label is `Simple_repeat` — 5,463 records, 6.35 Mb — which maps to
`repeat:tandem:simple`. The manifest note previously asserted EDTA "emits no
tandem/satellite labels"; that was true of the RepeatMasker-derived input and
is **not** true of the GFF3.

This does not require a scope change. `structural_te` admits
`repeat:TE`, so these votes fall outside EDTA's declared scope and the scope
gate excludes them from classification — which is the intended behaviour, not
a loss: EDTA is not a tandem-repeat finder, and FasTAN covers that call. The
records still contribute to existence and appear on EDTA's own per-tool track.

Verified on the smoke scaffold: 4.96 Mb resolves to `repeat:tandem*`,
29.14 Mb to `repeat:TE*`, and `n_support <= n_eligible` holds everywhere.
