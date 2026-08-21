# Satellome

`scope=tandem` · `rm_fields=no` · `ran=yes` (GCA_951799975.1)

## Source

Supplied locally as native BED5 (`GCA_951799975.1_fGobNig1.1_genomic.10kb.bed`,
26 KB). Not yet on GenomeArk; a prototype accommodation like Pantera and
LTRDeNovo — see `docs/INPUT_FORMAT.md` §6.

## Native format

```
chrom  chromStart  chromEnd  family  length
```

Five columns, no header. Column 5 is redundant — verified equal to
`chromEnd − chromStart` on all 573 rows. Coordinates are already 0-based
half-open BED; no shift is applied. Contig naming is mixed (`OX*` INSDC
accessions and `CATOHO*` WGS names in the same file); all 84 contigs resolve
via `chromAlias`, 0 unknown.

## Observed content

573 arrays, 18.80 Mb summed span, 168 families. Family names are
assembly-prefixed serials (`fGobNig19A`). Length range 10,005 – 575,115 bp
(median 15,000).

**The 10 kb floor is a property of this file, not of the tool.** The filename's
`10kb` suffix and the observed minimum (10,005 bp) say shorter arrays were
filtered before delivery. Absence of a Satellome call is therefore not evidence
of absence of satellite — FasTAN calls tandem arrays well below this floor.
Recorded in the manifest notes.

## Conversion decisions

- **Class**: every row emits `repeat_class_family=Satellite` →
  `repeat:tandem:satellite` (existing `Satellite*` rule, `tool=*`, high
  confidence). Unlike FasTAN's opt-in period heuristic this is not a converter
  guess: detecting satellite arrays is the tool's entire method, so the class
  is the tool's own assertion.
- **Scope `tandem`**: same reasoning as FasTAN. A satellite call is a real
  class that would conflict with `repeat:TE` at the hierarchy root; the scope
  gate keeps Satellome out of the eligibility denominator on TE loci and
  admits its vote only where the unrestricted tools already place the locus
  in `repeat:tandem`.
- **No scores**: the native output carries no score, identity or divergence;
  all RepeatMasker-style columns are `NA`, `rm_fields=no`. Verified the
  divergence track's `basesCovered` is unchanged by this tool.
- **Strand `.`**: satellite arrays carry no meaningful orientation here.
- **`hit_id`**: `<family>_<serial>` — unique per assembly, groupable by family.

Verified round-trip: 573/573 records, zero coordinate or name mismatches,
identical 18,802,940 bp summed span, all hit_ids unique.

## Interaction with the classify/detect label

Satellome triggered a refinement of the `nClassify` numerator introduced the
same day: at a locus three TE tools call `LINE:R2` and Satellome calls
`satellite`, the naive count of "tools asserting any class" rendered
`LINE 5/5`. The numerator now counts only tools whose classification is
path-compatible with the displayed consensus, so that locus reads `LINE 3/5`
with `classified by rm2,edta,pantera only` in the hover. On a disputed locus
(bare-`repeat` consensus) the numerator falls back to all real classifiers, so
`Class disputed 3/3` reads "three classifiers dispute".
