# fastLTR

`scope=ltr_only` · `rm_fields=yes` · `ran=yes` (GCA_951799975.1)

## Source

GenomeArk publishes a RepeatMasker re-annotation of the assembly using the
fastLTR library:

```
s3://genomeark/downstream_analyses/repeats/systematic_annotations/FastLTR-RepeatMasker/GCA_951799975.1.fa.mod.out
```

The sibling `FastLTR/` prefix holds only the classified library FASTA
(`fGobNig_all_rounds_classified.fasta`) with no coordinates, so the
RepeatMasker `.out` is the ingestible product.

## Observed content

16,772 hits, 10.71 Mb, 5 class labels, all under `LTR/`:

| label | hits |
|---|---|
| LTR/Gypsy | 7,254 |
| LTR/Unknown | 4,800 |
| LTR/Pao | 3,765 |
| LTR/ERV | 772 |
| LTR/Copia | 181 |

Every label is already covered by `config/class_map.tsv`, and all resolve
under `repeat:TE:ClassI:LTR` — confirming `scope=ltr_only` empirically rather
than by assumption.

## `rm_fields=yes`, not `no`

The manifest previously carried `rm_fields=no`, written when fastLTR was
expected to deliver structural predictions. This file is a RepeatMasker run,
so `perc_div` is genuine divergence from a library consensus and is
comparable with the other homology-based tools. It therefore enters the
`repeatDivergence` mean. Observed range on the smoke scaffold: 0–38.5%.

Had this stayed `no`, fastLTR's divergence would have been silently dropped
from the summary signal.

## Sequence names had to be recovered

The `.out` was produced against EDTA's `.fa.mod`, which rewrote every
sequence name to a short opaque id (`_J0000000`, `_J000001q`, …). No mapping
file is published alongside it, and the names match nothing in the assembly.
Note the filename: `.fa.mod.out`, where RepeatModeler's is `.fa.out`.

`scripts/infer_rm_rename.py` recovers the mapping from the `.out` itself.
RepeatMasker records `end` and `(left)` per hit, and their sum is the query
sequence's full length; that gives an exact length per renamed sequence
without needing the modified FASTA. All 126 sequences produced a
self-consistent length, and all 126 matched a sequence in `chrom.sizes`.

**120 of 126 are unique and safe. 6 are not.** Those 6 have round sizes
(2 kb, 5 kb, 6 kb, 7 kb) shared by up to 20 unplaced scaffolds. Three routes
to disambiguate were tried and none worked:

- the ids look base-62 sequential, but assembly FASTA order is neither
  `chrom.sizes` order nor NCBI assembly-report order (50 of 119 anchors
  violate monotonicity in both), so interpolation is unsound;
- elimination against sequences already claimed by unambiguous matches
  leaves every one of the 6 with multiple free candidates;
- no mapping file exists in the bucket.

The 6 affected sequences carry **7 hits totalling 7,629 bp — 0.04%** of
fastLTR's calls. They are dropped rather than guessed.

```bash
python scripts/infer_rm_rename.py inputs/fastltr.rm.out \
    data/GCA_951799975.1.chrom.sizes -o rename_raw.tsv
grep -v AMBIGUOUS rename_raw.tsv > inputs/fastltr.rename.tsv   # inspect first
python scripts/rmout2bed.py inputs/fastltr.rm.out \
    --rename inputs/fastltr.rename.tsv --drop-unmapped -o inputs/fastltr.bed
```

`infer_rm_rename.py` writes ambiguous rows with an `AMBIGUOUS:` marker so an
unedited mapping fails loudly instead of silently mis-assigning; drop or
resolve them before use. Without `--drop-unmapped`, unmapped records pass
through with their original `_J…` name and fail validation downstream.

**This is not fastLTR's doing** — it is an artefact of running RepeatMasker on
EDTA's modified FASTA. EDTA's own GFF3 for the same assembly carries real
sequence names. If the upstream run can publish the `.fa.mod` name mapping, or
run against the unmodified assembly, this whole step disappears and the 7 hits
are recovered.

## Empty-signal fix

A build in which no tool contributes consensus divergence produces an empty
`repeatDivergence` bedGraph. Kent's `bedGraphToBigWig` aborts on that with
`needLargeMem: trying to allocate 0 bytes`. Since contribTracks requires the
file to exist regardless, `bigfiles.bedgraph_to_bigwig` now routes empty input
to the pybigtools writer and emits a valid 0-interval bigWig.
