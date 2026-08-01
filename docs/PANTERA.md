# Pantera (panteraGA)

`scope=general_homology` · `rm_fields=yes` · `ran=yes` (GCA_951799975.1)

## Source — and why this one is not yet rebuildable from GenomeArk

panteraGA builds TE libraries from a FastGA alignment of two genomes. GenomeArk
publishes the **library**, not the annotation:

```
s3://genomeark/downstream_analyses/repeats/panteraGA/GCA_951799975.1/GCA_951799975.1.fGobNig-pantera-pass.fa.gz
```

That file is 0.19 MB — 196 consensus sequences, RepeatMasker-style headers
carrying the class (`>hAT_1-fGobNig#DNA/hAT`, `>Gypsy_15-fGobNig#LTR/Gypsy`).
There are no coordinates in it.

The ingestible product is the RepeatMasker run of that library against the
assembly, which produced the 1,213,581 hits used here. **That run is not on
GenomeArk**, so unlike the other four tools, `inputs/pantera.bed` is currently
supplied locally rather than fetched. Note that `panteraGA/` sits directly under
`downstream_analyses/repeats/`, not under `systematic_annotations/` where the
other three tools live; the upload convention is still in flux, so re-check the
prefix before automating a fetch.

To make it rebuildable, either the `.out` gets uploaded alongside the library,
or the pipeline gains a step that runs RepeatMasker with the published library.
The latter is reproducible from public data but costs a full RepeatMasker pass
over an 870 Mb genome.

## Observed content

1,213,581 hits, 229.1 Mb covered, 21 class labels, zero unmapped, zero
malformed coordinates, all contigs present in the assembly's `chrom.sizes`.

The vocabulary is the coarsest of the three homology tools — 40% bare `DNA` and
29% `Unknown`:

| label | hits |
|---|---|
| `DNA` | 490,696 |
| `Unknown` | 352,032 |
| `Simple_repeat` | 262,232 |
| `LINE/L2` | 30,650 |
| `DNA/TcMar-Tc1` | 21,696 |
| `LINE/RTE` | 18,866 |
| (15 more, each < 7,000) | |

`DNA` maps to `repeat:TE:ClassII` and `Unknown` to bare `repeat`, so Pantera
contributes strongly to existence calls and weakly to classification — which is
what the manifest's note records.

## Contribution

Pantera covers 229.1 Mb but only 2.2 Mb that no other tool calls. Its value is
independence, not reach: it has the lowest pairwise Jaccard of the three
homology callers (0.542 vs RepeatModeler2, 0.526 vs EDTA, against 0.705 between
those two), so per bp it supplies the most non-redundant support in the set.

## Ingestion

The BED16 is already in canonical form — real `OX*`/`CATOHO*` contig names, the
16 documented columns, `perc_div` populated — so it needs no converter:

```bash
python -m vgptrack.cli build --assembly GCA_951799975.1 \
    --sizes data/GCA_951799975.1.chrom.sizes \
    --alias data/GCA_951799975.1.chromAlias.txt \
    --bed pantera=inputs/pantera.bed  ...
```

## History

An earlier three-tool build used this same input (identical 1,213,581 / 229.1 Mb
figures). When the raw inputs were later lost from local disk and each tool's
input was rebuilt from GenomeArk, Pantera was the one tool whose annotation is
not published there — so it silently became an empty placeholder track while
fastLTR took its place. Nothing failed; the build simply shipped a zero-feature
track and a support ceiling of 4 instead of 5.

`--allow-missing` now exists because of that: a tool with `ran=yes` in
`config/tools.tsv` and no `--bed` on the command line aborts the build rather
than shipping empty. See [Guard against a forgotten
input](../README.md#guard-against-a-forgotten-input).