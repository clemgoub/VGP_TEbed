# Adding FasTAN to the hub

FasTAN (<https://github.com/thegenemyers/FASTAN>) is a tandem array finder. It
is the first non-TE tool in the pipeline, so it exercises paths the three TE
tools never did. This note covers the conversion, what changed in the pipeline,
and how to rebuild your local hub.

## 1. Convert the native output

FasTAN's BED is 5 columns, no header:

| col | meaning |
|-----|---------|
| 1-3 | chrom, chromStart, chromEnd (0-based half-open) |
| 4   | estimated average unit (period) size, bp |
| 5   | average identity among array units, per mille |

```bash
python scripts/fastan2bed.py fGobNig-tan.bed -o fastan.bed16.bed
```

Verified on `fGobNig-tan.bed`: 217,563 records in, 217,563 out, none dropped,
coordinates and score byte-identical to the input.

### Two deliberate abstentions

**Classification.** By default every record gets the single label `tandem`,
which `config/class_map.tsv` maps to `repeat:tandem`. FasTAN detects arrays; it
does not call them satellite vs simple vs low-complexity. Inferring a subclass
from period size would invent a classification FasTAN never made, and the
summary track would then report class *agreement* or *conflict* on evidence
that does not exist. `--classify-period` opts in to the conventional
`period<=6 -> Simple_repeat, >6 -> Satellite` cut when you want it.

**Library fields.** `SW_score`, `perc_del`, `perc_ins` and the three consensus
coordinates are emitted as `NA`, not `0`. FasTAN does no library alignment, so
those quantities do not exist; `0` would be a value the pipeline cannot
distinguish from a real measurement.

The estimated period is preserved in the feature name (`tandem_p31`), so it is
on the mouseover regardless of the classification flag.

## 2. Divergence is not the same quantity

FasTAN's identity yields a real divergence percentage, but it measures
**unit-to-unit divergence within an array** (array homogeneity) — not
divergence from a library consensus (an age proxy). Averaging the two would
make the genome-wide `repeatDivergence.bw` uninterpretable.

The manifest records this as `rm_fields=divergence_only`, and
`Tool.divergence_is_consensus` gates it: FasTAN is **excluded from the
aggregate divergence track** and its value still appears, self-labelled, on its
own per-tool track:

```
divergence = 5.4% (unit-to-unit within array, NOT vs consensus)
evidence   = tandem array structure
```

Verified: the aggregate track's mean is 12.02% for three tools and 12.00% with
FasTAN added — the residual is segment-boundary refinement, not blending.

## 3. Two bugs this surfaced

Both were latent and would have hit the next tool added regardless of which one
it was.

**Bit-indexed arrays sized by tool count.** `segment_sequence` allocated the
per-tool class array with `len(tools)` rows but indexed it by `tool.bit`. Those
are equal only while the tools that *ran* occupy the lowest bits contiguously.
`fastltr` (registered, not run) sitting below `fastan` broke that and raised
`IndexError`. Now sized by `max(bit)+1`.

**uint8 mask ceiling.** The per-base masks are `uint8` while `ToolSet` permits
64 tools, so a 9th tool would have silently truncated to bit 7 and produced
plausible-but-wrong support counts. Now raises with instructions to widen.

**Manifest row order is the bitmask order — append new tools at the end.**
Inserting mid-table renumbers every tool below it, and masks in already-built
tracks then decode to the wrong tools. This is documented in `tools.tsv`.

Also fixed: `igvsession.py` had a hardcoded tool colour/label table, so a tool
added to `tools.tsv` rendered grey and id-labelled in IGV while being correct in
the hub. It now reads the manifest.

## 4. Rebuild your hub

From the repo root, with the UCSC tools on `PATH`:

```bash
python scripts/fastan2bed.py /path/to/fGobNig-tan.bed -o smoke/fastan.bed

python -m vgptrack.cli build \
  --assembly GCA_951799975.1 \
  --sizes data/GCA_951799975.1.chrom.sizes \
  --alias data/GCA_951799975.1.chromAlias.txt \
  --bed rm2=smoke/rm2.bed \
  --bed edta=smoke/edta.bed \
  --bed pantera=smoke/pantera.bed \
  --bed fastan=smoke/fastan.bed \
  --out hub --work work \
  --description "Gobius niger (black goby)" \
  --twobit https://hgdownload.soe.ucsc.edu/hubs/GCA/951/799/975/GCA_951799975.1/GCA_951799975.1.2bit

python -m vgptrack.cli session --hub hub/GCA_951799975.1 \
  --out hub/igv_session.xml --locus "OX637595.1:1-200000"
```

Then in IGV: **Genome > Load Genome from File** (`GCA_951799975.1.genome.json`),
then **File > Open Session** (`igv_session.xml`).

### What you should see

On the three-chromosome slice (38,893 FasTAN arrays):

| | 3 tools | + FasTAN |
|---|---|---|
| repeat-covered | 84.27 Mb | 91.34 Mb |
| tandem-classified | 3.19 Mb | 11.03 Mb |
| summary features | 370,844 | 378,512 |

FasTAN contributes 7.07 Mb that no TE tool calls. Where it votes, consensus
resolves only to tandem classes, never to a TE class. It is excluded from the
support denominator at all 609,649 TE loci and included at every tandem locus,
and `n_support <= n_eligible` holds everywhere.

## 5. Regression

The three-tool test fixture builds **byte-identical** outputs before and after
all of the above — every `.bb` and `.bw` unchanged. The only difference is an
added `repeat_fastan.bb` placeholder track (0 features, valid bigBed) in
`trackDb.txt`, which is the intended behaviour for a registered tool with no
data on that assembly.
