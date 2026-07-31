# An integrated TE track for the Vertebrate Genomes Project

> ## ⚠️ Under active development — not for research use
>
> **This repository is a working prototype. Do not use these tracks or this
> pipeline as a source of repeat annotation for research.**
>
> - Output has **not been benchmarked** against a curated truth set. We do not
>   currently know its error rate.
> - The class vocabulary is **provisional** (`dfam-rm-2026.1-draft`) and will be
>   revised when the unified Dfam/Repbase classification is released. Class
>   labels may change.
> - Formats, field names and semantics **will change without notice** and
>   without migration paths.
> - Only **one assembly** has been processed end to end, with **three of four**
>   intended tools.
> - This is **not** the official VGP repeat annotation. It is a transparency
>   layer showing what existing tools currently say, published while the
>   official pipeline is being developed.
>
> If you need repeat annotation for analysis today, run a tool you can cite and
> validate yourself. If you want to help us make this trustworthy, see
> [Contributing](#contributing).

---

## Purpose

Many repeat-annotation tools are run on Vertebrate Genomes Project assemblies,
and they disagree — about where repeats are, how far they extend, and what they
are. This project publishes those disagreements **as evidence rather than as a
verdict**.

Nothing is filtered on the basis of agreement. A locus called by one tool alone
is displayed as prominently as one called by all of them. Where tools disagree,
the disagreement is encoded in the display rather than resolved away.

This is deliberately *not* a merged annotation. Merging requires deciding whose
boundary and whose classification wins, and that decision is exactly what the
official pipeline is being built to make. Until then, the honest product is the
evidence.

## What it produces

A UCSC-style assembly track hub (also loadable in IGV) containing:

| Track | What it shows |
|---|---|
| **repeatSummary** | consensus call per element — class, how many tools support it, how deeply they agree, and thick/thin geometry showing where they agree on boundaries |
| **repeatSupport**, **repeatSupportFrac**, **repeatDivergence** | per-base signals at full resolution, not display-merged |
| **toolUnique** | intervals called by exactly one tool — where the tools disagree most |
| **repeat_\<tool\>** | each tool's full unmodified output, with its own class label preserved verbatim |

![Expected view at three test loci](report/igv_expected_view.png)

*Three contrasting loci. Thick block = the core every eligible tool agreed on;
thin = full extent. Top: all three tools agree to superfamily. Middle: all three
call a repeat but contradict each other on class (grey), so the agreed core is a
57 bp sliver of a 13 kb feature. Bottom: a single-tool call.*

## Status

Processed end to end on one assembly: `GCA_951799975.1` (*Gobius niger*, black
goby), with RepeatModeler2, EDTA and Pantera. fastLTR is configured but has not
been run, and ships as an empty placeholder track.

| | |
|---|---|
| Assembly | 870.6 Mb |
| Covered by ≥1 tool | 431.3 Mb (49.5%) |
| Summary features | 1,919,967 |
| Single-tool intervals | 1,519,062 (91.1 Mb) |
| Class labels observed / mapped | 135 / 135 |
| Validation | `hubCheck` clean |

These numbers describe **what the tools said**, not what is true. See
[Known limitations](docs/SPECIFICATION.md#9-known-limitations).

## Quickstart

Requires Python ≥3.10 with pandas and numpy, plus the UCSC utilities
`bedToBigBed`, `bedGraphToBigWig`, `bigBedInfo` and `hubCheck` on your `PATH`
(from [hgdownload](https://hgdownload.soe.ucsc.edu/admin/exe/)).

Build a hub from the bundled 1.4 Mb test slice — takes about a second:

```bash
python -m vgptrack.cli build \
    --assembly TESTASM \
    --sizes tests/data/test.chrom.sizes \
    --bed rm2=tests/data/rm2.bed.gz \
    --bed edta=tests/data/edta.bed.gz \
    --bed pantera=tests/data/pantera.bed.gz \
    --out hub --description "test slice"
```

Expected output: `9,350 segments over 1,025,484 bp`, `3,739 display features`,
`hubCheck clean`.

On a real assembly, add `--alias` (tools use different sequence-naming
authorities, so this is effectively required) and `--twobit`:

```bash
python -m vgptrack.cli build \
    --assembly GCA_951799975.1 \
    --sizes GCA_951799975.1.chrom.sizes \
    --alias GCA_951799975.1.chromAlias.txt \
    --bed rm2=rm2.bed --bed edta=edta.bed --bed pantera=pantera.bed \
    --out hub --description "Gobius niger (black goby)" \
    --twobit https://hgdownload.soe.ucsc.edu/hubs/GCA/951/799/975/GCA_951799975.1/GCA_951799975.1.2bit
```

Adding a tool is a `--bed` argument plus a row in `config/tools.tsv`. Tools in
the manifest without a `--bed` become empty placeholder tracks automatically.

### Viewing in IGV

```bash
python -m vgptrack.cli session --hub hub/GCA_951799975.1 --out igv_session.xml
```

writes an IGV session plus a genome descriptor that streams the reference from
UCSC GenArk, so there is nothing to download. See [docs/TESTING.md](docs/TESTING.md).

## How it works

1. **Ingest** — each tool provides a standardized BED16
   ([INPUT_FORMAT.md](docs/INPUT_FORMAT.md)); coordinates
   validated against `chrom.sizes`, sequence names reconciled via chromAlias.
2. **Harmonize** — each tool's class label is mapped onto a shared hierarchy via
   `config/class_map.tsv`. The tool's original label is always preserved.
3. **Segment** — per-base state, run-length encoded. Support counts **distinct
   tools**, carried as a bitmask, so one tool's overlapping calls cannot inflate
   it.
4. **Aggregate into elements** — boundary disagreement is absorbed into BED12
   thick/thin geometry rather than merged away.
5. **Build tracks** — bigBed/bigWig via the UCSC utilities, validated with
   `hubCheck`.

Full detail in [docs/SPECIFICATION.md](docs/SPECIFICATION.md).

## Repository layout

```
vgptrack/          pipeline package
  ingest.py        BED16 reading, validation, class harmonization
  segment.py       per-base segmentation and agreement resolution
  summary.py       element aggregation and the summary track
  hub.py           per-tool tracks, discordance track, hub files
  bigfiles.py      UCSC utility wrappers
  igvsession.py    IGV session/genome writers
  cli.py           command-line driver
config/            tool manifest, class vocabulary, colour palette  (data, not code)
docs/              input format (normative), specification, testing guide, design notes
tests/data/        1.4 Mb slice for a runnable end-to-end example
scripts/           tool-output → BED16 converters
report/, qc/       concordance figures and ingest QC from the goby run
```

`config/class_map.tsv` is **data, not code**. It carries a version header and a
provisional-status marker; re-running against a new version rebuilds every track
and nothing in the pipeline hardcodes a class name.

## Contributing

Most useful right now:

- **A curated truth set** for any vertebrate assembly, so the output can finally
  be benchmarked rather than merely described.
- **Vocabulary review** — `config/class_map.tsv` maps 135 observed labels onto a
  Wicker-style hierarchy. Mapping errors silently become false disagreement.
- **Converters** for tools not yet supported (`scripts/`), producing the
  [BED16 input format](docs/INPUT_FORMAT.md). EDTA and fastLTR are wanted.
- **Design critique** — the open questions are collected in
  [docs/DESIGN_NOTES.md](docs/DESIGN_NOTES.md).

## Citing

Please **don't** cite this yet. There is nothing benchmarked to cite. If it is
useful in preparing work, link to this repository and state the commit.

## License

GPL-3.0. See [LICENSE](LICENSE).
