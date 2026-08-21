# Auto-generated two-level report

`report/generate_report.py --assembly <ACC> --repo .` emits one self-contained
HTML (`report/<ACC>_report.html`, figures embedded) from a built hub. It needs
only `work/segments.parquet`, `config/tools.tsv` and `data/<ACC>.chrom.sizes`;
the Dfam section appears when the family-evidence tables exist under
`report/data/` and degrades to a notice when they don't.

## Level 2 -- biological summary (confidence-tiered)

Repeat content is reported at four nested tiers, never as one number:

| tier | definition |
|---|---|
| union | any tool called a repeat (upper bound; includes WindowMasker singletons) |
| supported >=2 | two scope-eligible tools overlap |
| **consensus (headline)** | >=2 tools AND `conflict_depth < 0` |
| kind-resolved | consensus AND `agree_depth >= 2` |

Black goby (GCA_951799975.1): 58.4% -> 46.7% -> **31.3%** -> 29.1% of 871 Mb.
Class/order breakdowns and the divergence landscape use the headline tier only,
and the landscape further restricts to bp whose `mean_div` comes from tools with
true consensus divergence (`rm_fields=yes`); the caveat block is generated from
`config/tools.tsv` (scope + rm_fields), so it rewrites itself when the tool set
changes.

`agree_depth` counts agreed hierarchy levels: 1 = `repeat` only, 2 = kind
(TE vs tandem), 3 = class, 4 = order, 5 = superfamily. (An earlier draft
mislabelled these off by one -- verified against `_cascade` in
`vgptrack/segment.py`.)

## Level 1 -- technical (agreement / specificity / Dfam)

Length-weighted 9-tool Jaccard, corroborated-vs-solo coverage, conflict
involvement, and the agreement-depth spectrum, all computed fresh from
`segments.parquet` at generation time.

## Dfam low-hanging-fruit shortlist

Built by the family pipeline (currently the session notebook; to be scripted in
phase 2) into `report/data/`:

- `families.parquet` -- one row per (tool, family): copy stats, consensus-
  coverage profile (from `repeat_start/end/left`), corroboration fractions from
  the segmentation, gates, cluster id.
- `dfam_shortlist.tsv` / `dfam_shortlist_clusters.tsv` -- ranked candidates.

Six explicit gates, each reported separately (no hidden composite):
G1 >=5 near-full-length copies; G2 depth>=3 over >=99% of the consensus (the
Dfam seed-alignment requirement); G3 >=60% of bp on class-agreed segments
(`agree_depth >= 3`); G4 <=30% tandem-tool overlap; G5 median divergence <=25%;
G6 >=50% of bp on TE-consensus loci.

Two decisions that came from user review and are load-bearing:

1. **"Near-full-length" is relative to the tool's own consensus, which may be
   incomplete.** The `cons_completeness` column makes this explicit:
   `len_inconsistent` (implied consensus length varies >2% across hits -- often
   a truncated/chimeric library entry) and `weak_5p`/`weak_3p` (edge depth
   < 25% of interior depth). These flags select the phase-2 families whose
   consensus should be rebuilt from an MSA rather than taken from the library.
2. **The deposition unit is the cross-tool cluster, not the best single
   program's family.** Families are clustered by reciprocal >=50% genomic
   footprint overlap (>=1 kb joint); the shortlist ranks clusters, and phase 2
   derives a NEW seed alignment from the pooled copies of all members instead
   of choosing among programs. Black goby: 863 candidate clusters, 577 of them
   supported by >=2 programs.

**Unknown rescue:** a family labelled Unknown/bare-DNA whose loci sit on
confident consensus segments (>=2 tools, conflict-free, dominant path
`repeat:TE:*` at >=60% of bp) inherits that path in `rescued_class`. 1,056
families on black goby -- free classification improvements for the library.

## Phase 2 (deferred, plan on file)

Library FASTA reconciliation (mmseqs2 clustering x co-annotation
cross-validation; redundancy / fragmentation / satellite-leakage metrics per
program), assembly download, copy extraction, MAFFT MSA, consensus rebuild +
completeness assessment, Stockholm seed packets. Requires the three tool
libraries (`rm2_families.fa`, `edta_TElib.fa`, `pantera_lib.fa`) and the
assembly FASTA.

## Known input quirks (do not re-discover)

- 12 rows across rm2 (3) and pantera (9) have `repeat_start > repeat_end`
  (native RepeatMasker artifacts, all <=46 bp hits); dropped from family
  profiles, reported here.
- Simple-repeat `(MOTIF)n` entries dominate family counts (12.8k of 15.3k rm2
  names); excluded from Dfam gating via `is_simple`.
- EDTA BED16 carries no consensus coordinates (`NA`), so EDTA families get
  copy-level stats only and cannot pass G1/G2 directly -- they participate via
  clusters.
