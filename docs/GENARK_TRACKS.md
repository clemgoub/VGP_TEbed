# GenArk-derived tracks: TRF (simpleRepeat) and WindowMasker (WM + SDust)

Unlike the locally-supplied inputs (Pantera, LTRDeNovo, Satellome), both of
these download straight from the assembly's public GenArk directory — they are
the first inputs that are fully rebuildable for **any** VGP assembly with no
local files:

```
BASE=https://hgdownload.soe.ucsc.edu/hubs/GCA/951/799/975/GCA_951799975.1/bbi
curl -O $BASE/GCA_951799975.1_fGobNig1.1.simpleRepeat.bb
curl -O $BASE/GCA_951799975.1_fGobNig1.1.windowMasker.bb
./bin/bigBedToBed simpleRepeat.bb stdout | python scripts/simplerepeat2bed.py - -o inputs/trf.bed
./bin/bigBedToBed windowMasker.bb stdout | python scripts/windowmasker2bed.py - -o inputs/windowmasker.bed
```

## TRF — `trf` · `scope=tandem` · `rm_fields=divergence_only`

GenArk's `simpleRepeat` track is Tandem Repeats Finder output: 494,440 arrays,
166.6 Mb summed span, the full TRF table per row (period, copyNum,
consensusSize, perMatch, perIndel, score, base composition, entropy, unit
sequence).

**What it annotates**: tandem arrays of any period — 12% of rows are
microsatellite-range (period ≤ 6), 88% longer, up to multi-kb satellite units.

Conversion decisions:

- **Class `tandem` → `repeat:tandem` for every row.** Asserting
  `Simple_repeat` vs `Satellite` from period would be a converter guess, not a
  TRF assertion — identical reasoning to `fastan2bed.py`, whose
  `--classify-period` is opt-in. TRF agrees with FasTAN/Satellome at the
  tandem level and abstains deeper.
- `name = p<period>_x<copyNum>_<unit≤24bp>` — period and copy number reach
  every mouseover.
- `SW_score` = raw TRF alignment score. `perc_div = 100 − perMatch` is
  **unit-to-unit identity within the array**, not consensus divergence —
  hence `rm_fields=divergence_only`, excluded from the divergence track.
- `perIndel` maps to neither `perc_del` nor `perc_ins` alone; dropped (NA).
- Not carried: base composition, entropy, consensusSize, full unit sequence —
  all rederivable from the public track.

## WindowMasker — `windowmasker` · `scope=general_homology` · `rm_fields=no`

GenArk's `windowMasker` track ("WM + SDust"): 4,359,169 intervals, 352.2 Mb —
40% of the assembly. bigBed 3: bare intervals, no name, score or class.

**What it annotates**: windows of over-represented k-mers (WindowMasker
proper) plus low-complexity sequence (the SDust pass), merged
indistinguishably into one track. It asserts *repetitive*, nothing finer.

Conversion decisions:

- **Every row emits `Unknown` → bare `repeat`**: full existence support,
  abstention from every classification vote. WindowMasker is a pure detector —
  the class-map's UNINFORMATIVE rule is exactly its semantics. It can never
  create or break a classification conflict; it raises `n_support` (and, being
  unrestricted, `n_eligible`) wherever it overlaps.
- `scope=general_homology` because a k-mer masker can flag any repetitive
  sequence — TE, tandem or otherwise. Note the evidence is k-mer statistics,
  neither homology nor structure.
- WM-vs-SDust provenance is not recoverable from the public track; if it ever
  matters, WindowMasker must be re-run locally with separate outputs.

## Why these two matter beyond coverage

TRF is the third independent tandem detector (FasTAN, Satellome, TRF — three
different algorithms), so tandem loci can now reach three-classifier genuine
agreement. WindowMasker is the first tool with no classification opinion at
all, which exercises the abstention path at scale: it should raise support
ceilings without ever appearing in a conflict detail.

Both verified on the completed nine-tool build (12,110,870 segments over
508.7 Mb, hubCheck clean):

- 64,215 summary features (20.7 Mb) are `repeat:tandem*` with ≥3 classifying
  tools — e.g. `Satellite 3/5 | classified by rm2,fastan,trf only`.
- Zero conflict rows contain a WindowMasker class assertion; it appears in
  `perToolClass` only as `no class`, while raising `n_support` (max 8).
- The divergence track's coverage moved from 427,120,832 to 427,054,051 bp —
  neither new tool reports divergence; the alias-resolved union of
  divergence-reporting input intervals is 427,037,542 bp, so the finer
  segmentation reduced boundary smear of `mean_div` and moved coverage
  toward that ground truth (overstatement 83.3 kb → 16.5 kb).
