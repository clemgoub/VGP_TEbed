# Tandem subclass discriminability (GCA_951799975.1)

Question raised in slide review (2026-08-22): are `satellite` / `simple` /
`lowcomplexity` discriminant categories *as annotated*, given that TRF and
FasTAN cannot distinguish them (they vote bare `tandem`), Satellome only ever
says satellite, and only rm2/pantera assert all three from homology?

## Measurements (length-weighted, per-base painting of raw BED16 inputs)

Confirmation = fraction of the subclass's bp overlapped by an independent
tandem-array finder (TRF / FasTAN); period statistics are TRF/FasTAN periods
at the confirmed bases.

| rm2 subclass | bp | TRF confirm | FasTAN confirm | Satellome overlap | period at confirmed bp |
|---|---|---|---|---|---|
| Satellite | 3.43 Mb | 15.7% | 20.7% | 0.3% | median 84 bp; 23% >= 100 (TRF) / 22% (FasTAN) |
| Simple_repeat | 14.10 Mb | 74.3% | 74.0% | 14.8% | median 9 bp; 91% <= 49 (TRF) / 95% (FasTAN) |
| Low_complexity | 2.20 Mb | 60.5% | 63.2% | 21.7% | median 15 bp; 81% <= 49 (TRF) / 93% (FasTAN) |

Reference: Satellome arrays (18.8 Mb, all >= 10 kb) are 55%+ TRF-covered with
median period 31 bp, 31% >= 100 bp — the period signature homology "Satellite"
shows only at its confirmed minority of bases.

Composition of the Satellome footprint under rm2: 49.3% no rm2 call at all,
37.0% annotated as **TE families**, 11.1% Simple_repeat, 2.5% Low_complexity,
**0.1% Satellite**. The two satellite sources are near-disjoint: homology
"Satellite" and array-evidence satellite are different populations.

## Verdicts

1. **`Satellite` from homology is not a reliable subclass.** 79–84% of rm2
   Satellite bp is unconfirmed by either array finder, and it misses the
   actual large arrays almost entirely (0.3% Satellome overlap). Kept as
   existence evidence (`repeat:tandem`), subclass withheld.
2. **`Satellite` from Satellome stays.** Array evidence (>= 10 kb, TRF-style
   support) is what the subclass should mean. Rule scoped `tool=satellome`,
   `repeat:tandem:satellite`, high confidence.
3. **`Low_complexity` merges into `simple`.** Its period space at confirmed
   bases (median 15, 81-93% <= 49) is inside Simple_repeat's; it is *not*
   homopolymer-dominated (period 1–2 = 6% of confirmed bp), so it does not
   even match its own name's implication. One subclass fewer, nothing lost.
4. **`Simple_repeat` keeps subclass depth** — 74% independent confirmation,
   coherent low-period signature.

## Implementation (config/class_map.tsv, 2026-08-22)

```
Satellite*  satellome  repeat:tandem:satellite   high
Satellite*  *          repeat:tandem             medium   (demoted)
Simple_repeat*  *      repeat:tandem:simple      high
Low_complexity* *      repeat:tandem:simple      medium   (merged)
```

Rule resolution check: deeper-path-wins precedes tool-specificity in
`ClassMap.lookup`, so the satellome-scoped depth-3 rule beats the demoted
wildcard depth-2 rule for Satellome, and every other tool falls through to
the wildcard. Verified on all six affected (label, tool) pairs before the
rebuild.

Effect on the canonical hierarchy: `repeat:tandem:{satellite,simple}` are the
only tandem subclasses in play; `lowcomplexity` no longer occurs. FasTAN/TRF
continue to vote bare `tandem` (their abstention below kind is unchanged and
correct).

Caveats: measured on one assembly (a goby); the demotion is written as
assembly-independent because the mechanism (homology "Satellite" labels from a
de novo library vs physical array evidence) is not goby-specific, but the
numbers should be re-measured when more assemblies are in the hub. Regenerate
with the session notebook's painting pass or fold into a QC script in phase 2.

## Effect on the consensus build

See the rebuild comparison in the commit that carries this file: fixture
9,350 -> 9,348 segments (vote changes shift two boundaries); genome-wide
numbers in the commit message.
