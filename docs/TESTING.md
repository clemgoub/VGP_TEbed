# Testing the prototype hub

Two options, and they answer different questions. **Start with IGV** — it takes
two minutes and tests the data. Use a local UCSC browser only when you need to
test the *hub* rather than the tracks.

---

## Why IGV first

The decisive fact: UCSC already hosts this assembly as a GenArk hub, and its
sequence names and lengths are **identical to ours** — all 298 sequences,
verified name-and-length against
`hgdownload.soe.ucsc.edu/hubs/GCA/951/799/975/GCA_951799975.1/`.

So there is nothing to build. IGV can stream the reference from UCSC and read
our local bigBeds against it.

| | IGV | local UCSC browser |
|---|---|---|
| setup | 2 min, nothing to install | Docker/`hgMirror`, MariaDB, a few hours |
| tests the tracks | yes | yes |
| tests the *hub* (`trackDb`, superTrack, composite) | no | yes |
| mouseOver / `mouseOverField` | no — IGV shows all fields in a popup | yes |
| `itemRgb`, thick/thin BED12 | yes | yes |

IGV **cannot** load `hub/hub.txt`: its hub support expects `useOneFile on`
single-file hubs, and ours is deliberately multi-file because that is what the
UCSC contributed-tracks layout requires. This is not a defect in either — load
the bigBeds directly instead, which is what the session file below does.

---

## Option 1 — IGV (recommended first pass)

```
File > Genome > Load Genome from File...   ->  GCA_951799975.1.genome.json
File > Open Session...                     ->  igv_session.xml
```

`GCA_951799975.1.genome.json` streams the 2bit and chromAlias from UCSC — no
252 MB download. `igv_session.xml` loads all 8 tracks with the per-tool colours,
bigWig ranges (support 0–3, divergence 0–40 %) and display modes already set,
and opens at test locus A.

Both files use absolute paths to the hub, so keep the hub where it is or
regenerate the session after moving it.

### Test loci

Paste into IGV's location box. Each was chosen to exercise a different
behaviour, and `report/igv_expected_view.png` shows what each should look like.

| | locus | what you should see |
|---|---|---|
| **A** | `OX637595.1:3,614,000-3,620,400` | `DNA 3/3`, all three tools, thick core 4121/5369 bp, divergence 2.6 % |
| **B** | `OX637595.1:4,700,200-4,714,400` | `Class disputed 3/3` — grey, core only 57/13188 bp. All three tools call a repeat; they disagree TE-vs-tandem |
| **C** | `OX637595.1:4,714,800-4,720,000` | `DNA 1/3` — EDTA only, 4.4 kb, no full-support core, divergence 22.6 % |
| **D** | `OX637595.1:3,948,813-3,953,947` | `LTR 3/3`, core 5133/5134 bp — near-perfect boundary agreement, the contrast to B |

### What to check

1. **Colour** tracks the consensus class, and grey means disputed — not
   unclassified.
2. **Thick vs thin.** The thick block is the core all tools agree on; the thin
   flanks are where they disagree about boundaries. At B the block is almost
   entirely thin; at D almost entirely thick. That contrast is the main design
   claim.
3. **Click a feature.** IGV prints all 23 fields. Check `supportingTools`
   matches which per-tool rows are present, and `perToolClass` shows the
   disagreement explicitly at B.
4. **`repeatSupport`** should step down exactly where a per-tool row ends.
5. **Search.** Typing a class name in the location box should jump to a feature
   (`searchIndex name`).

Note IGV ignores `mouseOverField`, so the composed one-line hover is *not*
testable here — that needs Option 2.

---

## Option 2 — local UCSC browser

Worth it only to test hub structure: whether the superTrack groups correctly,
whether the composite's per-tool subtracks appear under it, whether the
mouseOver line reads well, and whether the documentation page renders.

A faster intermediate that needs no local browser: put the hub behind any HTTPS
URL and paste it into **My Data > Track Hubs > My Hubs** on
`genome.ucsc.edu`. UCSC fetches it by URL, so you get the real rendering with
no install. The assembly resolves because GenArk already carries it.

If you want a genuinely local instance, `hubCheck` already passes clean on the
full hub, so the remaining risk is display, not validity.

---

## Regenerating these files

```bash
PYTHONPATH=. python -m vgptrack.cli session \
    --hub hub/GCA_951799975.1 --out igv_session.xml
```
