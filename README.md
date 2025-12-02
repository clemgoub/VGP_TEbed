# An integrated TE track for the Vertebrate Genome Project

/// work in progress, subject to change ///

## Purpose

Development of an integrated track in `.bed` format to summarize TE annotations emanating from multiple tools such as `RepeatModeler2`, `EDTA`, `Pantera`, etc... used in the Vertebrate Genome Project. 
The track will display different level of information and will refer back to each individual tool specifics.

## Concept / Brainstorming

We are looking at annotating each base of a genome. For each base we ask: is it a TE? If yes, what type of TE?
Hierarchical structure of annotation. We have multiple tools, and each tools can have overlaping annotations.
1. The smalest unit is a hit made by one tool, represented by an interval.
2. The same tool can have overlaping hits
3. Mutliple tools can have overlaping hits

For a given base of the genome, we can ask how confident we are that this is a repeat. We can do that by asking how many tool agree that this is a repeat. But what if one tools has 2 or more overlaping hits for that base? Is is further evidence, or less or is it irrelevant (I think it is irrelevant and instead it says that there is redundancy in the library).

## Specification

describe here the format specification of the integrated track (it will be `.bed`)

The first 12 fields are constrained by the bed standard:

1. `chr`: chromosome
2. `start`: feature start  (0-based)
3. `end`: feature end (1-based)
4. `name`: unique instance name
5. `score`: *undefined*
6. `strand`: `+`, `-` or `.` if conflict
7. `thickStart`: The starting position at which the feature is drawn thickly (for example, the start codon in gene displays). When there is no thick part, thickStart and thickEnd are usually set to the chromStart position.
8. `thickEnd`: The ending position at which the feature is drawn thickly (for example the stop codon in gene displays).
9. `itemRgb`: An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
10. `blockCount`: The number of blocks (exons) in the BED line.
11. `blockSizes`: A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
12. `blockStarts`: A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount.

The rest is up for grabs. The current order is irrelevant, the goal is to list the information we need to display.

13. `support`: number of methods supporting a repeat in this region `<int>`
14. `methods`: list of methods supporting a repeat in this regions `<string>` (eg: `RM2;EDTA`)
15. 

## Input format requirement
>The goal is to be able to link each hit to a consensus in the library produced by each tool. This is straightfoward using Repeatmasker output (`.out`), but EDTA outputs are more complex as it blends structural and homology calls. 

RepeatMasker `.out` or custom `.bed` with the following specification:

### Mendatory 

1. `chr`: chromosome
2. `start`: feature start  (0-based)
3. `end`: feature end (1-based)
4. `name`: repeat name (consensus name or instance name for EDTA structural call)
5. `score`: SW from RepeatMasker run or `.`
6. `strand`: `+` or `-`
7. `divergence`: % substitutions in matching region compared to the consensus (or `1 - identity` EDTA ?) or `NULL` (EDTA structural)
8. `repeat_class`: in repeatmasker format 

### non-mendatory (based on absent from EDTA output)

8. `consensus_start`: hit start on consensus (0-based)
9. `consensus_end`: hit end on consensus (1-based)
`% deletion`: % of bases opposite a gap in the query sequence (deleted bp)
`% insertion`: % of bases opposite a gap in the repeat consensus (inserted bp)

## Integration tool

describe how the script that will produce the integrated track work.
