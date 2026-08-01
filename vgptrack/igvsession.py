"""Write IGV session + genome-descriptor files for a built hub assembly dir.

IGV cannot consume this hub's `hub.txt`: IGV's track-hub support expects
single-file (`useOneFile on`) hubs, whereas the UCSC contributed-tracks layout
this project targets is deliberately multi-file. The tracks themselves are
ordinary bigBed/bigWig and load fine, so instead of degrading the hub to suit
IGV we emit a session that points at the same files directly.

Track display settings are duplicated here rather than parsed out of
trackDb.txt: trackDb's vocabulary (superTrack, composite, parent) has no IGV
equivalent, so a translation would be mostly special cases. Both come from
TRACKS below, which is the single place to edit.
"""

from __future__ import annotations

import json
import os
import xml.etree.ElementTree as ET

# filename, IGV track name, displayMode, (min,max) for bigWigs or None, colour,
# renderer. Renderer is None for feature tracks (IGV picks its own).
TRACKS = [
    ("repeatSummary.bb",     "Repeat Summary (consensus)", "EXPANDED", None,    None,        None),
    ("repeatSupport.bw",     "Tool Support (n tools)",     "FULL",     (0, 3),  "60,60,160", "BAR_CHART"),
    ("repeatSupportFrac.bw", "Support Fraction",           "FULL",     (0, 1),  "60,120,60", "BAR_CHART"),
    ("repeatDivergence.bw",  "Divergence %",               "FULL",     (0, 40), "160,60,60", "HEATMAP"),
    ("toolUnique.bb",        "Single-tool Calls",          "SQUISHED", None,    None,        None),
]


def _heatmap_scale(lim, colour: str) -> str:
    """IGV ContinuousColorScale for divergence, saturated at the LOW end.

    Format is ``ContinuousColorScale;min;max;minColour;maxColour``. The colour
    stops are deliberately inverted relative to the intuitive reading: low
    divergence gets the full colour and high divergence ramps to white.

    Divergence measures decay from the family consensus, so a copy at 2% is a
    recent, probably still-active insertion while one at 35% is an ancient
    relic. Ramping saturation *down* as divergence rises therefore makes the
    young insertions the visually dense ones, which is the signal worth
    spotting when scanning a chromosome.

    IGV's heatmap renderer has no alpha channel, so the fade is produced by
    ramping to white rather than to transparent. That is visually equivalent
    on IGV's white background but will read as solid white, not see-through,
    if a track sits on a coloured backdrop.
    """
    lo, hi = lim
    return f"ContinuousColorScale;{float(lo)};{float(hi)};{colour};255,255,255"

# Per-tool tracks are discovered from the directory, so a hub built with a
# different tool set still gets a complete session. Colours and labels come from
# the manifest rather than a table here: duplicating them meant a tool added to
# config/tools.tsv rendered grey and id-labelled in IGV while being correct in
# the hub. These fall back to the manifest-free defaults if it cannot be read,
# because a session file is a convenience and must never block a build.
_FALLBACK_COLOURS = {
    "rm2": "31,119,180", "edta": "227,26,28",
    "pantera": "51,160,44", "fastltr": "255,127,0",
}
_FALLBACK_LABELS = {
    "rm2": "RepeatModeler2", "edta": "EDTA",
    "pantera": "Pantera", "fastltr": "fastLTR",
}


def _tool_style(tools_tsv="config/tools.tsv"):
    """(colours, labels) keyed by tool_id, read from the manifest."""
    try:
        from .vocab import ToolSet
        ts = ToolSet.load(tools_tsv)
    except Exception:
        return dict(_FALLBACK_COLOURS), dict(_FALLBACK_LABELS)
    return ({t.tool_id: t.color for t in ts},
            {t.tool_id: t.short_label for t in ts})

GENARK = ("https://hgdownload.soe.ucsc.edu/hubs/{a}/{b}/{c}/{d}/{acc}")


def genark_base(accession: str) -> str:
    """GenArk URL for an accession: GCA_951799975.1 -> .../GCA/951/799/975/..."""
    prefix, digits = accession.split("_", 1)
    num = digits.split(".")[0]
    return GENARK.format(a=prefix, b=num[0:3], c=num[3:6], d=num[6:9], acc=accession)


def write_genome_json(accession: str, path: str, name: str = "") -> dict:
    """IGV genome descriptor streaming the reference from UCSC GenArk.

    Only valid when the local tracks use the same sequence names as GenArk --
    the caller is expected to have verified that.
    """
    b = genark_base(accession)
    g = {
        "id": accession,
        "name": name or accession,
        "twoBitURL": f"{b}/{accession}.2bit",
        "chromSizesURL": f"{b}/{accession}.chrom.sizes.txt",
        "chromAliasBbURL": f"{b}/{accession}.chromAlias.bb",
        "wholeGenomeView": False,
    }
    with open(path, "w") as fh:
        json.dump(g, fh, indent=2)
    return g


def _is_empty_bb(path: str) -> bool:
    """True if a bigBed holds zero features.

    Read from the header's itemCount rather than shelling out to bigBedInfo, so
    session writing has no dependency on the Kent tools being installed.
    Returns False on any parse problem: mislabelling a populated track is worse
    than missing an empty one.
    """
    import struct
    try:
        with open(path, "rb") as fh:
            head = fh.read(64)
            if head[:4] == b"\xeb\xf2\x89\x87":
                endian = "<"
            elif head[:4] == b"\x87\x89\xf2\xeb":
                endian = ">"
            else:
                return False
            # itemCount is NOT in the 64-byte header; it is a uint64 written at
            # fullDataOffset (header bytes 16..24). Verified against bigBedInfo.
            full_data = struct.unpack(endian + "Q", head[16:24])[0]
            fh.seek(full_data)
            return struct.unpack(endian + "Q", fh.read(8))[0] == 0
    except Exception:
        return False


def write_session(hubdir: str, path: str, accession: str, locus: str = "",
                  tools_tsv: str = "config/tools.tsv") -> str:
    """IGV session XML referencing every track present in *hubdir*.

    Paths are absolute: IGV resolves session-relative paths against its own
    working directory, not the session file, so relative paths break as soon as
    the session is opened from anywhere else.
    """
    hubdir = os.path.abspath(hubdir)
    colours, labels = _tool_style(tools_tsv)
    present = [t for t in TRACKS if os.path.exists(os.path.join(hubdir, t[0]))]
    for fn in sorted(os.listdir(hubdir)):
        if fn.startswith("repeat_") and fn.endswith(".bb"):
            tool = fn[len("repeat_"):-len(".bb")]
            # A tool that did not run is present as a valid 0-feature bigBed so
            # that every assembly dir has identical filenames. Say so in the
            # track name, otherwise it reads in IGV as "this tool found nothing".
            disp_name = labels.get(tool, tool)
            label = (f"{disp_name} (not run)"
                     if _is_empty_bb(os.path.join(hubdir, fn)) else disp_name)
            present.append((fn, label, "SQUISHED", None, colours.get(tool), None))

    root = ET.Element("Session", genome=accession, version="8")
    if locus:
        root.set("locus", locus)
    res = ET.SubElement(root, "Resources")
    for fn, *_ in present:
        ET.SubElement(res, "Resource", path=os.path.join(hubdir, fn))
    panel = ET.SubElement(root, "Panel", name="DataPanel")
    for fn, name, disp, lim, col, rend in present:
        attrs = {"id": os.path.join(hubdir, fn), "name": name,
                 "displayMode": disp, "visible": "true"}
        if col:
            attrs["color"] = col
        if lim:
            attrs.update(autoScale="false", renderer=rend or "BAR_CHART")
        # A heatmap track needs an explicit colour scale, otherwise IGV falls
        # back to its default blue-white-red diverging scale, which implies a
        # meaningful midpoint that divergence does not have.
        if rend == "HEATMAP" and lim and col:
            attrs["colorScale"] = _heatmap_scale(lim, col)
        tr = ET.SubElement(panel, "Track", **attrs)
        if lim:
            ET.SubElement(tr, "DataRange", minimum=str(lim[0]),
                          maximum=str(lim[1]), type="LINEAR")
    ET.SubElement(root, "Panel", name="FeaturePanel")
    ET.indent(root)
    ET.ElementTree(root).write(path, encoding="UTF-8", xml_declaration=True)
    return path
