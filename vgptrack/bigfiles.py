"""Backend-agnostic writers for UCSC binary track formats.

Primary backend is the UCSC kent command-line suite, which is required for
autoSql (``-as``) support: the rich mouseover fields in the repeatSummary
track are declared through an .as file, and no pure-Python writer supports
that.  ``pybigtools`` is retained as a degraded fallback that can still
produce plain bigBed/bigWig (no typed extra fields).

Kent binaries for macOSX.arm64 / linux.x86_64 can be fetched with
``fetch_kent_tools()`` from https://hgdownload.soe.ucsc.edu/admin/exe/ .
"""

from __future__ import annotations

import os
import platform
import shutil
import subprocess
import urllib.request
from pathlib import Path

KENT_TOOLS = [
    "bedToBigBed",
    "bedGraphToBigWig",
    "bigBedInfo",
    "bigWigInfo",
    "bedSort",
    "bigBedToBed",
    "bigWigToBedGraph",
    "hubCheck",
]

_UCSC_EXE = "https://hgdownload.soe.ucsc.edu/admin/exe"


class BigFileError(RuntimeError):
    pass


def kent_platform() -> str:
    """UCSC's directory name for the current platform."""
    sysname, machine = platform.system(), platform.machine()
    if sysname == "Darwin":
        return "macOSX.arm64" if machine in ("arm64", "aarch64") else "macOSX.x86_64"
    if sysname == "Linux":
        return "linux.arm64" if machine in ("aarch64", "arm64") else "linux.x86_64"
    raise BigFileError(f"no UCSC kent build known for {sysname}/{machine}")


def find_kent(name: str, bindir: str | Path | None = None) -> Path | None:
    """Locate a kent binary in *bindir*, $VGPTRACK_BIN, or $PATH."""
    for cand in (bindir, os.environ.get("VGPTRACK_BIN"), "bin"):
        if cand:
            p = Path(cand) / name
            if p.is_file() and os.access(p, os.X_OK):
                return p.resolve()
    found = shutil.which(name)
    return Path(found) if found else None


def fetch_kent_tools(bindir: str | Path = "bin", tools=None) -> dict[str, str]:
    """Download the kent binaries needed by this pipeline. Returns name->status."""
    bindir = Path(bindir)
    bindir.mkdir(parents=True, exist_ok=True)
    base = f"{_UCSC_EXE}/{kent_platform()}"
    status = {}
    for t in tools or KENT_TOOLS:
        dest = bindir / t
        if dest.is_file() and os.access(dest, os.X_OK):
            status[t] = "present"
            continue
        try:
            urllib.request.urlretrieve(f"{base}/{t}", dest)
            dest.chmod(0o755)
            status[t] = "downloaded"
        except Exception as exc:  # noqa: BLE001
            status[t] = f"failed: {exc}"
    return status


def have_kent(bindir=None) -> bool:
    return all(find_kent(t, bindir) for t in ("bedToBigBed", "bedGraphToBigWig"))


def _run(cmd: list, what: str) -> None:
    proc = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    if proc.returncode != 0:
        raise BigFileError(
            f"{what} failed (exit {proc.returncode})\n"
            f"  cmd: {' '.join(str(c) for c in cmd)}\n"
            f"  stderr: {proc.stderr.strip()[:2000]}"
        )


def sort_bed(src, dest=None, bindir=None):
    """Sort a BED file by chrom then start, as bedToBigBed requires."""
    src = Path(src)
    dest = Path(dest) if dest else src
    exe = find_kent("bedSort", bindir)
    if exe:
        _run([exe, src, dest], "bedSort")
        return dest
    # Portable fallback: LC_ALL=C sort -k1,1 -k2,2n, header lines stripped.
    env = dict(os.environ, LC_ALL="C")
    with open(src) as fh:
        body = [l for l in fh if not l.startswith(("#", "track", "browser"))]
    proc = subprocess.run(
        ["sort", "-k1,1", "-k2,2n"], input="".join(body),
        capture_output=True, text=True, env=env,
    )
    if proc.returncode != 0:
        raise BigFileError(f"sort fallback failed: {proc.stderr[:500]}")
    Path(dest).write_text(proc.stdout)
    return dest


def bed_to_bigbed(bed, chrom_sizes, out, as_file=None, bed_type=None,
                  extra_index=None, bindir=None, sort=True, tab=True):
    """Convert a BED file to bigBed, preferring kent (autoSql-capable).

    ``bed_type`` is e.g. ``"bed12+14"``; required whenever *as_file* is given.
    ``extra_index`` is a list of field names to index for search (e.g. ["name"]).
    """
    bed, out = Path(bed), Path(out)
    if as_file and not bed_type:
        raise ValueError("bed_type is required when an autoSql file is supplied")
    exe = find_kent("bedToBigBed", bindir)
    if exe is None:
        if as_file:
            raise BigFileError(
                "bedToBigBed not found and autoSql fields were requested; "
                "run fetch_kent_tools() -- the pybigtools fallback cannot "
                "declare typed extra fields."
            )
        return _bed_to_bigbed_pybigtools(bed, chrom_sizes, out)
    work = bed
    if sort:
        work = out.parent / (bed.stem + ".sorted.bed")
        out.parent.mkdir(parents=True, exist_ok=True)
        sort_bed(bed, work, bindir)
    cmd = [exe]
    if tab:
        cmd.append("-tab")
    if bed_type:
        cmd.append(f"-type={bed_type}")
    if as_file:
        cmd.append(f"-as={as_file}")
    if extra_index:
        cmd.append("-extraIndex=" + ",".join(extra_index))
    cmd += [work, chrom_sizes, out]
    _run(cmd, "bedToBigBed")
    if sort and work != bed:
        work.unlink(missing_ok=True)
    return out


def bedgraph_to_bigwig(bedgraph, chrom_sizes, out, bindir=None, sort=True):
    """Convert a bedGraph to bigWig."""
    bedgraph, out = Path(bedgraph), Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    exe = find_kent("bedGraphToBigWig", bindir)
    if exe is None:
        return _bedgraph_to_bigwig_pybigtools(bedgraph, chrom_sizes, out)
    # An empty bedGraph is a legitimate state, not an error: no tool in this
    # build contributes consensus divergence (every tool is rm_fields=no or
    # divergence_only). Kent's bedGraphToBigWig aborts on it with
    # "needLargeMem: trying to allocate 0 bytes", so emit a valid 0-interval
    # bigWig instead -- contribTracks requires the file to exist regardless.
    if bedgraph.stat().st_size == 0:
        return _bedgraph_to_bigwig_pybigtools(bedgraph, chrom_sizes, out)
    work = bedgraph
    if sort:
        work = out.parent / (bedgraph.stem + ".sorted.bg")
        sort_bed(bedgraph, work, bindir)
    _run([exe, work, chrom_sizes, out], "bedGraphToBigWig")
    if sort and work != bedgraph:
        work.unlink(missing_ok=True)
    return out


def _read_sizes(chrom_sizes) -> dict[str, int]:
    sizes = {}
    with open(chrom_sizes) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split()
            sizes[f[0]] = int(f[1])
    return sizes


def _bed_to_bigbed_pybigtools(bed, chrom_sizes, out):
    import pybigtools

    sizes = _read_sizes(chrom_sizes)
    rows = []
    with open(bed) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            rows.append((f[0], int(f[1]), int(f[2]), "\t".join(f[3:])))
    rows.sort(key=lambda r: (r[0], r[1]))
    pybigtools.open(str(out), "w").write(sizes, iter(rows))
    return Path(out)


def _bedgraph_to_bigwig_pybigtools(bedgraph, chrom_sizes, out):
    import pybigtools

    sizes = _read_sizes(chrom_sizes)
    rows = []
    with open(bedgraph) as fh:
        for line in fh:
            if line.startswith(("#", "track")) or not line.strip():
                continue
            f = line.split()
            rows.append((f[0], int(f[1]), int(f[2]), float(f[3])))
    rows.sort(key=lambda r: (r[0], r[1]))
    pybigtools.open(str(out), "w").write(sizes, iter(rows))
    return Path(out)


def bigfile_info(path, bindir=None) -> dict:
    """Return summary stats for a .bb/.bw, for post-build validation."""
    path = Path(path)
    tool = "bigBedInfo" if path.suffix == ".bb" else "bigWigInfo"
    exe = find_kent(tool, bindir)
    if exe is None:
        raise BigFileError(f"{tool} not found; run fetch_kent_tools()")
    proc = subprocess.run([str(exe), str(path)], capture_output=True, text=True)
    if proc.returncode != 0:
        raise BigFileError(f"{tool} failed on {path}: {proc.stderr[:500]}")
    info = {}
    for line in proc.stdout.splitlines():
        if ":" in line:
            k, _, v = line.partition(":")
            info[k.strip()] = v.strip().replace(",", "")
    return info


def hub_check(hub_txt, bindir=None, tracks=True):
    """Run UCSC hubCheck. Returns (returncode, combined output)."""
    import subprocess
    exe = find_kent("hubCheck", bindir)
    if exe is None:
        return 0, "hubCheck not available -- skipped"
    cmd = [str(exe)] + ([] if tracks else ["-noTracks"]) + [str(hub_txt)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    return r.returncode, (r.stdout + r.stderr).strip() or "clean"


def bed_type_from_as(as_text: str, standard: int) -> str:
    """Derive the bedToBigBed `-type` argument from an autoSql definition.

    Hardcoding "bed12+10" drifts silently the moment a field is added to the
    schema -- bedToBigBed accepts an undercount without complaint, and the extra
    fields then go unindexed. Counting the autoSql is the single source of truth.
    """
    import re
    n = len(re.findall(r"^\s+\S+\s+(\w+)\s*(\[[^\]]*\])?\s*;", as_text, flags=re.M))
    if n < standard:
        raise BigFileError(f"autoSql declares {n} fields, fewer than the "
                           f"{standard} required by bed{standard}")
    return f"bed{standard}+{n - standard}"
