#!/usr/bin/env python3
"""
pull_and_extract.py
-------------------
Download a diverse set of public metagenomes from SRA, build de Bruijn graphs
with megahit_topo (the binary produced by this repo), and extract topology
features for each one.

Requirements
  - sra-tools (fasterq-dump or fastq-dump) on PATH
    https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit
  - megahit_topo binary built from this repo:
      cd metaTopoGraph/build && cmake .. && make -j

Output layout (under --build-dir, default build/metaTopoGraph/)
  datasets/
    <name>/
      fastq_files/      (deleted after extraction unless --keep-fastq)
  outputs/
    features_<name>.json
  dataset_properties.tsv
"""

import argparse
import csv
import gzip
import json
import os
import shutil
import subprocess
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path


# ---------------------------------------------------------------------------
# Biome queries.
# Accessions are NOT hard-coded.  Each entry is a plain NCBI SRA query.
# The script resolves them to real SRR accessions at runtime via E-utilities.
# ---------------------------------------------------------------------------
BIOME_QUERIES = [
    # (name,               biome,          ncbi_sra_query)
    ("human_gut",    "human_gut",
     'gut metagenome[organism] AND WGS[strategy] AND PAIRED[layout] AND illumina[platform]'),
    ("human_oral",   "oral",
     'oral metagenome[organism] AND WGS[strategy] AND PAIRED[layout] AND illumina[platform]'),
    ("human_skin",   "skin",
     'skin metagenome[organism] AND WGS[strategy] AND PAIRED[layout] AND illumina[platform]'),
    ("marine_surf",  "marine",
     'marine metagenome[organism] AND WGS[strategy] AND PAIRED[layout] AND illumina[platform] NOT amplicon[strategy]'),
    ("marine_deep",  "marine_deep",
     'marine sediment metagenome[organism] AND WGS[strategy] AND PAIRED[layout] AND illumina[platform]'),
    ("soil_grass",   "soil",
     'soil metagenome[organism] AND grassland AND WGS[strategy] AND PAIRED[layout] AND illumina[platform]'),
    ("soil_forest",  "soil",
     'soil metagenome[organism] AND forest AND WGS[strategy] AND PAIRED[layout] AND illumina[platform]'),
    ("freshwater",   "freshwater",
     'freshwater metagenome[organism] AND lake AND WGS[strategy] AND PAIRED[layout] AND illumina[platform]'),
    ("hot_spring",   "thermophilic",
     'hot springs metagenome[organism] AND WGS[strategy] AND PAIRED[layout] AND illumina[platform]'),
    ("wastewater",   "wastewater",
     'wastewater metagenome[organism] AND WGS[strategy] AND PAIRED[layout] AND illumina[platform]'),
    ("bovine_rumen", "rumen",
     'bovine gut metagenome[organism] AND WGS[strategy] AND PAIRED[layout] AND illumina[platform]'),
    ("permafrost",   "permafrost",
     'permafrost metagenome[organism] AND WGS[strategy] AND PAIRED[layout] AND illumina[platform]'),
    ("hydro_vent",   "vent",
     'hydrothermal vent metagenome[organism] AND WGS[strategy] AND PAIRED[layout] AND illumina[platform]'),
    ("acid_mine",    "acid_mine",
     'acid mine drainage metagenome[organism] AND WGS[strategy] AND PAIRED[layout] AND illumina[platform]'),
    ("coral_reef",   "coral",
     'coral metagenome[organism] AND WGS[strategy] AND PAIRED[layout] AND illumina[platform]'),
]

DEFAULT_MAX_SPOTS = 3_000_000  # paired spots -> ~900 MB at 150 bp PE

ESEARCH  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
ENA_FILEREPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"

# Seconds to wait between NCBI API calls.
# Without an API key NCBI allows 3 req/s; with a key 10 req/s.
# We use a conservative 0.4 s (2.5 req/s) to stay safely under the limit.
_NCBI_API_KEY: str = ""       # set by --ncbi-api-key or NCBI_API_KEY env var
_NCBI_SLEEP:   float = 0.4    # overridden to 0.12 when an API key is provided


def _ncbi_get(url: str, params: dict) -> bytes:
    """GET an NCBI E-utilities URL with rate limiting and optional API key."""
    if _NCBI_API_KEY:
        params["api_key"] = _NCBI_API_KEY
    full_url = f"{url}?{urllib.parse.urlencode(params)}"
    time.sleep(_NCBI_SLEEP)
    with urllib.request.urlopen(full_url, timeout=30) as r:
        return r.read()


# ---------------------------------------------------------------------------
# SRA accession resolution (live NCBI query, no hard-coded IDs)
# ---------------------------------------------------------------------------

def resolve_accession(query: str) -> tuple:
    """
    Query NCBI SRA for `query`, return the first run that looks like a
    metagenome of reasonable size: (accession, title, total_spots).
    Returns (None, None, None) on any failure.
    """
    # Step 1: get UIDs matching the query
    try:
        uids = json.loads(_ncbi_get(ESEARCH, {
            "db":      "sra",
            "term":    query,
            "retmax":  "20",
            "retmode": "json",
        }))["esearchresult"]["idlist"]
    except Exception as exc:
        print(f"  [WARN] NCBI esearch failed: {exc}", flush=True)
        return None, None, None

    if not uids:
        print("  [WARN] No SRA results for query.", flush=True)
        return None, None, None

    # Step 2: fetch summaries and pick the first run with usable size
    try:
        sums = json.loads(_ncbi_get(ESUMMARY, {
            "db":      "sra",
            "id":      ",".join(uids),
            "retmode": "json",
        }))
    except Exception as exc:
        print(f"  [WARN] NCBI esummary failed: {exc}", flush=True)
        return None, None, None

    for uid in sums["result"].get("uids", []):
        rec = sums["result"][str(uid)]
        runs_xml = rec.get("runs", "")
        if not runs_xml:
            continue
        # Each Run element: <Run acc="SRR..." total_spots="..." total_bases="..."/>
        try:
            root = ET.fromstring(f"<runs>{runs_xml}</runs>")
        except ET.ParseError:
            continue
        for run_elem in root.findall("Run"):
            acc   = run_elem.get("acc", "")
            bases = int(run_elem.get("total_bases", 0))
            spots = int(run_elem.get("total_spots", 0))
            # Skip runs that are way too small or huge before capping
            if bases < 50_000_000:      # < 50 MB raw bases: too sparse
                continue
            title = rec.get("title", acc)
            return acc, title, spots

    print("  [WARN] No qualifying run found in the top results.", flush=True)
    return None, None, None


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run_cmd(cmd: list, check: bool = True) -> int:
    print(f"  $ {' '.join(str(c) for c in cmd)}", flush=True)
    rc = subprocess.run(cmd, check=False).returncode
    if check and rc != 0:
        raise subprocess.CalledProcessError(rc, cmd)
    return rc


def gzip_file(src: str) -> str:
    dst = src + ".gz"
    with open(src, "rb") as fi, gzip.open(dst, "wb") as fo:
        shutil.copyfileobj(fi, fo)
    os.remove(src)
    return dst


def dir_size_gb(path: str) -> float:
    total = 0
    for root, _, files in os.walk(path):
        for f in files:
            try:
                total += os.path.getsize(os.path.join(root, f))
            except OSError:
                pass
    return total / 1e9


def count_fastq_reads(path: str) -> int:
    opener = gzip.open if path.endswith(".gz") else open
    n = 0
    try:
        with opener(path, "rt", errors="replace") as fh:
            for _ in fh:
                n += 1
    except Exception:
        pass
    return n // 4


def count_reads_in_dir(directory: str) -> int:
    total = 0
    for fname in os.listdir(directory):
        if fname.endswith(".fastq") or fname.endswith(".fastq.gz"):
            total += count_fastq_reads(os.path.join(directory, fname))
    return total


def find_r1_r2(fastq_dir: str):
    """
    fasterq-dump writes <ACC>_1.fastq and <ACC>_2.fastq for paired runs.
    Returns (r1_path, r2_path) or (r1_path, None) for SE.
    """
    fqs = sorted(
        os.path.join(fastq_dir, f) for f in os.listdir(fastq_dir)
        if f.endswith(".fastq") or f.endswith(".fastq.gz")
    )
    if not fqs:
        return None, None
    r1 = next((p for p in fqs if "_1.fastq" in p), None)
    r2 = next((p for p in fqs if "_2.fastq" in p), None)
    if r1 is None:               # SE or unexpected naming: take largest as R1
        fqs_by_size = sorted(fqs, key=os.path.getsize, reverse=True)
        r1 = fqs_by_size[0]
        r2 = fqs_by_size[1] if len(fqs_by_size) > 1 else None
    return r1, r2


# ---------------------------------------------------------------------------
# ENA direct FASTQ download (fast: HTTP stream, no sra-tools join phase)
# ---------------------------------------------------------------------------

def get_ena_fastq_urls(acc: str) -> list:
    """
    Query ENA portal API for direct gzip FASTQ download URLs.
    Returns a list of HTTPS URLs (R1 first, R2 second for PE), or [] on failure.
    ENA mirrors almost all SRA runs; the FTP paths convert 1:1 to HTTPS.
    """
    params = urllib.parse.urlencode({
        "accession": acc,
        "result":    "read_run",
        "fields":    "fastq_ftp",
        "format":    "tsv",
    })
    try:
        with urllib.request.urlopen(
                f"{ENA_FILEREPORT}?{params}", timeout=30) as r:
            text = r.read().decode(errors="replace")
        lines = [l for l in text.strip().split("\n") if l]
        if len(lines) < 2:
            return []
        # ENA always returns run_accession as column 0; find fastq_ftp by header name
        headers = lines[0].split("\t")
        try:
            col_idx = headers.index("fastq_ftp")
        except ValueError:
            return []
        ftp_field = lines[1].split("\t")[col_idx]
        urls = []
        for raw in ftp_field.split(";"):
            raw = raw.strip()
            if not raw:
                continue
            # Convert ftp://... or bare path to https://
            if raw.startswith("ftp://"):
                raw = "https://" + raw[len("ftp://"):]
            elif not raw.startswith("http"):
                raw = "https://" + raw
            urls.append(raw)
        return sorted(urls)   # _1 before _2
    except Exception as exc:
        print(f"  [INFO] ENA lookup failed ({exc})", flush=True)
        return []


def stream_fastq(url: str, dest: str, max_reads: int):
    """
    Stream a gzip FASTQ from `url`, write only the first `max_reads` reads
    to `dest` (also gzip).  Closes the connection early once the limit is hit
    -- no need to download the rest of the file.
    """
    max_lines = max_reads * 4
    req = urllib.request.Request(url, headers={"Accept-Encoding": "identity"})
    CHUNK = 1 << 20  # 1 MB read buffer
    with urllib.request.urlopen(req, timeout=120) as resp:
        raw_stream = gzip.GzipFile(fileobj=resp)
        with gzip.open(dest, "wt") as out:
            buf = b""
            lines_written = 0
            done = False
            while not done:
                chunk = raw_stream.read(CHUNK)
                if not chunk:
                    break
                buf += chunk
                while b"\n" in buf:
                    line, buf = buf.split(b"\n", 1)
                    out.write(line.decode(errors="replace") + "\n")
                    lines_written += 1
                    if lines_written >= max_lines:
                        done = True
                        break


def download_via_ena(acc: str, fastq_dir: str, max_spots: int) -> bool:
    """
    Download paired (or SE) FASTQ files directly from ENA, capped at
    max_spots reads per file.  Returns True on success.
    """
    urls = get_ena_fastq_urls(acc)
    if not urls:
        return False
    os.makedirs(fastq_dir, exist_ok=True)
    for url in urls:
        fname = os.path.basename(url.split("?")[0])
        if not fname.endswith(".gz"):
            fname += ".gz"
        dest = os.path.join(fastq_dir, fname)
        print(f"  Streaming {fname} from ENA (max {max_spots:,} reads) ...",
              flush=True)
        try:
            stream_fastq(url, dest, max_spots)
        except Exception as exc:
            print(f"  [WARN] ENA stream failed for {url}: {exc}", flush=True)
            if os.path.exists(dest):
                os.remove(dest)
            return False
    return True


# Detected at startup by _probe_fasterq_flags() -- only used as fallback


def _probe_fasterq_flags() -> str:
    """
    Run fasterq-dump --help and inspect which spot-limiting flags exist.
    Returns 'NX' if -N/-X are supported, 'maxSpotCount' if that flag exists,
    or 'none' if neither is found (we will truncate after download).
    """
    try:
        result = subprocess.run(
            ["fasterq-dump", "--help"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            timeout=10
        )
        text = (result.stdout + result.stderr).decode(errors="replace")
        if "-N " in text or "--minSpotId" in text or "-N," in text:
            return "NX"
        if "--maxSpotCount" in text:
            return "maxSpotCount"
        return "none"
    except Exception:
        return "none"


def _truncate_fastq(path: str, max_reads: int):
    """Keep only the first max_reads reads in a (possibly gzipped) FASTQ."""
    max_lines = max_reads * 4
    opener_r = gzip.open if path.endswith(".gz") else open
    tmp = path + ".tmp"
    opener_w = gzip.open if path.endswith(".gz") else open
    mode_r = "rt"
    mode_w = "wt"
    with opener_r(path, mode_r) as fin, opener_w(tmp, mode_w) as fout:
        for i, line in enumerate(fin):
            if i >= max_lines:
                break
            fout.write(line)
    os.replace(tmp, path)

# ---------------------------------------------------------------------------
# Download  (ENA stream first, fasterq-dump fallback)
# ---------------------------------------------------------------------------


def download_accession(acc: str, fastq_dir: str, max_spots: int,
                        threads: int, use_fasterq: bool) -> bool:
    # ---- Try ENA direct stream first (fast, no join phase) ----
    if download_via_ena(acc, fastq_dir, max_spots):
        return True
    print("  ENA unavailable, falling back to fasterq-dump ...", flush=True)
    if not use_fasterq:
        print("  [WARN] fasterq-dump not in PATH either.", flush=True)
        return False
    # ---- fasterq-dump fallback ----
    os.makedirs(fastq_dir, exist_ok=True)
    cmd = ["fasterq-dump", acc,
           "--outdir",  fastq_dir,
           "--split-3",
           "--threads", str(threads),
           "--progress"]
    if _FASTERQ_SPOT_FLAG == "NX":
        cmd += ["-N", "1", "-X", str(max_spots)]
    elif _FASTERQ_SPOT_FLAG == "maxSpotCount":
        cmd += ["--maxSpotCount", str(max_spots)]
    rc = run_cmd(cmd, check=False)
    if rc != 0:
        print(f"  [WARN] fasterq-dump failed (exit {rc})", flush=True)
        return False
    if _FASTERQ_SPOT_FLAG == "none":
        print(f"  Truncating to {max_spots:,} spots ...", flush=True)
        for fname in os.listdir(fastq_dir):
            if fname.endswith(".fastq") or fname.endswith(".fastq.gz"):
                _truncate_fastq(os.path.join(fastq_dir, fname), max_spots)
    return True


def compress_uncompressed(fastq_dir: str):
    for fname in list(os.listdir(fastq_dir)):
        path = os.path.join(fastq_dir, fname)
        if fname.endswith(".fastq") and not fname.endswith(".fastq.gz"):
            print(f"  Compressing {fname} ...", flush=True)
            gzip_file(path)


# ---------------------------------------------------------------------------
# TSV
# ---------------------------------------------------------------------------

def write_tsv(rows: list, tsv_path: str):
    with open(tsv_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["dataset", "accession", "title", "biome",
                    "size_gb", "num_reads", "status"])
        w.writerows(rows)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build-dir", default="build/metaTopoGraph",
                    help="Root output directory (default: build/metaTopoGraph)")
    ap.add_argument("--megahit-topo", default="",
                    help="Path to megahit_topo binary (auto-detected if not given)")
    ap.add_argument("--max-spots", type=int, default=DEFAULT_MAX_SPOTS,
                    help=f"Max SRA spots per run (default {DEFAULT_MAX_SPOTS:,})")
    ap.add_argument("--threads", type=int, default=4,
                    help="Threads for fasterq-dump and megahit_topo (default 4)")
    ap.add_argument("--mem", type=float, default=8.0,
                    help="Memory GB for SDBG build (default 8.0)")
    ap.add_argument("--sample", type=int, default=200_000,
                    help="Edges sampled per dataset for feature extraction (default 200000)")
    ap.add_argument("--min-count", type=int, default=2,
                    help="Min k-mer frequency for SDBG (default 2)")
    ap.add_argument("--keep-fastq", action="store_true",
                    help="Keep FASTQ files after extraction (default: delete them)")
    ap.add_argument("--skip-download", action="store_true",
                    help="Skip SRA download; use already-present FASTQ files")
    ap.add_argument("--datasets", nargs="+", metavar="NAME",
                    help="Run only these dataset names (default: all)")
    ap.add_argument("--ncbi-api-key", default="",
                    help="NCBI API key for higher E-utilities rate limit (10 req/s). "
                         "Also read from NCBI_API_KEY env var. "
                         "Get one free at https://www.ncbi.nlm.nih.gov/account/")
    args = ap.parse_args()

    build_dir    = Path(args.build_dir).resolve()
    datasets_dir = build_dir / "datasets"
    outputs_dir  = build_dir / "outputs"
    tsv_path     = build_dir / "dataset_properties.tsv"

    datasets_dir.mkdir(parents=True, exist_ok=True)
    outputs_dir.mkdir(parents=True, exist_ok=True)

    # ---- Configure NCBI rate limiting ----
    global _NCBI_API_KEY, _NCBI_SLEEP
    _NCBI_API_KEY = (args.ncbi_api_key
                     or os.environ.get("NCBI_API_KEY", ""))
    if _NCBI_API_KEY:
        _NCBI_SLEEP = 0.12   # 10 req/s with key
        print(f"NCBI API key : set (10 req/s)")
    else:
        _NCBI_SLEEP = 0.4    # 3 req/s without key
        print("NCBI API key : not set (2.5 req/s, get one free at ncbi.nlm.nih.gov/account)")

    # ---- Locate megahit_topo binary ----
    topo = args.megahit_topo
    if not topo:
        script_dir = Path(__file__).resolve().parent
        candidates = [
            script_dir / "build" / "megahit_topo",
            script_dir / "build" / "megahit_topo.exe",
        ]
        for c in candidates:
            if c.is_file():
                topo = str(c)
                break
    if not topo:
        topo = shutil.which("megahit_topo") or ""
    if not topo:
        sys.exit(
            "ERROR: megahit_topo not found.\n"
            "Build it first:\n"
            "  cd metaTopoGraph/build && cmake .. && make -j\n"
            "Or pass --megahit-topo /path/to/megahit_topo"
        )
    print(f"megahit_topo : {topo}")
    print(f"Output root  : {build_dir}")

    # ---- Check SRA tools and probe flags ----
    global _FASTERQ_SPOT_FLAG
    use_fasterq = bool(shutil.which("fasterq-dump"))
    use_legacy  = bool(shutil.which("fastq-dump"))
    if not args.skip_download and not use_fasterq and not use_legacy:
        sys.exit(
            "ERROR: sra-tools not found in PATH.\n"
            "Install from https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit"
        )
    if use_fasterq:
        _FASTERQ_SPOT_FLAG = _probe_fasterq_flags()
        print(f"fasterq-dump : found  spot-cap flag={_FASTERQ_SPOT_FLAG}")
    else:
        print("fasterq-dump : not found, using fastq-dump (slower)")

    # ---- Filter catalogue ----
    catalogue = BIOME_QUERIES
    if args.datasets:
        wanted = set(args.datasets)
        catalogue = [d for d in BIOME_QUERIES if d[0] in wanted]
        if not catalogue:
            names = [d[0] for d in BIOME_QUERIES]
            sys.exit(f"ERROR: No match. Available names:\n  {names}")

    rows = []

    for (name, biome, query) in catalogue:
        print(f"\n{'='*62}")
        print(f"  {name}  biome={biome}")
        print(f"{'='*62}")

        fastq_dir   = datasets_dir / name / "fastq_files"
        output_json = outputs_dir / f"features_{name}.json"

        # ---- Resolve SRA accession ----
        acc, title, total_spots = None, "", 0
        if not args.skip_download:
            print(f"  Querying NCBI SRA ...", flush=True)
            acc, title, total_spots = resolve_accession(query)
            if acc is None:
                print(f"  [SKIP] Could not resolve accession for {name}")
                rows.append([name, "n/a", "", biome, 0, 0, "resolve_failed"])
                continue
            print(f"  Resolved: {acc}  ({title})", flush=True)
            print(f"  Downloading (max {args.max_spots:,} spots) ...", flush=True)
            ok = download_accession(acc, str(fastq_dir), args.max_spots,
                                    args.threads, use_fasterq)
            if not ok:
                rows.append([name, acc, title, biome, 0, 0, "download_failed"])
                continue
        else:
            if not fastq_dir.exists():
                print(f"  [SKIP] {fastq_dir} does not exist")
                continue
            acc, title = name, "(existing)"

        # ---- Locate FASTQ files ----
        r1, r2 = find_r1_r2(str(fastq_dir))
        if r1 is None:
            print(f"  [WARN] No FASTQ files in {fastq_dir}")
            rows.append([name, acc, title, biome, 0, 0, "no_fastq"])
            continue

        size_gb = dir_size_gb(str(fastq_dir))
        n_reads = count_reads_in_dir(str(fastq_dir))
        print(f"  FASTQ: {size_gb:.2f} GB  {n_reads:,} reads", flush=True)

        # ---- Extract features ----
        print(f"  Extracting -> {output_json.name} ...", flush=True)
        cmd = [
            topo,
            "--reads",     r1,
            "--output",    str(output_json),
            "--sample",    str(args.sample),
            "--threads",   str(args.threads),
            "--mem",       str(args.mem),
            "--min-count", str(args.min_count),
        ]
        if r2:
            cmd += ["--reads2", r2]

        t0 = time.time()
        rc = run_cmd(cmd, check=False)
        elapsed = time.time() - t0

        if rc != 0:
            print(f"  [WARN] megahit_topo failed (exit {rc})", flush=True)
            rows.append([name, acc, title, biome,
                         f"{size_gb:.2f}", n_reads, "extraction_failed"])
            continue

        print(f"  OK ({elapsed:.1f}s)", flush=True)

        # ---- Delete FASTQ files (save disk) ----
        if not args.keep_fastq:
            print(f"  Deleting FASTQ files to save disk ...", flush=True)
            shutil.rmtree(str(fastq_dir), ignore_errors=True)

        rows.append([name, acc, title, biome,
                     f"{size_gb:.2f}", n_reads, "ok"])

    # ---- Write TSV ----
    write_tsv(rows, str(tsv_path))
    print(f"\nDataset summary -> {tsv_path}")

    ok_count   = sum(1 for r in rows if r[-1] == "ok")
    fail_count = len(rows) - ok_count
    print(f"Done: {ok_count}/{len(rows)} OK" +
          (f"  ({fail_count} failed)" if fail_count else ""))


if __name__ == "__main__":
    main()
