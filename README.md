<p align="center">
  <img src="icon.png" alt="metaTopoGraph" width="350"/>
</p>

Extracts 7 de Bruijn graph topology features from a MEGAHIT SDBG at k=21 by default.
Runs as a standalone binary from raw reads — builds the SDBG internally, extracts features, then cleans up.

## Features extracted

A single full scan over all edges in the SDBG produces:

```
n_prominent_jumps     — number of directed edge transitions where the multiplicity ratio
                        between adjacent edges is >= 3.0 (JUMP_THRESHOLD). Encodes error
                        branch density and coverage unevenness.

mult_min              — minimum edge multiplicity across the entire graph.
mult_max              — maximum edge multiplicity. High values indicate repeated k-mers
                        (rRNA operons, transposons, conserved genes).
mult_mean             — mean edge multiplicity. Depressed below true coverage by the
                        presence of error k-mers at multiplicity ~1.

mean_min_branch_ratio — at branching nodes (out-degree >= 2): mean of the minimum
                        mult[source] / mult[outgoing_j] ratio across all outgoing edges.
                        Near 1.0 means the real path is always present at branches.

mean_max_branch_ratio — at branching nodes (out-degree >= 2): mean of the maximum
                        mult[source] / mult[outgoing_j] ratio. Approximates the real
                        genome coverage C, because the dominant real path consistently
                        wins against the error branches (mult ~1).

n_tips                — edges with in-degree == 0 or out-degree == 0. Caused by
                        sequencing errors near read boundaries or genuine contig ends.
```

`mean_max_branch_ratio` and `mult_mean` together allow a Bayesian estimate of the
per-base sequencing error rate and the multiplicity threshold below which an edge is
more likely erroneous than real — with no external calibration needed.

## Build

Clone with submodules — MEGAHIT lives under `libs/megahit/` as a git submodule:

```bash
git clone --recurse-submodules https://github.com/feeka/metaTopoGraph.git
# or, if already cloned:
git submodule update --init --recursive
```

Use the provided build script (WSL / Linux):

```bash
bash build.sh
```

This builds MEGAHIT first, then `megahit_topo`, and symlinks the required
`megahit_core_no_hw_accel` binary into the topo build directory.

To build manually:

```bash
# 1. build MEGAHIT
cmake -S libs/megahit -B ~/megahit-build -DCMAKE_BUILD_TYPE=Release
make -C ~/megahit-build -j$(nproc)

# 2. build megahit_topo
cmake -S . -B ~/metaTopoGraph-build -DCMAKE_BUILD_TYPE=Release
make -C ~/metaTopoGraph-build -j$(nproc)
```

Requires: CMake >= 3.5, C++14 compiler, zlib, OpenMP.

## Usage

```bash
./megahit_topo --reads sample.fa --output features.json [--mem 16.0] [--min-count 2] [--threads 8]
```

```
--reads       FASTA/FASTQ input (required)
--output      Output JSON path (required, must end in .json)
--threads     OpenMP threads (default: all cores)
--mem         Memory for SDBG build in GB (default: 8.0)
--min-count   Min k-mer frequency; 2 filters most error k-mers at coverage >= 5x (default: 2)
```

`megahit_core_no_hw_accel` (or `megahit_core`) must be present in the same directory as
`megahit_topo`. The temporary SDBG is written to `<binary_dir>/megahit_topo_tmp_<pid>/`
and deleted on exit.

## Output format

```json
{
  "node": {
    "n_prominent_jumps": 1591614,
    "mult_min": 1,
    "mult_max": 1231,
    "mult_mean": 12.4418,
    "mean_min_branch_ratio": 0.982312,
    "mean_max_branch_ratio": 28.9574,
    "n_tips": 409078
  },
  "timing": {
    "node_ms": 5136.02
  }
}
```

## Python scripts

Three companion scripts cover the full dataset collection → analysis → visualisation workflow.
They require Python >= 3.9 and `pandas`, `numpy`, `scipy`, `matplotlib`, `seaborn`.

### `pull_and_extract.py` — bulk dataset collection

Downloads public metagenomes from NCBI SRA, builds the SDBG with `megahit_topo`, and
extracts features for each one. Supports 55 biome categories; each FASTQ is deleted
immediately after extraction to save disk.

```bash
# pull 1 dataset per category (55 total, default)
python pull_and_extract.py --build-dir /data/metagenomes --threads 8

# pull 2 per category, randomise spot count between 1.5M–3M
python pull_and_extract.py --n-datasets 110 --max-spots 3000000 --threads 8

# keep raw reads and the SDBG for inspection
python pull_and_extract.py --keep-fastq --keep-graph
```

```
--build-dir       Root output directory (default: build/metaTopoGraph)
--n-datasets N    Total datasets, spread evenly across categories (default: 1 per category)
--max-spots X     Upper spot cap; each dataset draws randomly from [X/2, X] (default: 3000000)
--kmer-size       k-mer length passed to SDBG build (default: 21)
--keep-fastq      Don't delete FASTQ after extraction
--keep-graph      Don't delete temporary SDBG
--threads         Threads for download and SDBG build (default: 4)
--ncbi-api-key    NCBI API key for higher rate limit (10 req/s)
```

### `analyze.py` — descriptive statistics and correlations

Reads every `features_*.json` in a folder and produces CSVs and plots:
summary statistics, Pearson/Spearman correlation matrices, and a hierarchical clustering heatmap.

```bash
python analyze.py --input build/metaTopoGraph/outputs --output results/
```

### `plot_features.py` — per-sample and multi-sample visualisation

```bash
# single sample
python plot_features.py features.json

# whole directory
python plot_features.py build/metaTopoGraph/outputs/
```

## Source layout

```
CMakeLists.txt
build.sh                 one-shot WSL/Linux build script
libs/megahit/            upstream MEGAHIT (SDBG data structures, read-only)
src/
  topo_features.h        NodeFeatures and TopoFeatures structs
  topo_extractor.h/.cpp  single-pass full-scan extraction
  topo_json_writer.h/.cpp  writes features.json
  main.cpp               CLI entry point
```
