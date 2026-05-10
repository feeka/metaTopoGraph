<p align="left">
  <img src="icon.png" alt="metaTopoGraph" width="250"/>
</p>

# metaTopoGraph

Builds a de Bruijn graph from metagenomic reads and extracts graph topology features. If you provide a reference FASTA, it also runs jellyfish to label each k-mer as a sequencing error or true genomic sequence, and finds the best multiplicity thresholds to separate them.

All temp files are cleaned up automatically.

## Requirements

- CMake ≥ 3.5, C++14, zlib, OpenMP
- [jellyfish 2](https://github.com/gmarcais/Jellyfish) on `PATH` — only needed with `--ref-fasta`

## Build

```bash
git clone --recurse-submodules https://github.com/feeka/metaTopoGraph.git
bash build.sh
```

Binaries land in `~/metaTopoGraph-build/`. `build.sh` handles MEGAHIT and `megahit_topo` in one shot.

## Usage

```bash
# graph features only
megahit_topo --reads R1.fq.gz --reads2 R2.fq.gz --output features.json

# graph features + k-mer error features + labels
megahit_topo --reads R1.fq.gz --reads2 R2.fq.gz \
             --ref-fasta reference.fa --output features.json
```

| Flag | Description |
|---|---|
| `--reads` | R1 (or SE) reads file — FASTA/FASTQ, plain or gzipped |
| `--reads2` | R2 reads file (paired-end) |
| `--ref-fasta` | Reference FASTA — enables k-mer features and labels |
| `--output` | Output JSON path |
| `--threads` | Thread count (default: all cores) |
| `--mem` | Memory for SDBG build in GB (default: 90% of RAM) |
| `--min-count` | Min k-mer frequency for SDBG build (default: 1) |
| `--kmer-size` | k-mer length (default: 21, must be odd) |
| `--keep-graph` | Don't delete the temp SDBG directory |
| `--graph` | Use a pre-built SDBG instead of reads |
| `--read-count` | Total reads — required with `--graph`, auto-detected otherwise |

## Output

```json
{
  "node": {
    "n_kmers": 20971520,
    "mult_mean": 12.44,
    "mult_max_ratio": 98.9,
    "mean_max_branch_ratio": 28.96,
    "max_branch_ratio": 1231.0,
    "mean_jump_ratio": 6.12,
    "tip_density": 0.0195,
    "branch_density": 0.0244
  },
  "kmer_features": {
    "err_hist_node":          [0.71, 0.14, 0.08, 0.04, 0.02, 0.01, 0.00],
    "err_hist_read":          [0.52, 0.20, 0.14, 0.08, 0.04, 0.02, 0.00],
    "overlap_density":        [0.94, 0.41, 0.10, 0.02, 0.01, 0.005, 0.001],
    "err_neighbor_mult":      [3.8, 7.3, 12.1, 14.9, 27.5, 39.2, 92.1],
    "mean_error_mult":        1.32,
    "mean_true_mult":         18.7,
    "neighbor_contrast_ratio": 14.4,
    "isolated_error_frac":    0.048
  },
  "labels": {
    "min_count": 3,
    "contrast_threshold": 5.2,
    "error_fraction": 0.044
  },
  "timing": { "node_ms": 5136.02, "kmer_ms": 12340.0 }
}
```

`kmer_features`, `labels`, and `kmer_ms` are omitted when `--ref-fasta` is not given.

## Features explained

### Graph topology (`node`)

| Field | What it is |
|---|---|
| `n_kmers` | Total k-mer edges in the SDBG (both strands counted) |
| `mult_mean` | Mean edge multiplicity ≈ sequencing coverage |
| `mult_max_ratio` | `mult_max / mult_mean` — how far the peak is above average |
| `mean_max_branch_ratio` | At each branching node: ratio of the incoming edge mult to the highest outgoing edge mult, averaged — approximates true genome coverage |
| `max_branch_ratio` | Same ratio, global maximum |
| `mean_jump_ratio` | Per-branch-node: `max(hi,lo) / min(hi,lo)` between incoming and outgoing — symmetric version of the above |
| `tip_density` | Dead-end edges / total edges |
| `branch_density` | Branching nodes / total edges |

### K-mer error features (`kmer_features`)

Multiplicity bins: `[1], [2–3], [4–7], [8–15], [16–31], [32–63], [64+]`.

| Field | What it is |
|---|---|
| `err_hist_node[7]` | Fraction of error k-mers in each bin (by unique count) |
| `err_hist_read[7]` | Fraction of error k-mer occurrences in each bin (weighted by multiplicity) |
| `overlap_density[7]` | P(error) per bin — how contaminated each multiplicity range is |
| `err_neighbor_mult[7]` | Mean max-neighbour multiplicity for error k-mers in each bin |
| `mean_error_mult` | Mean multiplicity of error k-mers |
| `mean_true_mult` | Mean multiplicity of true k-mers |
| `neighbor_contrast_ratio` | `mean(max_nbr_mult of errors) / mean_error_mult` |
| `isolated_error_frac` | Fraction of error k-mers with no neighbour present in reads |

### Labels

Two thresholds are found by scanning all possible values and picking the one with the best F1 against the jellyfish ground truth:

- **`min_count`** — flag a k-mer as error if `mult < min_count`
- **`contrast_threshold`** — AND condition on top: also require `max_neighbour_mult / mult > contrast_threshold`
- **`error_fraction`** — fraction of unique canonical k-mers not found in the reference

> `error_fraction` counts unique k-mers, not read occurrences. One base error corrupts up to k overlapping k-mers, so this number is much higher than the raw base error rate. Coverage-weighted rate ≈ `(error_fraction × mean_error_mult) / mult_mean`.

## Source

```
build.sh
CMakeLists.txt
libs/megahit/           MEGAHIT (submodule)
src/
  topo_features.h       structs: NodeFeatures, KmerFeatures, Labels, TopoFeatures
  topo_extractor.h/.cpp graph topology pass + jellyfish k-mer extraction
  topo_json_writer.h/.cpp  JSON output
  main.cpp              CLI
  README.md             implementation notes
```


## Features

### Graph topology (always extracted)

- **`n_kmers`** — total valid k-mer edges in the SDBG (both strands)
- **`mult_mean`** — mean edge multiplicity ≈ sequencing coverage
- **`mult_max_ratio`** — mult_max / mult_mean — normalised coverage headroom
- **`mean_max_branch_ratio`** — mean of max(mult_i / mult_out_j) at branching nodes ≈ C_real
- **`max_branch_ratio`** — global maximum of (mult_i / mult_out_j) across all branches
- **`mean_jump_ratio`** — mean per-branch-node of max(hi,lo)/min(hi,lo) — symmetric ratio
- **`tip_density`** — n_tips / n_kmers
- **`branch_density`** — n_branch_nodes / n_kmers

### K-mer error features (extracted when `--ref-fasta` is supplied)

Multiplicity bins: `[1], [2–3], [4–7], [8–15], [16–31], [32–63], [64+]` (N_BINS = 7).

- **`err_hist_node[7]`** — fraction of error k-mers (by count) in each bin
- **`err_hist_read[7]`** — fraction of error k-mer occurrences (sum of mult) in each bin
- **`overlap_density[7]`** — P(error | mult in bin) per bin
- **`err_neighbor_mult[7]`** — mean max-neighbour multiplicity for error k-mers per bin
- **`mean_error_mult`** — mean multiplicity of sequencing-error k-mers
- **`mean_true_mult`** — mean multiplicity of true (reference) k-mers
- **`neighbor_contrast_ratio`** — mean(max_nbr_mult of errors) / mean_error_mult
- **`isolated_error_frac`** — fraction of error k-mers with no neighbour present in reads

### Labels (computed jointly with k-mer features)

- **`min_count`** — optimal multiplicity threshold: predict_error = (mult < min_count)
- **`contrast_threshold`** — optimal neighbour contrast threshold (AND condition): predict_error = (mult < min_count) AND (max_nbr/mult > contrast_threshold)
- **`error_fraction`** — fraction of unique canonical read k-mers absent from the reference (unique k-mer fraction, not coverage-weighted; see note below)

Both thresholds are chosen independently to maximise F1 against the jellyfish ground truth.

> **Note on `error_fraction`**: this counts unique canonical k-mers, not read occurrences.
> A single sequencing error corrupts up to k overlapping k-mers, each appearing only once,
> so this value is significantly higher than the raw base error rate.
> The coverage-weighted error rate is approximately
> `(error_fraction × mean_error_mult) / mult_mean`.

## Dependencies

- CMake ≥ 3.5, C++14 compiler, zlib, OpenMP
- **jellyfish 2** must be on `PATH` when `--ref-fasta` is used

## Build

Clone with submodules — MEGAHIT lives under `libs/megahit/` as a git submodule:

```bash
git clone --recurse-submodules https://github.com/feeka/metaTopoGraph.git
# or, if already cloned:
git submodule update --init --recursive
```

**Quick build** (recommended):

```bash
bash build.sh
```

The script builds MEGAHIT, then `megahit_topo`, placing both binaries under `~/metaTopoGraph-build/`.

**Manual build**:

```bash
# 1. build MEGAHIT
cmake -S libs/megahit -B ~/megahit-build -DCMAKE_BUILD_TYPE=Release
make -C ~/megahit-build -j$(nproc)

# 2. build megahit_topo
cmake -S . -B ~/metaTopoGraph-build -DCMAKE_BUILD_TYPE=Release
make -C ~/metaTopoGraph-build -j$(nproc)
```

## Usage

```bash
# topology features only — single-end reads
./megahit_topo --reads sample.fa --output features.json

# topology features only — paired-end reads
./megahit_topo --reads R1.fq.gz --reads2 R2.fq.gz --output features.json

# topology + k-mer features + labels — paired-end with reference
./megahit_topo --reads R1.fq.gz --reads2 R2.fq.gz \
               --ref-fasta reference.fa --output features.json

# pre-built SDBG (topology only)
./megahit_topo --graph sdbg_prefix --output features.json --read-count 2000000
```

```
--reads        R1 (or single-end) FASTA/FASTQ reads file (required unless --graph)
--reads2       R2 FASTA/FASTQ file for paired-end input
--ref-fasta    Reference FASTA for jellyfish error-k-mer classification
               Enables k-mer features and labels in the output.
--output       Output JSON path (required)
--threads      OpenMP threads (default: all cores)
--mem          Memory for SDBG build in GB (default: 90% of RAM)
--min-count    Min k-mer frequency for SDBG build (default: 1)
--kmer-size    k-mer length (default: 21, must be odd)
--keep-graph   Keep the temporary SDBG directory after extraction
--read-count   Total read count for --graph mode (auto-detected in --reads mode)
```

## Output format

```json
{
  "node": {
    "n_kmers": 20971520,
    "mult_mean": 12.44,
    "mult_max_ratio": 98.9,
    "mean_max_branch_ratio": 28.96,
    "max_branch_ratio": 1231.0,
    "mean_jump_ratio": 6.12,
    "tip_density": 0.0195,
    "branch_density": 0.0244
  },
  "kmer_features": {
    "err_hist_node":        [0.71, 0.14, 0.08, 0.04, 0.02, 0.01, 0.00],
    "err_hist_read":        [0.52, 0.20, 0.14, 0.08, 0.04, 0.02, 0.00],
    "overlap_density":      [0.94, 0.41, 0.10, 0.02, 0.01, 0.005, 0.001],
    "err_neighbor_mult":    [3.8, 7.3, 12.1, 14.9, 27.5, 39.2, 92.1],
    "mean_error_mult":      1.32,
    "mean_true_mult":       18.7,
    "neighbor_contrast_ratio": 14.4,
    "isolated_error_frac":  0.048
  },
  "labels": {
    "min_count": 3,
    "contrast_threshold": 5.2,
    "error_fraction": 0.044
  },
  "timing": {
    "node_ms": 5136.02,
    "kmer_ms": 12340.0
  }
}
```

`kmer_features`, `labels`, and `kmer_ms` are absent when `--ref-fasta` is not provided.

## Temporary files

In `--reads` mode the tool creates a single temp directory (`megahit_topo_tmp_<pid>/`) next to the binary.
All SDBG files and jellyfish intermediates (`ref.jf`, `reads.jf`, `reads_dump.fa`, `query.txt`) live there and are removed automatically on exit unless `--keep-graph` is passed.

## Source layout

```
CMakeLists.txt
libs/megahit/              upstream MEGAHIT (SDBG data structures, read-only)
src/
  topo_features.h          NodeFeatures, KmerFeatures, Labels, TopoFeatures structs
  topo_extractor.h/.cpp    graph topology extraction + jellyfish k-mer extraction
  topo_json_writer.h/.cpp  JSON output
  main.cpp                 CLI entry point
  README.md                implementation reference
```
