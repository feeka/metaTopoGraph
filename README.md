<p align="center">
  <img src="icon.png" alt="metaTopoGraph" width="350"/>
</p>

Extracts de Bruijn graph topology features from metagenomic reads.
Builds a MEGAHIT SDBG internally from raw reads, extracts features in a single pass, then cleans up.

## Features

```
n_reads               total reads provided (R1 + R2)
n_kmers               total valid k-mer edges in the graph (E_real + E_err)
mult_min              minimum edge multiplicity
mult_max              maximum edge multiplicity
mult_mean             mean edge multiplicity
mean_min_branch_ratio mean of min(mul_i / mul_out_j) at branching nodes (out-degree >= 2)
mean_max_branch_ratio mean of max(mul_i / mul_out_j) at branching nodes; approximates real genome coverage
max_branch_ratio      global maximum of (mul_i / mul_out_j) across all branches
min_branch_ratio      global minimum of (mul_i / mul_out_j) across all branches
n_tips                edges with in-degree == 0 or out-degree == 0
n_branch_nodes        edges with out-degree >= 2
```

## Build

Clone with submodules — MEGAHIT lives under `libs/megahit/` as a git submodule:

```bash
git clone --recurse-submodules https://github.com/feeka/metaTopoGraph.git
# or, if already cloned:
git submodule update --init --recursive
```

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
# single-end
./megahit_topo --reads sample.fa --output features.json

# paired-end
./megahit_topo --reads R1.fq.gz --reads2 R2.fq.gz --output features.json

# pre-built SDBG
./megahit_topo --graph sdbg_prefix --output features.json --read-count 2000000
```

```
--reads        R1 (or single-end) FASTA/FASTQ reads file
--reads2       R2 FASTA/FASTQ file for paired-end input
--output       Output JSON path (required)
--threads      OpenMP threads (default: all cores)
--mem          Memory for SDBG build in GB (default: 90% of RAM)
--min-count    Min k-mer frequency (default: 1)
--kmer-size    k-mer length (default: 21, must be odd)
--keep-graph   Keep the temporary SDBG directory after extraction
--read-count   Total read count for --graph mode (auto-detected in --reads mode)
```

## Output format

```json
{
  "node": {
    "n_reads": 2000000,
    "n_kmers": 20971520,
    "mult_min": 1,
    "mult_max": 1231,
    "mult_mean": 12.4418,
    "mean_min_branch_ratio": 0.982312,
    "mean_max_branch_ratio": 28.9574,
    "max_branch_ratio": 1231.0,
    "min_branch_ratio": 1.0,
    "n_tips": 409078,
    "n_branch_nodes": 512345
  },
  "timing": {
    "node_ms": 5136.02
  }
}
```

## Source layout

```
CMakeLists.txt
libs/megahit/              upstream MEGAHIT (SDBG data structures, read-only)
src/
  topo_features.h          NodeFeatures and TopoFeatures structs
  topo_extractor.h/.cpp    single-pass full-scan extraction
  topo_json_writer.h/.cpp  JSON output
  main.cpp                 CLI entry point
  README.md                implementation reference
  error_threshold_derivation.md  derivation of the Bayesian error threshold
```
