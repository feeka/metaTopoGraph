# metaTopoGraph

Extracts 18 de Bruijn graph topology features from a MEGAHIT SDBG at k=21.  
Runs as a standalone binary after `megahit_core` has written its SDBG to disk — no re-assembly.

## Features extracted

| # | Name | Group | Why |
|---|------|-------|-----|
| 1 | valley_position | histogram | error/signal boundary in k-mer hist — [[Sun 2018]](https://doi.org/10.1093/bioinformatics/btx637) |
| 2 | valley_depth | histogram | how well errors separate from real coverage — [[Sun 2018]](https://doi.org/10.1093/bioinformatics/btx637) |
| 3 | mult_1_fraction | histogram | fraction of singleton k-mers ≈ error rate proxy — [[Sun 2018]](https://doi.org/10.1093/bioinformatics/btx637) |
| 4 | mean_node_multiplicity | histogram | average coverage depth estimate |
| 5 | multiplicity_cv | histogram | coverage evenness — high in complex metagenomes |
| 6 | n_signal_modes | histogram | number of distinct coverage peaks ≈ strain count |
| 7 | primary_mode_depth | histogram | dominant sequencing depth |
| 8 | high_mult_tail_ratio | histogram | repeat content — spikes in repeat-rich genomes |
| 9 | tip_density | node | short dead-end paths caused by sequencing errors — [[Zerbino 2008]](https://pmc.ncbi.nlm.nih.gov/articles/PMC2336801/) |
| 10 | branching_node_fraction | node | graph complexity — encodes diversity + repeats — [[Rizzi 2019]](https://doi.org/10.1007/s40484-019-0181-x) |
| 11 | linear_node_fraction | node | fraction of unambiguous path — complement of branching |
| 12 | high_degree_node_fraction | node | dense junctions — repeats and chimeras |
| 13 | mult_at_tips | node | whether tips are low-cov errors or high-cov repeats — [[Zerbino 2008]](https://pmc.ncbi.nlm.nih.gov/articles/PMC2336801/) |
| 14 | mult_at_branches | node | whether junctions are repeat-driven or variant-driven |
| 15 | mean_tip_length | walk | longer tips → more complex errors — [[Zerbino 2008]](https://pmc.ncbi.nlm.nih.gov/articles/PMC2336801/) |
| 16 | bubble_density | walk | bubbles per branch ≈ variant / repeat density — [[Iqbal 2012]](https://pmc.ncbi.nlm.nih.gov/articles/PMC3272472/) |
| 17 | error_bubble_fraction | walk | bubbles from errors (mult ratio > 5) — [[Zerbino 2008]](https://pmc.ncbi.nlm.nih.gov/articles/PMC2336801/) |
| 18 | balanced_bubble_fraction | walk | bubbles from real variants / SNPs — [[Iqbal 2012]](https://pmc.ncbi.nlm.nih.gov/articles/PMC3272472/) |

## Build

```bash
cd metaTopoGraph/build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4
```

Requires: CMake ≥ 3.5, C++14 compiler, zlib, OpenMP.

## Usage

```bash
# From a pre-built SDBG (megahit_core must have already run)
./megahit_topo --graph <outdir>/k21 --output features.json [--sample 100000] [--threads 8]

# From raw reads (builds SDBG at k=21, extracts features, then deletes the SDBG)
./megahit_topo --reads sample.fa --output features.json [--mem 16.0] [--min-count 2] [--threads 8]
```

| Flag | Description | Default |
|---|---|---|
| `--graph` | SDBG file prefix | — |
| `--reads` | FASTA/FASTQ input (triggers SDBG build) | — |
| `--output` | Output JSON path | required |
| `--sample` | Edges sampled for node/walk/bubble phases | 100 000 |
| `--threads` | OpenMP threads | all cores |
| `--mem` | Memory for SDBG build in GB (reads mode only) | 8.0 |
| `--min-count` | Min k-mer frequency (reads mode only) | 2 |

In reads mode, `megahit_core_no_hw_accel` (or `megahit_core`) must be present in the same directory as `megahit_topo`. The temporary SDBG is written to `<binary_dir>/megahit_topo_tmp_<pid>/` and deleted on exit.

## Source layout

```
CMakeLists.txt
libs/megahit/          upstream MEGAHIT (SDBG data structures, read-only)
src/
  topo_features.h      struct definitions for all 18 features
  topo_extractor.h/.cpp  extraction logic (histogram, node, walk, bubble)
  topo_json_writer.h/.cpp  writes features.json
  main.cpp             CLI entry point
```
