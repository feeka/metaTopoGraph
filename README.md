# metaTopoGraph

Extracts 18 de Bruijn graph topology features from a MEGAHIT SDBG at k=21.  
Runs as a standalone binary after `megahit_core` has written its SDBG to disk — no re-assembly.

## Features extracted

| # | Name | Group |
|---|------|-------|
| 1 | valley_position | histogram |
| 2 | valley_depth | histogram |
| 3 | mult_1_fraction | histogram |
| 4 | mean_node_multiplicity | histogram |
| 5 | multiplicity_cv | histogram |
| 6 | n_signal_modes | histogram |
| 7 | primary_mode_depth | histogram |
| 8 | high_mult_tail_ratio | histogram |
| 9 | tip_density | node |
| 10 | branching_node_fraction | node |
| 11 | linear_node_fraction | node |
| 12 | high_degree_node_fraction | node |
| 13 | mult_at_tips | node |
| 14 | mult_at_branches | node |
| 15 | mean_tip_length | walk |
| 16 | bubble_density | walk |
| 17 | error_bubble_fraction | walk |
| 18 | balanced_bubble_fraction | walk |

## Build

```bash
cd metaTopoGraph/build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4
```

Requires: CMake ≥ 3.5, C++14 compiler, zlib, OpenMP.

## Usage

```bash
# 1. Run megahit_core as usual to build the SDBG (produces <outdir>/k21.sdbg_info etc.)
# 2. Extract topology features
./megahit_topo --graph <outdir>/k21 --output features.json [--sample 100000] [--threads 8]
```

`--graph` is the SDBG file prefix (same string passed to MEGAHIT's sdbg build step).  
`--sample` controls how many edges are visited during node/walk/bubble phases (default 100 000).  
`--threads` sets the number of OpenMP threads for the histogram phase (default: all available cores).

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
