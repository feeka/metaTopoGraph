# MetaTopoGraph — Integration Plan for MEGAHIT

This plan is for execution by another Claude instance in VS Code. The goal is to extend MEGAHIT to extract de Bruijn graph topology features at k=21 from the existing SDBG data structures.

---

## Goal

Add a feature extraction stage to MEGAHIT that runs after the SDBG is built at k=21. The stage extracts 18 graph topology features from the existing SDBG data structures.

## Scope and Non-Goals

**In scope:**
- Add a new compilation unit for topology feature extraction
- Hook into the existing SDBG construction pipeline at k=21
- Reuse MEGAHIT's existing node/edge traversal infrastructure

**Not in scope:**
- Changing MEGAHIT's default assembly behavior in any way
- Adding ML, classification, or prediction
- Modifying SDBG construction itself
- Producing contigs or running the iteration stage

---

## Where This Lives in the Source Tree

The MEGAHIT source layout is:

```
src/
├── assembly/        SDBG-based assembly logic (sdbg_pruning.cpp etc.)
├── sdbg/            Succinct de Bruijn graph data structure
├── kmlib/           Internal k-mer utilities (bit vectors etc.)
├── sequence/        Sequence handling
├── sorting/         External sorting
├── localasm/        Local assembly
├── iterate/         k-mer iteration
├── utils/           Histogram, helpers
├── tools/           Auxiliary binaries
├── main_assemble.cpp
├── main_sdbg_build.cpp
└── megahit          Python driver script
```

In the project directory(metaTopoGraph):

```
...
libs/megahit
src/
├── topo_extractor.h
├── topo_extractor.cpp
├── topo_features.h        # struct definitions
├── topo_json_writer.h
├── topo_json_writer.cpp
├── main.cpp
CMakeLists.txt
```

---

## Workflow

For topology extraction the flow becomes:

1. Build library (unchanged)
2. Build SDBG at k=21 (unchanged)
3. **Run topology extraction (new)** → produces features.json → exit
4. (Skip everything else)

The new step reads the SDBG that step 2 already wrote to disk. No duplicate work.

---

## The 18 Features

These are the features the C++ implementation extracts from the SDBG.

### From the SDBG metadata and edge multiplicity histogram

The SDBG class already exposes `EdgeMultiplicity(edge_id)` and the histogram is computed during pruning (`sdbg_pruning.cpp` builds a `Histgram<mul_t>`).

| # | Name | Source |
|---|---|---|
| 1 | valley_position | argmin of smoothed histogram between mult 2 and 50 |
| 2 | valley_depth | count_at_valley / count_at_error_peak |
| 3 | mult_1_fraction | hist[1] / total_edges |
| 4 | mean_node_multiplicity | weighted mean of multiplicity histogram |
| 5 | multiplicity_cv | std / mean of multiplicity |
| 6 | n_signal_modes | peak count in smoothed histogram right of valley |
| 7 | primary_mode_depth | argmax right of valley |
| 8 | high_mult_tail_ratio | sum(hist[3*primary:]) / total |

### From SDBG node-level traversal

The SDBG class supports `Outdegree(node)`, `Indegree(node)`, and edge iteration. Sample n nodes uniformly using SDBG's existing iteration machinery.

| # | Name | Source |
|---|---|---|
| 9 | tip_density | nodes where one-side degree is 0 |
| 10 | branching_node_fraction | nodes with in-deg ≥ 2 OR out-deg ≥ 2 |
| 11 | linear_node_fraction | nodes with in-deg=1 AND out-deg=1 |
| 12 | high_degree_node_fraction | nodes with in-deg + out-deg ≥ 4 |
| 13 | mult_at_tips | mean multiplicity of tip nodes |
| 14 | mult_at_branches | mean multiplicity of branching nodes |

### From bounded local walks

MEGAHIT's `assembly/sdbg_pruning.cpp` already contains tip-walking logic in the `Trim` function. Reuse the same traversal pattern.

| # | Name | Source |
|---|---|---|
| 15 | mean_tip_length | walk from each sampled tip until branch or max length |
| 16 | bubble_density | bubbles found per sampled branching node |
| 17 | error_bubble_fraction | bubbles where multiplicity ratio > 5 |
| 18 | balanced_bubble_fraction | bubbles where ratio < 2 AND both branches above valley_position |

---

## Implementation Order

Execute in this order. Each step compiles and is testable before moving on.

### Step 1: Create the build infrastructure

- Create `src/` directory and `CMakeLists.txt`
- Update top-level `CMakeLists.txt` to add `src/` as a subdirectory and to define a new executable `megahit_topo` from `main.cpp`
- Verify the project still compiles with `cmake .. && make -j4`

### Step 2: Define the feature data structures

File `src/topo_features.h`:

- Define `struct HistogramFeatures` with the 8 histogram-derived fields
- Define `struct NodeFeatures` with the 6 node-level fields
- Define `struct WalkFeatures` with the 4 walk-derived fields
- Define `struct TopoFeatures` containing all three plus a timing struct
- Use `double` for all floating point fields
- Use `uint64_t` for counts

### Step 3: JSON writer

File `src/topo_json_writer.cpp`:

- Function `void WriteFeatures(const TopoFeatures& features, const std::string& output_path)`
- Hand-write JSON output to avoid adding a JSON library dependency. The schema is small and fixed.
- Test by populating a `TopoFeatures` struct manually and verifying the file content

### Step 4: Histogram features

File `src/topo_extractor.cpp`, function `ExtractHistogramFeatures`:

- Take an `SDBG&` reference as input
- Build the multiplicity histogram by iterating valid edges and calling `dbg.EdgeMultiplicity(i)`
- This logic mirrors `sdbg_pruning::InferMinDepth` — reuse its histogram-building style
- Implement features 1–8
- For valley finding: smooth the histogram with a simple gaussian kernel (write a small inline gaussian smoother, no external dependency), then find the local minimum between mult 2 and 50
- Test against a few hand-constructed histograms

### Step 5: Node-level features

File `src/topo_extractor.cpp`, function `ExtractNodeFeatures`:

- Take `SDBG&` and a sample size (default 100000)
- Use `dbg.size()` and a stride-based sampling strategy: visit every `dbg.size() / sample_size`-th node
- For each sampled node, query `dbg.Outdegree(node)` and `dbg.Indegree(node)`
- For each sampled node, query `dbg.NodeMultiplicity(node)` (or compute from edge multiplicities — confirm which method MEGAHIT exposes)
- Classify and aggregate features 9–14
- Skip nodes where `dbg.IsValidNode(node)` returns false

### Step 6: Tip walking

File `src/topo_extractor.cpp`, function `WalkTipFromNode`:

- Take a tip node, return its walked length in base pairs
- Walk in the non-dead direction using the same pattern as `sdbg_pruning::Trim` — call `dbg.PrevSimplePathEdge` or `dbg.NextSimplePathEdge` (verify exact method names in current MEGAHIT)
- Stop when reaching a branch (degree ≥ 2 on either side), another dead end, or `max_length = 5000` bp
- Function `ExtractTipWalkFeatures` averages over sampled tips → feature 15

### Step 7: Bubble detection

File `src/topo_extractor.cpp`, function `FindBubblesFromBranch`:

- For a node with out-degree 2, follow both outgoing branches with bounded BFS
- Maximum bubble length: 500 bp
- Maximum walk depth: 5 edges per branch
- A bubble is detected if both walks reach the same node within those bounds
- For each detected bubble, record the multiplicity of both branches
- Classify the bubble:
  - Error bubble: max(b1_mult, b2_mult) / min(b1_mult, b2_mult) > 5
  - Balanced bubble: ratio < 2 AND both branches > valley_position
- Aggregate over sampled branching nodes → features 16–18

### Step 8: Top-level extractor

File `src/topo_extractor.cpp`, function `RunExtraction`:

- Signature: `TopoFeatures RunExtraction(SDBG& dbg, const ExtractionOptions& opts)`
- Calls all four extraction functions in order
- Records timing for each phase using `std::chrono::steady_clock`
- Returns a populated `TopoFeatures`
