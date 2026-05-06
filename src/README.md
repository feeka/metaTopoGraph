# Feature Extraction — Implementation Reference

All 18 features are computed from MEGAHIT's Succinct de Bruijn Graph (SDBG) by
`topo_extractor.cpp` in four sequential phases. Each phase is independently
OpenMP-parallelised. Extraction options are in `ExtractionOptions` (defaults
shown); features and timing are returned in `TopoFeatures`.

---

## Phase 1 — Multiplicity Histogram (features 1–8)

**Source:** `ExtractHistogramFeatures(SDBG&)`

**Method:** Full edge scan. Every valid edge in the SDBG is visited once. The
scan is parallelised with per-thread private histograms (size 1 001, cap at
multiplicity 1 000) that are merged after the parallel section to avoid false
sharing. A discrete Gaussian kernel (σ = 5, radius = 15) is then convolved over
the raw counts to smooth noise before valley/mode detection.

| # | Field | Formula / Method |
|---|-------|-----------------|
| 1 | `valley_position` | argmin of the smoothed histogram in the range [2, max(50, rough_peak/2)] where *rough_peak* is the argmax of the smoothed histogram for multiplicity ≥ 10. The floor of 50 prevents the search from stopping too early at moderate coverage; the half-peak heuristic is an empirical guess with no theoretical basis. |
| 2 | `valley_depth` | `smoothed[valley_pos] / smoothed[1]`. Measures how clearly the error–signal dip separates the two populations. Values near 0 indicate a clean separation; values near 1 indicate no visible valley. |
| 3 | `mult_1_fraction` | `hist[1] / total_edges`. Raw fraction; computed directly from the unsmoothed histogram before Gaussian smoothing. |
| 4 | `mean_node_multiplicity` | `Σ(i · hist[i]) / total` for i = 1..1000. Weighted mean of the raw histogram. |
| 5 | `multiplicity_cv` | `sqrt(Σ(i² · hist[i])/total − mean²) / mean`. Coefficient of variation of the multiplicity distribution. Uses the same raw histogram pass as feature 4. |
| 6 | `n_signal_modes` | Count of local maxima in the *smoothed* histogram to the right of `valley_position` whose height ≥ `mode_prominence × global_max_right_of_valley` (default prominence = 0.10). Checked by the strict condition `s[i] > s[i−1] && s[i] > s[i+1]`. |
| 7 | `primary_mode_depth` | argmax of the smoothed histogram for multiplicity > `valley_position`. Equals the expected genome coverage for a haploid sample; multiples of it indicate polyploidy or strain mixing. |
| 8 | `high_mult_tail_ratio` | `Σ hist[i] / total` for i ≥ 3 × `primary_mode_depth`. Captures repeat/mobile-element signal in the high-multiplicity tail. |

---

## Phase 2 — Node-Level Degree Census (features 9–14)

**Source:** `ExtractNodeFeatures(SDBG&, ExtractionOptions&)`

**Method:** Strided sampling. Edges are visited at stride `n / sample_size`
(default sample_size = 100 000) giving an approximately uniform sample across
the edge space. Parallelised with OpenMP `reduction` clauses over six counters
and two floating-point accumulators. Schedule: `dynamic, 256` to handle
variable validity density.

For each sampled valid edge the SDBG provides `EdgeIndegree` and
`EdgeOutdegree` (number of predecessors / successors of that edge's source
node).

| # | Field | Condition | Formula |
|---|-------|-----------|---------|
| 9  | `tip_density` | `outdeg == 0 OR indeg == 0` | `n_tips / n_valid` |
| 10 | `branching_node_fraction` | `indeg >= 2 OR outdeg >= 2` | `n_branch / n_valid` |
| 11 | `linear_node_fraction` | `indeg == 1 AND outdeg == 1` | `n_linear / n_valid` |
| 12 | `high_degree_node_fraction` | `indeg + outdeg >= 4` | `n_high / n_valid` |
| 13 | `mult_at_tips` | tip edges only | `Σ multiplicity / n_tips` |
| 14 | `mult_at_branches` | branching edges only | `Σ multiplicity / n_branch` |

---

## Phase 3 — Tip-Path Walk (feature 15)

**Source:** `ExtractMeanTipLength(SDBG&, ExtractionOptions&)`

**Method:** Same strided sample as Phase 2. Only dead-end forward edges
(`outdeg == 0`, `indeg > 0`) are selected as walk start points. From each
start edge the walk proceeds *backward* via `SDBG::UniquePrevEdge` — which
returns the unique predecessor edge if it exists and the current edge has
in-degree 1, otherwise `kNullID`. The walk terminates on `kNullID` (junction
or chain start) or after `max_tip_length` steps (default 5 000). Parallelised
with `reduction(+:n_tips, total_steps)`, schedule `dynamic, 128`.

| # | Field | Formula |
|---|-------|---------|
| 15 | `mean_tip_length` | `total_steps / n_tips` (edges) |

---

## Phase 4 — Bubble Detection (features 16–18)

**Source:** `ExtractBubbleFeatures(SDBG&, ExtractionOptions&, double valley_position)`

**Method:** Same strided sample. For each sampled edge with out-degree ≥ 2
(`OutgoingEdges` fills an array of up to 8 successors) all pairs of outgoing
branches are tested by `WalkBranchPair`.

**`WalkBranchPair(branch0, branch1, max_depth)`:** advances both cursors
simultaneously using `FirstOutgoingEdge` (takes the first valid outgoing edge
without requiring unique outdegree, allowing the walk to cross high-degree
merge nodes). Convergence is detected when `SDBG::Forward(cur0) == SDBG::Forward(cur1)`,
i.e. both paths arrive at the same node (Forward returns the canonical last-edge
identifier of the target node). The walk budget is `max_bubble_depth` steps
(default 200; SNP bubbles at k = 21 need ≥ 22 edges). Representative
multiplicities are the multiplicities of the two *initial* branch edges.

Parallelised with `reduction` over four counters, schedule `dynamic, 64`
because bubble walk cost is highly variable.

| # | Field | Condition | Formula |
|---|-------|-----------|---------|
| 16 | `bubble_density` | any converging pair | `n_bubbles / n_branching` |
| 17 | `error_bubble_fraction` | `max(m0,m1) / min(m0,m1) > 5` | `n_error / n_bubbles` |
| 18 | `balanced_bubble_fraction` | ratio < 2 AND both arms' mult > `valley_position` | `n_balanced / n_bubbles` |

A ratio > 5 distinguishes sequencing-error branches (one arm is nearly absent)
from true genomic variants. A ratio < 2 with both arms above the valley
identifies heterozygous or strain-level SNP bubbles.

---

## Helper functions

| Function | Purpose |
|----------|---------|
| `GaussianSmooth(hist, σ, radius)` | Convolves a uint64 histogram with a normalised discrete Gaussian kernel; used in Phase 1 with σ = 5, radius = 15. |
| `WalkBranchPair(dbg, b0, b1, max_depth)` | Dual-cursor forward walk; returns `BubbleResult{found, mult0, mult1}`. |
| `FirstOutgoingEdge(dbg, edge)` | Returns the first valid outgoing edge (non-unique outdegree); used inside `WalkBranchPair`. |

---

## Parallelisation summary

| Phase | OpenMP pattern | Schedule |
|-------|---------------|----------|
| 1 — histogram | per-thread private histograms, manual merge | `static` |
| 2 — node census | `reduction` over 6 counters + 2 doubles | `dynamic, 256` |
| 3 — tip walk | `reduction` over 2 counters | `dynamic, 128` |
| 4 — bubble | `reduction` over 4 counters | `dynamic, 64` |

Dynamic scheduling is used for Phases 2–4 because valid-edge density and walk
depth vary across the edge-ID space. Phase 1 uses static scheduling because the
per-edge cost (one multiplicity lookup + histogram increment) is uniform.

---

## Key parameters (`ExtractionOptions`)

| Parameter | Default | Effect |
|-----------|---------|--------|
| `sample_size` | 100 000 | Edges visited in Phases 2–4; stride = `dbg.size() / sample_size`. |
| `max_tip_length` | 5 000 edges | Maximum backward walk length for feature 15. |
| `max_bubble_depth` | 200 edges | Maximum forward walk depth per branch in Phase 4. |
| `mode_prominence` | 0.10 | Minimum height fraction of the global right-of-valley maximum for a local peak to be counted in feature 6. |
