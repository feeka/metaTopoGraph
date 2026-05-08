# Feature Extraction — Implementation Reference

All 7 features are computed from MEGAHIT's Succinct de Bruijn Graph (SDBG) by
`topo_extractor.cpp` in a single full scan. `ExtractionOptions` is empty (no
parameters); features and timing are returned in `TopoFeatures`.

---

## Single-pass full scan

**Source:** `ExtractNodeFeatures(SDBG&)`

**Method:** Every valid edge in the SDBG is visited exactly once. The scan is
parallelised with `#pragma omp parallel for` and per-thread accumulator structs
that are merged after the parallel section to avoid false sharing on counters
and floating-point accumulators. Schedule: `static` — per-edge cost is uniform
(one multiplicity lookup + degree queries).

For each valid edge, `EdgeMultiplicity`, `EdgeIndegree`, and `EdgeOutdegree`
are queried, plus `OutgoingEdges` for branching edges.

```
n_prominent_jumps
  Condition : for every outgoing edge j of edge i,
              max(mul_i, mul_j) / min(mul_i, mul_j) >= JUMP_THRESHOLD (3.0)
  Increment : +1 per (i, j) pair that meets the condition
  Meaning   : counts directed transitions between adjacent edges where
              multiplicity changes abruptly. At high coverage, most jumps
              are real-path vs error-branch transitions.

mult_min
  The minimum EdgeMultiplicity across all valid edges.

mult_max
  The maximum EdgeMultiplicity across all valid edges.
  High values indicate repeated k-mers (rRNA operons, conserved genes).

mult_mean
  sum(EdgeMultiplicity) / n_valid_edges
  Pulled below true coverage by error k-mers at multiplicity ~1.
  Together with mean_max_branch_ratio, allows Bayesian estimation of
  the per-base error rate with no external calibration.

mean_min_branch_ratio
  Condition : edge has out-degree >= 2
  Per edge  : compute mul[i] / mul[out_j] for every outgoing j,
              record the minimum ratio
  Feature   : mean of those minimums across all branching edges
  Meaning   : near 1.0 means the real continuation path is always
              present at branches (real path never drops out).

mean_max_branch_ratio
  Condition : edge has out-degree >= 2
  Per edge  : compute mul[i] / mul[out_j] for every outgoing j,
              record the maximum ratio
  Feature   : mean of those maximums across all branching edges
  Meaning   : approximates real genome coverage C, because the dominant
              real path consistently wins over error branches (mul ~1).
              Directly usable as lambda_r in the Bayesian threshold formula.

n_tips
  Condition : EdgeIndegree == 0 OR EdgeOutdegree == 0
  Meaning   : dead-end edges. At low coverage (lambda_r < 5) these are
              mostly coverage gaps; at high coverage mostly sequencing
              errors near read boundaries.
```

---

## Bayesian error threshold

Given only the JSON output and the total read count, the multiplicity threshold
$m^*$ below which an edge is more likely erroneous than real is:

```
lambda_r  = mean_max_branch_ratio              (real genome coverage)
ratio     = (lambda_r - mult_mean) / (mult_mean - 1)   (E_err / E_real)
lambda_e  = lambda_r * p / 3                   (error k-mer Poisson mean)
  where p = ratio / (21 * lambda_r)

m* = ( ln(ratio) + (lambda_r - lambda_e) ) / ln(lambda_r / lambda_e)
```

Edges with multiplicity <= floor(m*) are more likely erroneous. This threshold
is only meaningful for lambda_r > ~5; below that, real and error k-mers are
indistinguishable by multiplicity alone.

---

## Parallelisation

```
OpenMP pattern  : parallel for with per-thread accumulator structs, manual merge
Schedule        : static (uniform per-edge cost)
Accumulators    : n_prominent_jumps (uint64), n_tips (uint64),
                  mult_sum (double), mult_min (double), mult_max (double),
                  branch_min_ratio_sum (double), branch_max_ratio_sum (double),
                  n_branching (uint64)
```

---

## Constants

```
JUMP_THRESHOLD  3.0   defined in topo_extractor.h
                      ratio between adjacent edge multiplicities that
                      counts as a prominent jump
```
