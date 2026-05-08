# Feature Extraction — Implementation Reference

All features are computed from MEGAHIT's Succinct de Bruijn Graph (SDBG) by
`topo_extractor.cpp` in a single pass over all valid edges.

---

## Features

```
n_reads
  Total reads provided (passed in from CLI, 0 if using --graph mode without --read-count).

n_kmers
  Total valid edges in the SDBG. Equals E_real + E_err under the two-population model.
  Combined with mult_mean and mean_max_branch_ratio, E_err can be estimated as:
    E_err = n_kmers * (C_real - mult_mean) / (C_real - 1)

mult_min / mult_max
  Min and max EdgeMultiplicity across all valid edges.

mult_mean
  sum(EdgeMultiplicity) / n_kmers. Pulled below real coverage by error k-mers at
  multiplicity ~1.

mean_min_branch_ratio
  At each branching edge (out-degree >= 2): compute mul_i / mul_out_j for all outgoing
  j, take the minimum. Feature = mean of those minimums. Near 1.0 means the real
  continuation path is always present at branches.

mean_max_branch_ratio
  Same but maximum ratio per branch. Approximates real genome coverage C, because the
  dominant real path consistently outweighs error branches.

max_branch_ratio
  Global maximum of (mul_i / mul_out_j) across all branches. The single most extreme
  multiplicity spike in the graph.

min_branch_ratio
  Global minimum of (mul_i / mul_out_j) across all branches. The gentlest branch —
  where both outgoing paths have the most similar coverage.

n_tips
  Edges with in-degree == 0 or out-degree == 0. Dead-ends caused by errors near read
  boundaries or genuine contig ends.

n_branch_nodes
  Edges with out-degree >= 2. The branching node count.
```

---

## Parallelisation

```
Pattern    : omp parallel for with per-thread ThrAcc structs, manual merge after
Schedule   : dynamic, chunk 512
Accumulators per thread:
  mult_min, mult_max, mult_sum, total
  n_tips, n_branch
  min_ratio_sum, max_ratio_sum
  global_max_r, global_min_r
```
