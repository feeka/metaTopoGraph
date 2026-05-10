# Implementation Notes

## Graph topology — `ExtractNodeFeatures`

Single OpenMP pass over all SDBG edges. For each valid edge, reads multiplicity, in-degree, out-degree. At branching nodes (out-degree ≥ 2), computes ratios between the current edge and each outgoing edge.

| Field | How it's computed |
|---|---|
| `n_kmers` | Count of valid edges |
| `mult_mean` | `sum(mult) / n_kmers` |
| `mult_max_ratio` | `mult_max / mult_mean` |
| `mean_max_branch_ratio` | Mean over branching nodes of `max(mult / mult_out_j)` — approximates true coverage |
| `max_branch_ratio` | Global max of that ratio |
| `mean_jump_ratio` | Mean over branching nodes of `max(hi, lo) / min(hi, lo)` between edge and its best outgoing neighbour |
| `tip_density` | `n_tips / n_kmers` |
| `branch_density` | `n_branch_nodes / n_kmers` |

Parallelisation: `omp parallel for schedule(dynamic, 512)` with per-thread accumulators, merged after.

---

## K-mer features — `ExtractKmerFeatures`

Runs jellyfish via `system()` calls, then classifies k-mers and accumulates statistics.

### Jellyfish steps

```bash
jellyfish count -m K -s 200M -C -t T -o ref.jf   <reference>
jellyfish count -m K -s 100M -C -t T -o reads.jf  <reads>
jellyfish dump reads.jf                            > reads_dump.fa
jellyfish query ref.jf -s reads_dump.fa            > query.txt
```

A k-mer is **error** if its count in `query.txt` is 0 (not in reference). All files go into `opts.tmp_dir` and are deleted by `main.cpp` on exit.

### Multiplicity bins

```
bin 0: mult == 1
bin 1: mult in [2, 3]
bin 2: mult in [4, 7]
bin 3: mult in [8, 15]
bin 4: mult in [16, 31]
bin 5: mult in [32, 63]
bin 6: mult >= 64
```

Assigned by `mult_bin(m)`: right-shift m until 1, count shifts, cap at N_BINS-1.

### Neighbour lookup

For a canonical k-mer `s`, the 8 neighbours are:
- Right extensions: `canonical(s[1:] + X)` for X in {A,C,G,T}
- Left extensions: `canonical(X + s[:-1])` for X in {A,C,G,T}

Each is looked up in the in-memory dump map. Canonical = lexicographically smaller of kmer and its reverse complement, matching jellyfish `-C`.

---

## Labels — threshold search

### `min_count`

Build a per-multiplicity histogram of error and total counts. Scan τ from 2 upward — at each step, accept all k-mers with `mult == τ-1` as predicted errors. Track TP/FP/FN, compute F1, keep the best τ. O(max_mult).

### `contrast_threshold`

Fix `min_count`. Take only k-mers with `mult < min_count`. For each, compute `contrast = max_neighbour_mult / mult`. Sort by contrast ascending. Walk through the list raising the threshold — each step removes one k-mer from the predicted-error set. Track F1, keep best threshold. O(n log n).

Note: the two thresholds are optimised sequentially, not jointly. This is intentional — it's fast and per-sample labels don't need to be globally optimal, just informative training targets.
