# Feature Extraction — Implementation Reference

All features are computed by `topo_extractor.cpp`.
Graph topology features come from a single parallel pass over all SDBG edges.
K-mer error features are computed by shelling out to **jellyfish 2** (must be on PATH).

---

## Graph Topology Features (`NodeFeatures`)

Computed by `ExtractNodeFeatures(SDBG& dbg, uint64_t total_reads)`.

```
n_kmers
  Total valid edges in the SDBG. Both strands are stored separately, so this equals
  2 × canonical k-mer count from jellyfish -C.

mult_mean
  mean(EdgeMultiplicity) over all valid edges. Pulled below the true genome coverage
  by low-multiplicity error k-mers.

mult_max_ratio
  mult_max / mult_mean — normalised coverage headroom. Invariant to sequencing depth.

mean_max_branch_ratio
  At each branching edge (out-degree ≥ 2): compute mult_i / mult_out_j for all
  outgoing j, take the maximum. Feature = mean of those maxima across all branches.
  Approximates the real genome coverage C because the dominant true-path consistently
  outweighs error branches.

max_branch_ratio
  Global maximum of (mult_i / mult_out_j) across all branches. Single most extreme
  multiplicity spike — captures the worst local artefact.

mean_jump_ratio
  Mean per branch-node of max(hi, lo) / min(hi, lo) where hi/lo are the edge
  multiplicity and the max outgoing multiplicity. Symmetric; robust to direction.

tip_density
  n_tips / n_kmers, where tips are edges with in-degree == 0 or out-degree == 0.
  Dead-ends caused by errors near read boundaries or genuine contig terminations.

branch_density
  n_branch_nodes / n_kmers. Density of branching points in the graph.
```

---

## K-mer Error Features (`KmerFeatures`)

Computed by `ExtractKmerFeatures(const ExtractionOptions& opts, Labels& labels_out)`.
Requires `opts.ref_fasta`, `opts.reads1`, and `opts.tmp_dir` to be set.

### Jellyfish workflow (all run as system() calls)

```
jellyfish count -m K -s 200M -C -t T -o ref.jf   <reference_fasta>
jellyfish count -m K -s 100M -C -t T -o reads.jf  <reads1> [reads2]
jellyfish dump reads.jf                            → reads_dump.fa  (FASTA: >count / KMER)
jellyfish query ref.jf -s reads_dump.fa            → query.txt      (KMER count_in_ref)
```

K-mer is labelled **error** if count_in_ref == 0, **true** otherwise.

All intermediate files are written to `opts.tmp_dir` and cleaned up by main.cpp on exit.

### Multiplicity bins

```
bin 0 → mult == 1
bin 1 → mult in [2, 3]
bin 2 → mult in [4, 7]
bin 3 → mult in [8, 15]
bin 4 → mult in [16, 31]
bin 5 → mult in [32, 63]
bin 6 → mult ≥ 64
```

### Features (32 total)

```
err_hist_node[7]
  Fraction of error k-mers (by count) in each bin. Node-weighted histogram.

err_hist_read[7]
  Fraction of error k-mer occurrences (Σ mult) in each bin. Read-weighted histogram.

overlap_density[7]
  P(error | mult in bin) = n_errors_in_bin / n_total_in_bin per bin.

err_neighbor_mult[7]
  Mean max-neighbour multiplicity for error k-mers per bin.
  Neighbour multiplicity = max mult over all 8 adjacent canonical k-mers found in reads.

mean_error_mult        mean multiplicity of error k-mers
mean_true_mult         mean multiplicity of true k-mers
neighbor_contrast_ratio  mean(max_nbr_mult of errors) / mean_error_mult
isolated_error_frac    fraction of error k-mers with no neighbour present in reads
```

### Neighbour lookup

Given a canonical k-mer `s` of length k, the 8 candidate canonical neighbours are:
- right extensions: `canonical(s[1:] + X)` for X in {A,C,G,T}
- left  extensions: `canonical(X + s[:-1])` for X in {A,C,G,T}

Each is looked up in the in-memory dump map. This correctly handles both strands because
canonical() chooses the lexicographically smaller of the k-mer and its reverse complement,
matching jellyfish -C behaviour.

---

## Labels (`Labels`)

Optimal thresholds jointly maximising F1 against jellyfish ground truth.

```
min_count
  Optimal multiplicity threshold τ₁.
  Predicted error rule (single threshold): mult < τ₁.
  Chosen by scanning τ₁ ∈ {2, …, mult_max+1} and picking the one with highest F1.

contrast_threshold
  Optimal neighbour contrast threshold τ₂ (AND condition).
  Predicted error rule (joint): mult < τ₁ AND (max_nbr_mult / mult) > τ₂.
  Chosen by sorting low-mult k-mers by contrast ascending, scanning τ₂ upward,
  and picking the one with highest F1 over all k-mers.

error_fraction
  Fraction of all unique canonical read k-mers that are sequencing errors.
```

---

## Parallelisation

Graph topology extraction uses OpenMP:

```
Pattern    : omp parallel for with per-thread ThrAcc structs, manual merge after
Schedule   : dynamic, chunk 512
Thread accumulators: mult_min, mult_max, mult_sum, total, n_tips, n_branch,
                     max_ratio_sum, jump_ratio_sum, global_max_r, global_max_mag
```

Jellyfish extraction is single-threaded in C++ (I/O-bound); jellyfish itself uses
`opts.n_threads` threads internally via its `-t` flag.
