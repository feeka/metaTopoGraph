#pragma once

#include <cstdint>

// Number of log2-width multiplicity bins:
//   bin 0 → mult == 1
//   bin 1 → mult in [2, 3]
//   bin 2 → mult in [4, 7]
//   bin 3 → mult in [8, 15]
//   bin 4 → mult in [16, 31]
//   bin 5 → mult in [32, 63]
//   bin 6 → mult >= 64
static constexpr int N_BINS = 7;

// Features derived from a full pass over all graph edges.
struct NodeFeatures {
    uint64_t n_kmers;                // total valid k-mer edges in SDBG (both strands)
    double   mult_mean;              // mean edge multiplicity  ≈ coverage
    double   mult_max_ratio;         // mult_max / mult_mean  — normalized headroom
    double   mean_max_branch_ratio;  // mean of max(mult_i / mult_out_j) for branch nodes  ≈ C_real
    double   max_branch_ratio;       // global maximum of that ratio across all branches
    double   mean_jump_ratio;        // mean per branch-node of max(hi,lo)/min(hi,lo) — symmetric
    double   tip_density;            // n_tips / n_kmers
    double   branch_density;         // n_branch_nodes / n_kmers
};

// Features derived from jellyfish k-mer classification (error vs. true).
// All arrays are indexed by multiplicity bin (N_BINS entries).
struct KmerFeatures {
    // Histogram of error k-mers across multiplicity bins.
    double err_hist_node[N_BINS];    // fraction of error k-mers (by count) in each bin
    double err_hist_read[N_BINS];    // fraction of error k-mer occurrences (sum of mult) in each bin

    // P(error | mult in bin) per bin.
    double overlap_density[N_BINS];

    // Mean max-neighbor multiplicity for error k-mers per bin.
    // Neighbor mult is max over all 8 adjacent k-mers present in the reads.
    double err_neighbor_mult[N_BINS];

    // Scalar summaries.
    double mean_error_mult;          // mean multiplicity of error k-mers
    double mean_true_mult;           // mean multiplicity of true k-mers
    double neighbor_contrast_ratio;  // mean(max_nbr_mult of errors) / mean_error_mult
    double isolated_error_frac;      // fraction of error k-mers with no neighbor in reads
};

// Prediction labels: jointly optimal thresholds for classifying error k-mers,
// optimised by F1 maximisation against jellyfish ground truth.
//   Rule: predict_error = (mult < min_count) AND (max_nbr_mult/mult > contrast_threshold)
struct Labels {
    int    min_count;            // optimal multiplicity threshold (predict error if mult < min_count)
    double contrast_threshold;   // optimal neighbor contrast threshold (AND condition)
    double error_fraction;       // fraction of all read k-mers that are sequencing errors
};

// Wall-clock timing (milliseconds).
struct ExtractionTiming {
    double node_ms = 0.0;
    double kmer_ms = 0.0;  // 0 if --ref-fasta not provided
};

// Top-level aggregate returned by RunExtraction.
struct TopoFeatures {
    NodeFeatures     node;
    KmerFeatures     kmer     = {};   // zero-initialised; valid only if has_kmer_features
    Labels           labels   = {};   // zero-initialised; valid only if has_kmer_features
    ExtractionTiming timing;
    bool             has_kmer_features = false;
};
