#pragma once

#include <cstdint>

// Features derived from a full pass over all graph edges.
struct NodeFeatures {
    uint64_t n_reads;                // total reads fed to SDBG (0 if not provided)
    uint64_t n_kmers;                // total valid k-mer edges (E_real + E_err)
    double   mult_min;               // minimum edge multiplicity
    double   mult_max;               // maximum edge multiplicity
    double   mult_mean;              // mean edge multiplicity
    double   mean_max_branch_ratio;  // mean of max(mult[i]/mult[out_j]) for out-degree>=2 edges; approx C_real
    double   max_branch_ratio;       // global maximum of (mul_i/mul_out_j) over all branches
    double   mean_jump_ratio;        // mean per branch-node of max(max(mul_i,mul_j)/min(mul_i,mul_j)); symmetric magnitude
    uint64_t n_tips;                 // edges with in-degree==0 or out-degree==0
    uint64_t n_branch_nodes;         // edges with out-degree>=2
};

// Wall-clock timing (milliseconds).
struct ExtractionTiming {
    double node_ms;
};

// Top-level aggregate returned by RunExtraction.
struct TopoFeatures {
    NodeFeatures     node;
    ExtractionTiming timing;
};
