#pragma once

#include <cstdint>

// Features derived from a full pass over all graph edges.
struct NodeFeatures {
    uint64_t n_prominent_jumps;      // directed edge transitions where mult ratio >= JUMP_THRESHOLD
    double   mult_min;               // minimum edge multiplicity
    double   mult_max;               // maximum edge multiplicity
    double   mult_mean;              // mean edge multiplicity
    double   mean_min_branch_ratio;  // mean of min(mult[i]/mult[out_j]) for out-degree>=2 edges
    double   mean_max_branch_ratio;  // mean of max(mult[i]/mult[out_j]) for out-degree>=2 edges
    uint64_t n_tips;                 // edges with in-degree==0 or out-degree==0
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
