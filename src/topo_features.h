#pragma once

#include <cstdint>

// Features derived from the edge multiplicity histogram (features 1-8).
struct HistogramFeatures {
    double valley_position;         // argmin of smoothed hist between mult 2..50
    double valley_depth;            // count_at_valley / count_at_error_peak
    double mult_1_fraction;         // hist[1] / total_edges
    double mean_node_multiplicity;  // weighted mean of multiplicity histogram
    double multiplicity_cv;         // std / mean of multiplicity
    double n_signal_modes;          // peak count in smoothed hist right of valley
    double primary_mode_depth;      // argmax right of valley
    double high_mult_tail_ratio;    // sum(hist[3*primary:]) / total
};

// Features derived from sampled node-level degree queries (features 9-14).
struct NodeFeatures {
    double tip_density;               // nodes where one-side degree == 0
    double branching_node_fraction;   // nodes with in-deg >= 2 OR out-deg >= 2
    double linear_node_fraction;      // nodes with in-deg == 1 AND out-deg == 1
    double high_degree_node_fraction; // nodes with in-deg + out-deg >= 4
    double mult_at_tips;              // mean multiplicity of tip edges
    double mult_at_branches;          // mean multiplicity of branching edges
};

// Features derived from bounded local walks (features 15-18).
struct WalkFeatures {
    double mean_tip_length;          // mean walk length from sampled tip edges
    double bubble_density;           // bubbles found per sampled branching edge
    double error_bubble_fraction;    // bubble fraction where mult ratio > 5
    double balanced_bubble_fraction; // bubble fraction where ratio < 2 and both > valley
};

// Wall-clock timing for each extraction phase (milliseconds).
struct ExtractionTiming {
    double histogram_ms;
    double node_ms;
    double walk_ms;
    double bubble_ms;
};

// Top-level aggregate returned by RunExtraction.
struct TopoFeatures {
    HistogramFeatures hist;
    NodeFeatures      node;
    WalkFeatures      walk;
    ExtractionTiming  timing;
};
