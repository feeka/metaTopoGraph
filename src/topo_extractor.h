#pragma once

#include <cstdint>
#include "sdbg/sdbg.h"
#include "topo_features.h"

// Options controlling sampling depth and walk bounds.
struct ExtractionOptions {
    uint64_t sample_size       = 100000; // number of edges to visit during sampling
    uint64_t max_tip_length    = 5000;   // max steps when walking a tip path
    uint64_t max_bubble_depth  = 5;      // max edge depth per branch in bubble search
};

// Intermediate result for the bubble extraction phase.
struct BubbleFeatures {
    double bubble_density;           // bubbles / sampled branching edges
    double error_bubble_fraction;    // fraction where mult ratio > 5
    double balanced_bubble_fraction; // fraction where ratio < 2 and both > valley
};

// --- Phase 1: histogram features (features 1-8) ---
HistogramFeatures ExtractHistogramFeatures(SDBG& dbg);

// --- Phase 2: node-level features (features 9-14) ---
NodeFeatures ExtractNodeFeatures(SDBG& dbg, const ExtractionOptions& opts);

// --- Phase 3: tip walk (feature 15) ---
// Returns mean tip length in edges over all sampled tip edges.
double ExtractMeanTipLength(SDBG& dbg, const ExtractionOptions& opts);

// --- Phase 4: bubble detection (features 16-18) ---
// valley_position is the histogram valley computed in phase 1; used to
// classify balanced vs. error bubbles.
BubbleFeatures ExtractBubbleFeatures(SDBG& dbg,
                                     const ExtractionOptions& opts,
                                     double valley_position);

// --- Top-level driver ---
TopoFeatures RunExtraction(SDBG& dbg, const ExtractionOptions& opts);
