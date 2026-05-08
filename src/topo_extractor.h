#pragma once

#include <cstdint>
#include "sdbg/sdbg.h"
#include "topo_features.h"

struct ExtractionOptions {
    uint64_t total_reads = 0;  // total read count (R1+R2); 0 = not provided
};

// Extract all features in two passes over all graph edges.
NodeFeatures ExtractNodeFeatures(SDBG& dbg, uint64_t total_reads);

// Top-level driver.
TopoFeatures RunExtraction(SDBG& dbg, const ExtractionOptions& opts);
