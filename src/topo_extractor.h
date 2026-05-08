#pragma once

#include <cstdint>
#include "sdbg/sdbg.h"
#include "topo_features.h"

// Options (reserved for future use).
struct ExtractionOptions {};

// Extract all features in a single full pass over all graph edges.
NodeFeatures ExtractNodeFeatures(SDBG& dbg);

// Top-level driver.
TopoFeatures RunExtraction(SDBG& dbg, const ExtractionOptions& opts);
