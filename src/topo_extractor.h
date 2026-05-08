#pragma once

#include <cstdint>
#include "sdbg/sdbg.h"
#include "topo_features.h"

// Multiplicity ratio between adjacent edges that counts as a prominent jump.
constexpr double JUMP_THRESHOLD = 3.0;

// Options (reserved for future use).
struct ExtractionOptions {};

// Extract all features in a single full pass over all graph edges.
NodeFeatures ExtractNodeFeatures(SDBG& dbg);

// Top-level driver.
TopoFeatures RunExtraction(SDBG& dbg, const ExtractionOptions& opts);
