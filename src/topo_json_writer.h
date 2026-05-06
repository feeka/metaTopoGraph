#pragma once

#include <string>
#include "topo_features.h"

// Write all 18 features to a JSON file at output_path.
// Throws std::runtime_error if the file cannot be opened.
void WriteFeatures(const TopoFeatures& features, const std::string& output_path);
