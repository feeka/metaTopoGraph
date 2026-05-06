#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

#include "sdbg/sdbg.h"
#include "topo_extractor.h"
#include "topo_json_writer.h"

static void PrintUsage(const char* prog) {
    std::cerr << "Usage: " << prog
              << " --graph <sdbg_prefix> --output <features.json>"
                 " [--sample <N>]\n"
              << "\n"
              << "  --graph   Path prefix for the SDBG files written by megahit_core\n"
              << "            (the tool will open <prefix>.sdbg_info and <prefix>.sdbg.*)\n"
              << "  --output  Path for the JSON output file\n"
              << "  --sample  Number of edges to sample (default 100000)\n";
}

int main(int argc, char** argv) {
    std::string graph_prefix;
    std::string output_path;
    ExtractionOptions opts;

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "--graph" && i + 1 < argc) {
            graph_prefix = argv[++i];
        } else if (arg == "--output" && i + 1 < argc) {
            output_path = argv[++i];
        } else if (arg == "--sample" && i + 1 < argc) {
            opts.sample_size = static_cast<uint64_t>(std::atoll(argv[++i]));
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            PrintUsage(argv[0]);
            return 1;
        }
    }

    if (graph_prefix.empty() || output_path.empty()) {
        PrintUsage(argv[0]);
        return 1;
    }

    SDBG dbg;
    try {
        std::cerr << "Loading SDBG from: " << graph_prefix << " ...\n";
        dbg.LoadFromFile(graph_prefix.c_str());
        std::cerr << "  edges: " << dbg.size() << "  k=" << dbg.k() << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Error loading SDBG: " << e.what() << "\n";
        return 1;
    }

    TopoFeatures features;
    try {
        std::cerr << "Extracting features (sample_size=" << opts.sample_size << ") ...\n";
        features = RunExtraction(dbg, opts);
    } catch (const std::exception& e) {
        std::cerr << "Error during extraction: " << e.what() << "\n";
        return 1;
    }

    std::cerr << "  histogram:  " << features.timing.histogram_ms << " ms\n";
    std::cerr << "  node:       " << features.timing.node_ms      << " ms\n";
    std::cerr << "  tip walk:   " << features.timing.walk_ms      << " ms\n";
    std::cerr << "  bubbles:    " << features.timing.bubble_ms    << " ms\n";

    try {
        WriteFeatures(features, output_path);
        std::cerr << "Features written to: " << output_path << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Error writing output: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
