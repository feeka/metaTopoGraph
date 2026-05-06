#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "sdbg/sdbg.h"
#include "topo_extractor.h"
#include "topo_json_writer.h"

static void PrintUsage(const char* prog) {
    std::cerr << "Usage: " << prog
              << " --graph <sdbg_prefix> --output <features.json>"
                 " [--sample <N>] [--threads <T>]\n"
              << "\n"
              << "  --graph    Path prefix for the SDBG files written by megahit_core\n"
              << "             (the tool will open <prefix>.sdbg_info and <prefix>.sdbg.*)\n"
              << "  --output   Path for the JSON output file\n"
              << "  --sample   Number of edges to sample (default 100000)\n"
              << "  --threads  Number of OpenMP threads (default: all available)\n";
}

int main(int argc, char** argv) {
    std::string graph_prefix;
    std::string output_path;
    ExtractionOptions opts;
    int n_threads = 0;  // 0 = use OpenMP default

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "--graph" && i + 1 < argc) {
            graph_prefix = argv[++i];
        } else if (arg == "--output" && i + 1 < argc) {
            output_path = argv[++i];
        } else if (arg == "--sample" && i + 1 < argc) {
            opts.sample_size = static_cast<uint64_t>(std::atoll(argv[++i]));
        } else if (arg == "--threads" && i + 1 < argc) {
            n_threads = std::atoi(argv[++i]);
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            PrintUsage(argv[0]);
            return 1;
        }
    }

#ifdef _OPENMP
    if (n_threads > 0) {
        omp_set_num_threads(n_threads);
    }
    std::cerr << "Threads: " << omp_get_max_threads() << "\n";
#else
    std::cerr << "Threads: 1 (OpenMP not available)\n";
#endif

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
