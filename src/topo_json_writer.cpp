#include "topo_json_writer.h"

#include <fstream>
#include <stdexcept>

void WriteFeatures(const TopoFeatures& f, const std::string& output_path) {
    std::ofstream out(output_path.c_str());
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open output file: " + output_path);
    }

    out << "{\n";

    out << "  \"node\": {\n";
    out << "    \"n_prominent_jumps\": "      << f.node.n_prominent_jumps      << ",\n";
    out << "    \"mult_min\": "               << f.node.mult_min               << ",\n";
    out << "    \"mult_max\": "               << f.node.mult_max               << ",\n";
    out << "    \"mult_mean\": "              << f.node.mult_mean              << ",\n";
    out << "    \"mean_min_branch_ratio\": "  << f.node.mean_min_branch_ratio  << ",\n";
    out << "    \"mean_max_branch_ratio\": "  << f.node.mean_max_branch_ratio  << ",\n";
    out << "    \"n_tips\": "                 << f.node.n_tips                 << ",\n";
    out << "    \"error_threshold\": "         << f.node.error_threshold         << "\n";
    out << "  },\n";

    out << "  \"timing\": {\n";
    out << "    \"node_ms\": " << f.timing.node_ms << "\n";
    out << "  }\n";

    out << "}\n";
}
