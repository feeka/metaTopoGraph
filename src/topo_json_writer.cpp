#include "topo_json_writer.h"

#include <fstream>
#include <stdexcept>

void WriteFeatures(const TopoFeatures& f, const std::string& output_path) {
    std::ofstream out(output_path.c_str());
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open output file: " + output_path);
    }

    out << "{\n";

    out << "  \"histogram\": {\n";
    out << "    \"valley_position\": "        << f.hist.valley_position        << ",\n";
    out << "    \"valley_depth\": "           << f.hist.valley_depth           << ",\n";
    out << "    \"mult_1_fraction\": "        << f.hist.mult_1_fraction        << ",\n";
    out << "    \"mean_node_multiplicity\": " << f.hist.mean_node_multiplicity << ",\n";
    out << "    \"multiplicity_cv\": "        << f.hist.multiplicity_cv        << ",\n";
    out << "    \"n_signal_modes\": "         << f.hist.n_signal_modes         << ",\n";
    out << "    \"primary_mode_depth\": "     << f.hist.primary_mode_depth     << ",\n";
    out << "    \"high_mult_tail_ratio\": "   << f.hist.high_mult_tail_ratio   << "\n";
    out << "  },\n";

    out << "  \"node\": {\n";
    out << "    \"tip_density\": "                  << f.node.tip_density                  << ",\n";
    out << "    \"branching_node_fraction\": "      << f.node.branching_node_fraction      << ",\n";
    out << "    \"linear_node_fraction\": "         << f.node.linear_node_fraction         << ",\n";
    out << "    \"high_degree_node_fraction\": "    << f.node.high_degree_node_fraction    << ",\n";
    out << "    \"mult_at_tips\": "                 << f.node.mult_at_tips                 << ",\n";
    out << "    \"mult_at_branches\": "             << f.node.mult_at_branches             << "\n";
    out << "  },\n";

    out << "  \"walk\": {\n";
    out << "    \"mean_tip_length\": "           << f.walk.mean_tip_length           << ",\n";
    out << "    \"bubble_density\": "            << f.walk.bubble_density            << ",\n";
    out << "    \"error_bubble_fraction\": "     << f.walk.error_bubble_fraction     << ",\n";
    out << "    \"balanced_bubble_fraction\": "  << f.walk.balanced_bubble_fraction  << "\n";
    out << "  },\n";

    out << "  \"timing\": {\n";
    out << "    \"histogram_ms\": " << f.timing.histogram_ms << ",\n";
    out << "    \"node_ms\": "      << f.timing.node_ms      << ",\n";
    out << "    \"walk_ms\": "      << f.timing.walk_ms      << ",\n";
    out << "    \"bubble_ms\": "    << f.timing.bubble_ms    << "\n";
    out << "  }\n";

    out << "}\n";
}
