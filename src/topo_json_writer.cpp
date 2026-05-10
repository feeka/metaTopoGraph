#include "topo_json_writer.h"

#include <fstream>
#include <stdexcept>

// Helper: write a JSON array of N_BINS doubles on one line.
static void WriteArray(std::ofstream& out, const double* arr, int n) {
    out << "[";
    for (int i = 0; i < n; ++i) {
        if (i > 0) out << ", ";
        out << arr[i];
    }
    out << "]";
}

void WriteFeatures(const TopoFeatures& f, const std::string& output_path) {
    std::ofstream out(output_path.c_str());
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open output file: " + output_path);
    }

    out << "{\n";

    // ---- node (graph topology) features ----
    out << "  \"node\": {\n";
    out << "    \"n_kmers\": "               << f.node.n_kmers               << ",\n";
    out << "    \"mult_mean\": "             << f.node.mult_mean             << ",\n";
    out << "    \"mult_max_ratio\": "        << f.node.mult_max_ratio        << ",\n";
    out << "    \"mean_max_branch_ratio\": " << f.node.mean_max_branch_ratio << ",\n";
    out << "    \"max_branch_ratio\": "      << f.node.max_branch_ratio      << ",\n";
    out << "    \"mean_jump_ratio\": "       << f.node.mean_jump_ratio       << ",\n";
    out << "    \"tip_density\": "           << f.node.tip_density           << ",\n";
    out << "    \"branch_density\": "        << f.node.branch_density        << "\n";
    out << "  },\n";

    // ---- k-mer features (only present when --ref-fasta was supplied) ----
    if (f.has_kmer_features) {
        const KmerFeatures& k = f.kmer;
        out << "  \"kmer_features\": {\n";
        out << "    \"err_hist_node\": ";        WriteArray(out, k.err_hist_node,      N_BINS); out << ",\n";
        out << "    \"err_hist_read\": ";        WriteArray(out, k.err_hist_read,      N_BINS); out << ",\n";
        out << "    \"overlap_density\": ";      WriteArray(out, k.overlap_density,    N_BINS); out << ",\n";
        out << "    \"err_neighbor_mult\": ";    WriteArray(out, k.err_neighbor_mult,  N_BINS); out << ",\n";
        out << "    \"mean_error_mult\": "       << k.mean_error_mult         << ",\n";
        out << "    \"mean_true_mult\": "        << k.mean_true_mult          << ",\n";
        out << "    \"neighbor_contrast_ratio\": " << k.neighbor_contrast_ratio << ",\n";
        out << "    \"isolated_error_frac\": "  << k.isolated_error_frac     << "\n";
        out << "  },\n";

        // ---- labels (prediction targets) ----
        out << "  \"labels\": {\n";
        out << "    \"min_count\": "          << f.labels.min_count          << ",\n";
        out << "    \"contrast_threshold\": " << f.labels.contrast_threshold << ",\n";
        out << "    \"error_fraction\": "     << f.labels.error_fraction     << "\n";
        out << "  },\n";
    }

    // ---- timing ----
    out << "  \"timing\": {\n";
    out << "    \"node_ms\": " << f.timing.node_ms;
    if (f.has_kmer_features)
        out << ",\n    \"kmer_ms\": " << f.timing.kmer_ms << "\n";
    else
        out << "\n";
    out << "  }\n";

    out << "}\n";
}
