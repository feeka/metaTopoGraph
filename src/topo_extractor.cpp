#include "topo_extractor.h"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "sdbg/sdbg_def.h"

// ---------------------------------------------------------------------------
// Full-scan node feature extraction (single pass)
// ---------------------------------------------------------------------------

NodeFeatures ExtractNodeFeatures(SDBG& dbg, uint64_t total_reads) {

#ifdef _OPENMP
    const int nthreads = omp_get_max_threads();
#else
    const int nthreads = 1;
#endif

    struct ThrAcc {
        double   mult_min      =  std::numeric_limits<double>::max();
        double   mult_max      = -std::numeric_limits<double>::max();
        double   mult_sum      = 0.0;
        uint64_t total         = 0;
        uint64_t n_tips        = 0;
        uint64_t n_branch      = 0;
        double   max_ratio_sum  = 0.0;
        double   jump_ratio_sum = 0.0;
        double   global_max_r   = -std::numeric_limits<double>::max();
        double   global_max_mag = -std::numeric_limits<double>::max();
    };
    std::vector<ThrAcc> acc(nthreads);

    // ---- Pass 1: all stats except n_prominent_jumps ----
#pragma omp parallel for schedule(dynamic, 512)
    for (uint64_t i = 0; i < dbg.size(); ++i) {
        if (!dbg.IsValidEdge(i)) continue;
#ifdef _OPENMP
        const int tid = omp_get_thread_num();
#else
        const int tid = 0;
#endif
        ThrAcc& a = acc[tid];

        const double mul    = dbg.EdgeMultiplicity(i);
        const int    indeg  = dbg.EdgeIndegree(i);
        const int    outdeg = dbg.EdgeOutdegree(i);

        a.mult_sum += mul;
        ++a.total;
        if (mul < a.mult_min) a.mult_min = mul;
        if (mul > a.mult_max) a.mult_max = mul;

        if (indeg == 0 || outdeg == 0) ++a.n_tips;

        if (outdeg >= 2) {
            uint64_t outs[8];
            int od = dbg.OutgoingEdges(i, outs);
            double max_r   = -std::numeric_limits<double>::max();
            double max_mag = -std::numeric_limits<double>::max();
            for (int j = 0; j < od; ++j) {
                const double mul_j = dbg.EdgeMultiplicity(outs[j]);
                const double r     = mul_j > 0.0 ? mul / mul_j : mul;
                const double hi    = (mul >= mul_j) ? mul : mul_j;
                const double lo    = (mul >= mul_j) ? mul_j : mul;
                const double mag   = (lo > 0.0) ? hi / lo : hi;
                if (r   > max_r)   max_r   = r;
                if (mag > max_mag) max_mag = mag;
            }
            a.max_ratio_sum  += max_r;
            a.jump_ratio_sum += max_mag;
            if (max_r   > a.global_max_r)   a.global_max_r   = max_r;
            if (max_mag > a.global_max_mag) a.global_max_mag = max_mag;
            ++a.n_branch;
        }
    }

    // Merge pass 1
    double   mult_min  =  std::numeric_limits<double>::max();
    double   mult_max  = -std::numeric_limits<double>::max();
    double   mult_sum  = 0.0;
    uint64_t total     = 0;
    uint64_t n_tips    = 0;
    uint64_t n_branch  = 0;
    double   max_r_sum    = 0.0;
    double   jump_r_sum   = 0.0;
    double   global_max_r   = -std::numeric_limits<double>::max();
    double   global_max_mag = -std::numeric_limits<double>::max();

    for (const ThrAcc& a : acc) {
        if (a.total == 0) continue;
        if (a.mult_min < mult_min) mult_min = a.mult_min;
        if (a.mult_max > mult_max) mult_max = a.mult_max;
        mult_sum  += a.mult_sum;
        total     += a.total;
        n_tips    += a.n_tips;
        n_branch     += a.n_branch;
        max_r_sum    += a.max_ratio_sum;
        jump_r_sum   += a.jump_ratio_sum;
        if (a.global_max_r   > global_max_r)   global_max_r   = a.global_max_r;
        if (a.global_max_mag > global_max_mag) global_max_mag = a.global_max_mag;
    }

    if (total == 0) return NodeFeatures{};

    const double mean_mult = mult_sum / total;
    const double C_real    = n_branch > 0 ? max_r_sum / n_branch : 0.0;

    NodeFeatures nf        = {};
    nf.n_kmers             = total;
    nf.mult_mean           = mean_mult;
    nf.mult_max_ratio      = (mean_mult > 0.0 && mult_max != -std::numeric_limits<double>::max())
                             ? mult_max / mean_mult : 0.0;
    nf.mean_max_branch_ratio = C_real;
    nf.max_branch_ratio    = (n_branch > 0 && global_max_r != -std::numeric_limits<double>::max())
                             ? global_max_r : 0.0;
    nf.mean_jump_ratio     = (n_branch > 0 && global_max_mag != -std::numeric_limits<double>::max())
                             ? jump_r_sum / n_branch : 0.0;
    nf.tip_density         = total > 0 ? static_cast<double>(n_tips)   / total : 0.0;
    nf.branch_density      = total > 0 ? static_cast<double>(n_branch) / total : 0.0;
    return nf;
}

// ---------------------------------------------------------------------------
// DNA utilities for jellyfish k-mer neighbour lookup
// ---------------------------------------------------------------------------

namespace {

static char dna_comp(char c) {
    switch (c) {
        case 'A': return 'T'; case 'T': return 'A';
        case 'C': return 'G'; case 'G': return 'C';
        case 'a': return 't'; case 't': return 'a';
        case 'c': return 'g'; case 'g': return 'c';
        default:  return 'N';
    }
}

static std::string dna_revcomp(const std::string& s) {
    const size_t n = s.size();
    std::string r(n, 'N');
    for (size_t i = 0; i < n; ++i)
        r[i] = dna_comp(s[n - 1 - i]);
    return r;
}

static std::string dna_canonical(const std::string& s) {
    std::string rc = dna_revcomp(s);
    return s < rc ? s : rc;
}

// Log2-floor bin index in [0, N_BINS-1] for multiplicity m >= 1.
static int mult_bin(uint32_t m) {
    int b = 0;
    uint32_t v = m;
    while (v > 1 && b < N_BINS - 1) { v >>= 1; ++b; }
    return b;
}

// Run a shell command; throw on non-zero exit.
static void RunShell(const std::string& cmd, const char* desc) {
    int ret = system(cmd.c_str());
    if (ret != 0)
        throw std::runtime_error(
            std::string(desc) + " failed (exit " + std::to_string(ret) + ")");
}

// Load a jellyfish FASTA dump (>count\nKMER\n ...) into a kmer→count map.
static std::unordered_map<std::string, uint32_t>
LoadJfDump(const std::string& path) {
    std::unordered_map<std::string, uint32_t> m;
    FILE* f = fopen(path.c_str(), "r");
    if (!f) throw std::runtime_error("Cannot open jellyfish dump: " + path);
    char hdr[64], seq[256];
    while (fgets(hdr, sizeof(hdr), f)) {
        if (hdr[0] != '>') continue;
        if (!fgets(seq, sizeof(seq), f)) break;
        uint32_t cnt = (uint32_t)std::atoi(hdr + 1);
        std::string kmer(seq);
        while (!kmer.empty() && (kmer.back() == '\n' || kmer.back() == '\r'))
            kmer.pop_back();
        if (!kmer.empty()) m[kmer] = cnt;
    }
    fclose(f);
    return m;
}

// Return canonical forms of all 8 potential neighbours of a canonical k-mer.
static std::vector<std::string> KmerNeighbors(const std::string& kmer) {
    const size_t k = kmer.size();
    const char bases[] = {'A', 'C', 'G', 'T'};
    std::vector<std::string> nbrs;
    nbrs.reserve(8);
    const std::string suf = kmer.substr(1);
    const std::string pre = kmer.substr(0, k - 1);
    for (char b : bases) nbrs.push_back(dna_canonical(suf + b));
    for (char b : bases) nbrs.push_back(dna_canonical(std::string(1, b) + pre));
    return nbrs;
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// K-mer feature extraction via jellyfish
// ---------------------------------------------------------------------------

KmerFeatures ExtractKmerFeatures(const ExtractionOptions& opts, Labels& labels_out) {
    const std::string& ref_fasta = opts.ref_fasta;
    const std::string& reads1    = opts.reads1;
    const std::string& reads2    = opts.reads2;
    const std::string& tmp_dir   = opts.tmp_dir;
    const int k = opts.kmer_k;
    const int t = opts.n_threads;

    const std::string ref_jf    = tmp_dir + "/ref.jf";
    const std::string reads_jf  = tmp_dir + "/reads.jf";
    const std::string dump_fa   = tmp_dir + "/reads_dump.fa";
    const std::string query_out = tmp_dir + "/query.txt";

    // ---- jellyfish count reference ----
    {
        std::ostringstream cmd;
        cmd << "jellyfish count -m " << k << " -s 200M -t " << t
            << " -C -o \"" << ref_jf << "\" \"" << ref_fasta << "\"";
        std::cerr << "jf count ref: " << cmd.str() << "\n";
        RunShell(cmd.str(), "jellyfish count ref");
    }
    // ---- jellyfish count reads ----
    {
        std::ostringstream cmd;
        cmd << "jellyfish count -m " << k << " -s 100M -t " << t
            << " -C -o \"" << reads_jf << "\" \"" << reads1 << "\"";
        if (!reads2.empty()) cmd << " \"" << reads2 << "\"";
        std::cerr << "jf count reads: " << cmd.str() << "\n";
        RunShell(cmd.str(), "jellyfish count reads");
    }
    // ---- jellyfish dump reads k-mers ----
    {
        std::ostringstream cmd;
        cmd << "jellyfish dump \"" << reads_jf << "\" > \"" << dump_fa << "\"";
        RunShell(cmd.str(), "jellyfish dump");
    }
    // ---- load dump into memory ----
    std::cerr << "Loading jellyfish dump...\n";
    auto kmer_mult = LoadJfDump(dump_fa);
    std::cerr << "  " << kmer_mult.size() << " unique k-mers in reads\n";

    // ---- jellyfish query reads k-mers against reference ----
    {
        std::ostringstream cmd;
        cmd << "jellyfish query \"" << ref_jf << "\" -s \"" << dump_fa
            << "\" > \"" << query_out << "\"";
        RunShell(cmd.str(), "jellyfish query");
    }
    // ---- parse query output: kmer → is_true ----
    std::unordered_map<std::string, bool> is_true;
    is_true.reserve(kmer_mult.size());
    {
        FILE* f = fopen(query_out.c_str(), "r");
        if (!f)
            throw std::runtime_error("Cannot open jellyfish query output: " + query_out);
        char line[512];
        while (fgets(line, sizeof(line), f)) {
            std::string l(line);
            size_t sp = l.find(' ');
            if (sp == std::string::npos) continue;
            std::string kmer    = l.substr(0, sp);
            uint32_t ref_cnt    = (uint32_t)std::atoi(l.c_str() + sp + 1);
            is_true[kmer]       = (ref_cnt > 0);
        }
        fclose(f);
    }

    // ---- per-k-mer: classify and get max neighbour multiplicity ----
    struct Entry {
        uint32_t mult;
        bool     is_error;
        uint32_t max_nbr_mult;   // 0 if no neighbour present in reads
    };
    std::vector<Entry> entries;
    entries.reserve(kmer_mult.size());

    for (auto& kv : kmer_mult) {
        const std::string& kmer = kv.first;
        uint32_t           mult = kv.second;
        auto it = is_true.find(kmer);
        bool err = (it == is_true.end() || !it->second);

        uint32_t max_nbr = 0;
        for (auto& nbr : KmerNeighbors(kmer)) {
            auto jt = kmer_mult.find(nbr);
            if (jt != kmer_mult.end() && jt->second > max_nbr)
                max_nbr = jt->second;
        }
        entries.push_back({mult, err, max_nbr});
    }
    // free memory no longer needed
    is_true.clear();
    kmer_mult.clear();

    // ---- accumulate bin statistics ----
    KmerFeatures kf = {};

    uint64_t bin_err_node  [N_BINS] = {};
    uint64_t bin_err_read  [N_BINS] = {};
    uint64_t bin_total_node[N_BINS] = {};
    uint64_t bin_total_read[N_BINS] = {};
    double   bin_nbr_sum   [N_BINS] = {};
    uint64_t bin_nbr_cnt   [N_BINS] = {};

    double   sum_err_mult  = 0.0;
    double   sum_true_mult = 0.0;
    uint64_t n_err         = 0;
    uint64_t n_true        = 0;
    double   sum_nbr_err   = 0.0;
    uint64_t n_nbr_err     = 0;
    uint64_t n_isolated    = 0;

    for (const Entry& e : entries) {
        int b = mult_bin(e.mult);
        bin_total_node[b]++;
        bin_total_read[b] += e.mult;
        if (e.is_error) {
            bin_err_node[b]++;
            bin_err_read[b] += e.mult;
            sum_err_mult += e.mult;
            ++n_err;
            if (e.max_nbr_mult > 0) {
                bin_nbr_sum[b] += e.max_nbr_mult;
                ++bin_nbr_cnt[b];
                sum_nbr_err += e.max_nbr_mult;
                ++n_nbr_err;
            } else {
                ++n_isolated;
            }
        } else {
            sum_true_mult += e.mult;
            ++n_true;
        }
    }

    uint64_t total_err_read = 0;
    for (int b = 0; b < N_BINS; ++b) total_err_read += bin_err_read[b];

    for (int b = 0; b < N_BINS; ++b) {
        kf.err_hist_node[b]     = n_err         > 0 ? (double)bin_err_node[b]  / n_err         : 0.0;
        kf.err_hist_read[b]     = total_err_read > 0 ? (double)bin_err_read[b]  / total_err_read : 0.0;
        kf.overlap_density[b]   = bin_total_node[b] > 0 ? (double)bin_err_node[b] / bin_total_node[b] : 0.0;
        kf.err_neighbor_mult[b] = bin_nbr_cnt[b] > 0 ? bin_nbr_sum[b] / bin_nbr_cnt[b] : 0.0;
    }

    kf.mean_error_mult        = n_err  > 0 ? sum_err_mult  / n_err  : 0.0;
    kf.mean_true_mult         = n_true > 0 ? sum_true_mult / n_true : 0.0;
    double mean_nbr           = n_nbr_err > 0 ? sum_nbr_err / n_nbr_err : 0.0;
    kf.neighbor_contrast_ratio = kf.mean_error_mult > 0.0
                                 ? mean_nbr / kf.mean_error_mult : 0.0;
    kf.isolated_error_frac    = n_err > 0 ? (double)n_isolated / n_err : 0.0;

    // ---- label: optimal min_count via F1 maximisation ----
    {
        uint32_t max_mult = 0;
        for (const Entry& e : entries) if (e.mult > max_mult) max_mult = e.mult;

        std::vector<uint64_t> err_at  (max_mult + 2, 0);
        std::vector<uint64_t> total_at(max_mult + 2, 0);
        for (const Entry& e : entries) {
            total_at[e.mult]++;
            if (e.is_error) err_at[e.mult]++;
        }

        // predict_error = (mult < tau) → scan tau from 2 upward
        uint64_t TP = 0, FP = 0, FN = n_err;
        double best_f1 = 0.0;
        int best_tau = 1;  // tau=1 means nothing is predicted error, F1=0
        for (uint32_t tau = 2; tau <= max_mult + 1; ++tau) {
            TP += err_at[tau - 1];
            FP += total_at[tau - 1] - err_at[tau - 1];
            FN -= err_at[tau - 1];
            double denom = 2.0 * TP + FP + FN;
            double f1    = denom > 0.0 ? 2.0 * TP / denom : 0.0;
            if (f1 > best_f1) { best_f1 = f1; best_tau = (int)tau; }
        }
        labels_out.min_count      = best_tau;
        labels_out.error_fraction = !entries.empty()
                                    ? (double)n_err / entries.size() : 0.0;
    }

    // ---- label: optimal contrast_threshold via F1 (AND condition) ----
    {
        const int tau1 = labels_out.min_count;

        struct LowEntry { double contrast; bool is_error; };
        std::vector<LowEntry> low;
        uint64_t TP = 0, FP = 0, FN = 0;

        for (const Entry& e : entries) {
            if ((int)e.mult < tau1) {
                double c = (e.mult > 0 && e.max_nbr_mult > 0)
                           ? (double)e.max_nbr_mult / e.mult : 0.0;
                low.push_back({c, e.is_error});
                if (e.is_error) ++TP; else ++FP;
            } else {
                if (e.is_error) ++FN;
            }
        }

        // Sort ascending by contrast so we can scan increasing tau2.
        std::sort(low.begin(), low.end(),
                  [](const LowEntry& a, const LowEntry& b)
                  { return a.contrast < b.contrast; });

        // At tau2=0: all low-mult k-mers predicted error (baseline).
        double best_f1   = 0.0;
        double best_tau2 = 0.0;
        {
            double denom = 2.0 * TP + FP + FN;
            double f1    = denom > 0.0 ? 2.0 * TP / denom : 0.0;
            if (f1 > best_f1) { best_f1 = f1; best_tau2 = 0.0; }
        }
        // Raise tau2: each step removes the lowest-contrast k-mer from predicted-error.
        for (const LowEntry& le : low) {
            if (le.is_error) { --TP; ++FN; } else { --FP; }
            double denom = 2.0 * TP + FP + FN;
            double f1    = denom > 0.0 ? 2.0 * TP / denom : 0.0;
            if (f1 > best_f1) { best_f1 = f1; best_tau2 = le.contrast; }
        }
        labels_out.contrast_threshold = best_tau2;
    }

    return kf;
}

// ---------------------------------------------------------------------------
// Top-level driver
// ---------------------------------------------------------------------------

TopoFeatures RunExtraction(SDBG& dbg, const ExtractionOptions& opts) {
    using Clock = std::chrono::steady_clock;
    auto t0 = Clock::now();

    TopoFeatures tf = {};
    tf.node = ExtractNodeFeatures(dbg, opts.total_reads);

    auto t1 = Clock::now();
    tf.timing.node_ms =
        std::chrono::duration<double, std::milli>(t1 - t0).count();

    if (!opts.ref_fasta.empty() && !opts.reads1.empty()) {
        std::cerr << "Extracting k-mer features via jellyfish...\n";
        tf.kmer            = ExtractKmerFeatures(opts, tf.labels);
        tf.has_kmer_features = true;
        auto t2 = Clock::now();
        tf.timing.kmer_ms  =
            std::chrono::duration<double, std::milli>(t2 - t1).count();
    }

    return tf;
}