#include "topo_extractor.h"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <limits>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "sdbg/sdbg_def.h"

// ---------------------------------------------------------------------------
// Compute Bayesian error threshold m* from graph statistics.
// Returns the multiplicity at which real and error Poisson posteriors are equal.
// Edges with mult <= m* are more likely erroneous than real.
// Returns 0.0 when coverage is too low to determine a meaningful threshold.
// ---------------------------------------------------------------------------

static double ComputeErrorThreshold(double C_real, double mult_mean, double kmer_k) {
    // Requires: C_real > 1, 1 < mult_mean < C_real
    if (C_real <= 1.0 || mult_mean <= 1.0 || mult_mean >= C_real)
        return 0.0;

    // E_err/E_real from the mixture mean equation:
    //   mult_mean = (C_real * E_real + 1 * E_err) / (E_real + E_err)
    const double ratio = (C_real - mult_mean) / (mult_mean - 1.0);
    if (ratio <= 0.0) return 0.0;

    // p = ratio / (k * C_real),  lambda_e = C_real * p / 3
    const double lambda_e = ratio / (3.0 * kmer_k);
    if (lambda_e <= 0.0 || lambda_e >= C_real) return 0.0;

    // m* = (ln(ratio) + (C_real - lambda_e)) / ln(C_real / lambda_e)
    const double m_star = (std::log(ratio) + (C_real - lambda_e))
                        / std::log(C_real / lambda_e);
    return m_star < 0.0 ? 0.0 : m_star;
}

// ---------------------------------------------------------------------------
// Full-scan node feature extraction (two passes)
// ---------------------------------------------------------------------------

NodeFeatures ExtractNodeFeatures(SDBG& dbg) {
    const double kmer_k = static_cast<double>(dbg.k());

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
        double   min_ratio_sum = 0.0;
        double   max_ratio_sum = 0.0;
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
            double min_r =  std::numeric_limits<double>::max();
            double max_r = -std::numeric_limits<double>::max();
            for (int j = 0; j < od; ++j) {
                const double mul_j = dbg.EdgeMultiplicity(outs[j]);
                const double r     = mul_j > 0.0 ? mul / mul_j : mul;
                if (r < min_r) min_r = r;
                if (r > max_r) max_r = r;
            }
            a.min_ratio_sum += min_r;
            a.max_ratio_sum += max_r;
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
    double   min_r_sum = 0.0;
    double   max_r_sum = 0.0;

    for (const ThrAcc& a : acc) {
        if (a.total == 0) continue;
        if (a.mult_min < mult_min) mult_min = a.mult_min;
        if (a.mult_max > mult_max) mult_max = a.mult_max;
        mult_sum  += a.mult_sum;
        total     += a.total;
        n_tips    += a.n_tips;
        n_branch  += a.n_branch;
        min_r_sum += a.min_ratio_sum;
        max_r_sum += a.max_ratio_sum;
    }

    if (total == 0) return NodeFeatures{};

    const double C_real    = n_branch > 0 ? max_r_sum / n_branch : 0.0;
    const double mean_mult = mult_sum / total;
    const double m_star    = ComputeErrorThreshold(C_real, mean_mult, kmer_k);

    // ---- Pass 2: count prominent jumps using adaptive threshold m* ----
    uint64_t n_jumps = 0;
    if (m_star > 0.0) {
        std::vector<uint64_t> jump_acc(nthreads, 0);

#pragma omp parallel for schedule(dynamic, 512)
        for (uint64_t i = 0; i < dbg.size(); ++i) {
            if (!dbg.IsValidEdge(i)) continue;
            const int outdeg = dbg.EdgeOutdegree(i);
            if (outdeg == 0) continue;
#ifdef _OPENMP
            const int tid = omp_get_thread_num();
#else
            const int tid = 0;
#endif
            const double mul = dbg.EdgeMultiplicity(i);
            uint64_t outs[8];
            int od = dbg.OutgoingEdges(i, outs);
            for (int j = 0; j < od; ++j) {
                const double mul_j = dbg.EdgeMultiplicity(outs[j]);
                if (std::min(mul, mul_j) <= m_star)
                    ++jump_acc[tid];
            }
        }
        for (uint64_t v : jump_acc) n_jumps += v;
    }

    NodeFeatures nf = {};
    nf.mult_min              = (mult_min ==  std::numeric_limits<double>::max()) ? 0.0 : mult_min;
    nf.mult_max              = (mult_max == -std::numeric_limits<double>::max()) ? 0.0 : mult_max;
    nf.mult_mean             = mean_mult;
    nf.n_prominent_jumps     = n_jumps;
    nf.n_tips                = n_tips;
    nf.mean_min_branch_ratio = n_branch > 0 ? min_r_sum / n_branch : 0.0;
    nf.mean_max_branch_ratio = C_real;
    nf.error_threshold       = m_star;
    return nf;
}

// ---------------------------------------------------------------------------
// Top-level driver
// ---------------------------------------------------------------------------

TopoFeatures RunExtraction(SDBG& dbg, const ExtractionOptions& /*opts*/) {
    using Clock = std::chrono::steady_clock;
    auto t0 = Clock::now();

    TopoFeatures tf = {};
    tf.node = ExtractNodeFeatures(dbg);

    auto t1 = Clock::now();
    tf.timing.node_ms =
        std::chrono::duration<double, std::milli>(t1 - t0).count();

    return tf;
}