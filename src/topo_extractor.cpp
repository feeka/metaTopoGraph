#include "topo_extractor.h"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <limits>
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

    const double C_real    = n_branch > 0 ? max_r_sum / n_branch : 0.0;
    const double mean_mult = mult_sum / total;

    NodeFeatures nf = {};
    nf.mult_min              = (mult_min ==  std::numeric_limits<double>::max()) ? 0.0 : mult_min;
    nf.mult_max              = (mult_max == -std::numeric_limits<double>::max()) ? 0.0 : mult_max;
    nf.mult_mean             = mean_mult;
    nf.n_reads               = total_reads;
    nf.n_kmers               = total;
    nf.n_tips                = n_tips;
    nf.n_branch_nodes        = n_branch;
    nf.mean_max_branch_ratio = C_real;
    nf.max_branch_ratio      = (n_branch > 0 && global_max_r   != -std::numeric_limits<double>::max()) ? global_max_r   : 0.0;
    nf.mean_jump_ratio       = (n_branch > 0 && global_max_mag != -std::numeric_limits<double>::max()) ? jump_r_sum / n_branch : 0.0;
    return nf;
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

    return tf;
}