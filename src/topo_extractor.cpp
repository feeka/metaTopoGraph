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
// Full-scan node feature extraction
// ---------------------------------------------------------------------------

NodeFeatures ExtractNodeFeatures(SDBG& dbg) {
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
        uint64_t n_jumps       = 0;
        uint64_t n_tips        = 0;
        uint64_t n_branch      = 0;
        double   min_ratio_sum = 0.0;
        double   max_ratio_sum = 0.0;
    };
    std::vector<ThrAcc> acc(nthreads);

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

        // Global stats
        a.mult_sum += mul;
        ++a.total;
        if (mul < a.mult_min) a.mult_min = mul;
        if (mul > a.mult_max) a.mult_max = mul;

        // Tip: dead-end on either side
        if (indeg == 0 || outdeg == 0) ++a.n_tips;

        // Outgoing neighbor analysis
        if (outdeg > 0) {
            uint64_t outs[8];
            int od = dbg.OutgoingEdges(i, outs);

            // Prominent jumps: every directed edge (i â†’ out_j) where ratio >= JUMP_THRESHOLD
            for (int j = 0; j < od; ++j) {
                const double mul_j = dbg.EdgeMultiplicity(outs[j]);
                const double hi    = std::max(mul, mul_j);
                const double lo    = std::min(mul, mul_j);
                const double ratio = lo > 0.0 ? hi / lo : hi;
                if (ratio >= JUMP_THRESHOLD) ++a.n_jumps;
            }

            // Branch proportion ratios: only for out-degree >= 2
            if (od >= 2) {
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
    }

    // Merge per-thread accumulators
    double   mult_min  =  std::numeric_limits<double>::max();
    double   mult_max  = -std::numeric_limits<double>::max();
    double   mult_sum  = 0.0;
    uint64_t total     = 0;
    uint64_t n_jumps   = 0;
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
        n_jumps   += a.n_jumps;
        n_tips    += a.n_tips;
        n_branch  += a.n_branch;
        min_r_sum += a.min_ratio_sum;
        max_r_sum += a.max_ratio_sum;
    }

    NodeFeatures nf = {};
    if (total == 0) return nf;

    nf.mult_min  = (mult_min ==  std::numeric_limits<double>::max()) ? 0.0 : mult_min;
    nf.mult_max  = (mult_max == -std::numeric_limits<double>::max()) ? 0.0 : mult_max;
    nf.mult_mean = mult_sum / total;

    nf.n_prominent_jumps     = n_jumps;
    nf.n_tips                = n_tips;
    nf.mean_min_branch_ratio = n_branch > 0 ? min_r_sum / n_branch : 0.0;
    nf.mean_max_branch_ratio = n_branch > 0 ? max_r_sum / n_branch : 0.0;

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

