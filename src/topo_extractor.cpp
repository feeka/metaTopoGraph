#include "topo_extractor.h"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "sdbg/sdbg_def.h"

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

// Convolve hist[0..size-1] with a discrete Gaussian kernel of the given sigma
// and radius. Returns a smoothed vector of the same length.
static std::vector<double> GaussianSmooth(const std::vector<uint64_t>& hist,
                                          double sigma, int radius) {
    const int n = static_cast<int>(hist.size());
    std::vector<double> kernel(2 * radius + 1);
    double ksum = 0.0;
    for (int i = -radius; i <= radius; ++i) {
        kernel[i + radius] = std::exp(-i * i / (2.0 * sigma * sigma));
        ksum += kernel[i + radius];
    }
    for (double& w : kernel) w /= ksum;

    std::vector<double> out(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double s = 0.0;
        for (int k = -radius; k <= radius; ++k) {
            int idx = i + k;
            if (idx >= 0 && idx < n) {
                s += kernel[k + radius] * static_cast<double>(hist[idx]);
            }
        }
        out[i] = s;
    }
    return out;
}

// Result of one bubble search from a single branching edge.
struct BubbleResult {
    bool   found;
    double mult0;   // representative multiplicity of branch 0
    double mult1;   // representative multiplicity of branch 1
};

// Walk two branches from a branching edge's two outgoing edges.
// Convergence: Forward(cur0) == Forward(cur1) means both edges target the
// same node (Forward returns the unique last-edge-of-node identifier).
// Aborts if either cursor reaches a node with out-degree != 1: a real bubble
// has strictly linear paths between source and sink, so any mid-walk branch
// means we are not looking at a simple bubble.
static BubbleResult WalkBranchPair(SDBG& dbg,
                                   uint64_t branch0, uint64_t branch1,
                                   uint64_t max_depth) {
    const double mult0 = dbg.EdgeMultiplicity(branch0);
    const double mult1 = dbg.EdgeMultiplicity(branch1);

    uint64_t cur0 = branch0;
    uint64_t cur1 = branch1;

    for (uint64_t depth = 0; depth < max_depth; ++depth) {
        // Forward(e) = last edge of the target node of e — unique per node.
        if (dbg.Forward(cur0) == dbg.Forward(cur1))
            return {true, mult0, mult1};

        uint64_t outs0[8], outs1[8];
        int od0 = dbg.OutgoingEdges(cur0, outs0);
        int od1 = dbg.OutgoingEdges(cur1, outs1);

        // Abort if either path hits a branching node or a dead-end mid-walk.
        if (od0 != 1 || od1 != 1) break;

        cur0 = outs0[0];
        cur1 = outs1[0];
    }
    return {false, 0.0, 0.0};
}

// ---------------------------------------------------------------------------
// Phase 1: histogram features
// ---------------------------------------------------------------------------

HistogramFeatures ExtractHistogramFeatures(SDBG& dbg) {
    // Build raw multiplicity histogram up to HIST_CAP.
    const uint32_t HIST_CAP = 1000;
    std::vector<uint64_t> hist(HIST_CAP + 1, 0);
    uint64_t total = 0;

    // Phase 1 scans every edge — parallelise with per-thread private histograms
    // to avoid false sharing on the shared array.
#ifdef _OPENMP
    const int nthreads = omp_get_max_threads();
#else
    const int nthreads = 1;
#endif
    std::vector<std::vector<uint64_t>> phist(nthreads,
                                             std::vector<uint64_t>(HIST_CAP + 1, 0));
    std::vector<uint64_t> ptotal(nthreads, 0);

#pragma omp parallel for schedule(static)
    for (uint64_t i = 0; i < dbg.size(); ++i) {
        if (dbg.IsValidEdge(i)) {
#ifdef _OPENMP
            const int tid = omp_get_thread_num();
#else
            const int tid = 0;
#endif
            mul_t m = dbg.EdgeMultiplicity(i);
            uint32_t idx = std::min(static_cast<uint32_t>(m), HIST_CAP);
            phist[tid][idx]++;
            ptotal[tid]++;
        }
    }

    for (int t = 0; t < nthreads; ++t) {
        for (uint32_t j = 0; j <= HIST_CAP; ++j) hist[j] += phist[t][j];
        total += ptotal[t];
    }

    HistogramFeatures hf = {};
    if (total == 0) return hf;

    // Feature 3 — mult_1_fraction
    hf.mult_1_fraction = static_cast<double>(hist[1]) / total;

    // Features 4 & 5 — mean and CV
    double sum = 0.0, sum_sq = 0.0;
    for (uint32_t i = 1; i <= HIST_CAP; ++i) {
        double fi = static_cast<double>(i);
        sum    += fi * hist[i];
        sum_sq += fi * fi * hist[i];
    }
    hf.mean_node_multiplicity = sum / total;
    {
        double variance = sum_sq / total - hf.mean_node_multiplicity * hf.mean_node_multiplicity;
        hf.multiplicity_cv = hf.mean_node_multiplicity > 0.0
            ? std::sqrt(std::max(0.0, variance)) / hf.mean_node_multiplicity
            : 0.0;
    }

    // Use wider smoothing (sigma=5) for valley and mode detection so that noise
    // ripples don't produce spurious peaks.  The tighter sigma=1.5 pass was
    // already used for mean/CV above.
    std::vector<double> smoothed = GaussianSmooth(hist, 5.0, 15);

    // Feature 1 — valley_position: local minimum between mult 2 and the rough
    // signal peak.  We first find the rough peak (argmax right of mult=10) to
    // set a dynamic upper bound, avoiding the fixed-50 boundary artefact.
    uint32_t rough_peak     = 10;
    double   rough_peak_val = smoothed[10];
    for (uint32_t i = 11; i <= HIST_CAP; ++i) {
        if (smoothed[i] > rough_peak_val) {
            rough_peak_val = smoothed[i];
            rough_peak     = i;
        }
    }
    // Search for valley between mult=2 and halfway to the rough peak.
    const uint32_t valley_upper = std::max(50u, rough_peak / 2);
    uint32_t valley_pos = 2;
    double   valley_val = smoothed[2];
    for (uint32_t i = 3; i <= std::min(valley_upper, HIST_CAP); ++i) {
        if (smoothed[i] < valley_val) {
            valley_val = smoothed[i];
            valley_pos = i;
        }
    }
    hf.valley_position = valley_pos;

    // Feature 2 — valley_depth = count_at_valley / max(smoothed[1..valley_pos]).
    // Using the actual smoothed peak in [1, valley_pos) rather than just smoothed[1]
    // because the error-read peak may not fall exactly at multiplicity 1.
    {
        double error_peak = smoothed[1];
        for (uint32_t i = 2; i < valley_pos; ++i)
            if (smoothed[i] > error_peak) error_peak = smoothed[i];
        if (error_peak <= 0.0) error_peak = 1.0;
        hf.valley_depth = valley_val / error_peak;
    }

    // Features 6 & 7 — n_signal_modes and primary_mode_depth
    // Find the global maximum right of valley first (the primary mode).
    // Only count secondary peaks that exceed mode_prominence_frac of that max.
    // This prevents counting noise bumps as distinct modes.
    const double mode_prominence_frac = 0.10;
    uint32_t primary_mode     = valley_pos;
    double   primary_mode_val = 0.0;
    const uint32_t right_end  = HIST_CAP - 1;

    for (uint32_t i = valley_pos + 1; i <= right_end; ++i) {
        if (smoothed[i] > primary_mode_val) {
            primary_mode_val = smoothed[i];
            primary_mode     = i;
        }
    }

    int n_modes = 0;
    const double prominence_threshold = mode_prominence_frac * primary_mode_val;
    for (uint32_t i = valley_pos + 1; i < right_end; ++i) {
        if (smoothed[i] > smoothed[i - 1] &&
            smoothed[i] > smoothed[i + 1] &&
            smoothed[i] >= prominence_threshold) {
            ++n_modes;
        }
    }
    hf.n_signal_modes    = n_modes;
    hf.primary_mode_depth = primary_mode;

    // Feature 8 — high_mult_tail_ratio = sum(hist[3*primary:]) / total
    {
        uint64_t tail = 0;
        uint32_t tail_start = std::min(3u * primary_mode, HIST_CAP);
        for (uint32_t i = tail_start; i <= HIST_CAP; ++i) tail += hist[i];
        hf.high_mult_tail_ratio = static_cast<double>(tail) / total;
    }

    return hf;
}

// ---------------------------------------------------------------------------
// Phase 2: node-level features
// ---------------------------------------------------------------------------

NodeFeatures ExtractNodeFeatures(SDBG& dbg, const ExtractionOptions& opts) {
    const uint64_t n      = dbg.size();
    const uint64_t stride = std::max(uint64_t(1), n / opts.sample_size);

    // Build the strided index list so OpenMP can divide it evenly.
    std::vector<uint64_t> indices;
    indices.reserve(opts.sample_size + 1);
    for (uint64_t i = 0; i < n; i += stride) indices.push_back(i);
    const int64_t m = static_cast<int64_t>(indices.size());

    uint64_t n_valid = 0, n_tips = 0, n_branch = 0, n_linear = 0, n_high = 0;
    double   tip_mult_sum = 0.0, branch_mult_sum = 0.0;

#pragma omp parallel for schedule(dynamic,256) \
    reduction(+:n_valid,n_tips,n_branch,n_linear,n_high,tip_mult_sum,branch_mult_sum)
    for (int64_t j = 0; j < m; ++j) {
        uint64_t i = indices[j];
        if (!dbg.IsValidEdge(i)) continue;
        ++n_valid;

        const int    indeg  = dbg.EdgeIndegree(i);
        const int    outdeg = dbg.EdgeOutdegree(i);
        const double mul    = dbg.EdgeMultiplicity(i);

        if (outdeg == 0 || indeg == 0) {
            ++n_tips;
            tip_mult_sum += mul;
        }
        if (indeg >= 2 || outdeg >= 2) {
            ++n_branch;
            branch_mult_sum += mul;
        }
        if (indeg == 1 && outdeg == 1) ++n_linear;
        if ((indeg + outdeg) >= 4)     ++n_high;
    }

    NodeFeatures nf = {};
    if (n_valid == 0) return nf;

    nf.tip_density               = static_cast<double>(n_tips)   / n_valid;
    nf.branching_node_fraction   = static_cast<double>(n_branch) / n_valid;
    nf.linear_node_fraction      = static_cast<double>(n_linear) / n_valid;
    nf.high_degree_node_fraction = static_cast<double>(n_high)   / n_valid;
    nf.mult_at_tips     = n_tips   > 0 ? tip_mult_sum    / n_tips   : 0.0;
    nf.mult_at_branches = n_branch > 0 ? branch_mult_sum / n_branch : 0.0;
    return nf;
}

// ---------------------------------------------------------------------------
// Phase 3: tip walk (feature 15)
// ---------------------------------------------------------------------------

double ExtractMeanTipLength(SDBG& dbg, const ExtractionOptions& opts) {
    // Walk from sampled tip edges (outdeg == 0) backward along the path
    // until we reach a branch, another dead-end, or max_tip_length steps.
    const uint64_t n      = dbg.size();
    const uint64_t stride = std::max(uint64_t(1), n / opts.sample_size);

    std::vector<uint64_t> indices;
    indices.reserve(opts.sample_size + 1);
    for (uint64_t i = 0; i < n; i += stride) indices.push_back(i);
    const int64_t m = static_cast<int64_t>(indices.size());

    uint64_t n_tips      = 0;
    uint64_t total_steps = 0;

#pragma omp parallel for schedule(dynamic,128) reduction(+:n_tips,total_steps)
    for (int64_t j = 0; j < m; ++j) {
        uint64_t i = indices[j];
        if (!dbg.IsValidEdge(i)) continue;
        if (dbg.EdgeOutdegree(i) != 0) continue;
        if (dbg.EdgeIndegree(i)  == 0) continue;

        ++n_tips;
        uint64_t cur    = i;
        uint64_t length = 0;

        for (uint64_t step = 0; step < opts.max_tip_length; ++step) {
            uint64_t prev = dbg.UniquePrevEdge(cur);
            if (prev == SDBG::kNullID) break;
            cur = prev;
            ++length;
        }
        total_steps += length;
    }

    return {n_tips > 0 ? static_cast<double>(total_steps) / n_tips : 0.0,
            n_tips};
}

// ---------------------------------------------------------------------------
// Phase 4: bubble detection (features 16-18)
// ---------------------------------------------------------------------------

BubbleFeatures ExtractBubbleFeatures(SDBG& dbg,
                                     const ExtractionOptions& opts,
                                     double valley_position) {
    const uint64_t n      = dbg.size();
    const uint64_t stride = std::max(uint64_t(1), n / opts.sample_size);

    std::vector<uint64_t> indices;
    indices.reserve(opts.sample_size + 1);
    for (uint64_t i = 0; i < n; i += stride) indices.push_back(i);
    const int64_t m = static_cast<int64_t>(indices.size());

    uint64_t n_branching  = 0;
    uint64_t n_bubbles    = 0;
    uint64_t n_error      = 0;
    uint64_t n_balanced   = 0;

    // Each bubble walk is independent — dynamic schedule handles variable cost.
    // Check all branch pairs (not just od==2): handles degree-3+ branching nodes.
#pragma omp parallel for schedule(dynamic,64) \
    reduction(+:n_branching,n_bubbles,n_error,n_balanced)
    for (int64_t j = 0; j < m; ++j) {
        uint64_t i = indices[j];
        if (!dbg.IsValidEdge(i)) continue;

        uint64_t outgoings[8];
        int od = dbg.OutgoingEdges(i, outgoings);
        if (od < 2) continue;

        ++n_branching;

        // Check all pairs of outgoing branches.
        for (int b0 = 0; b0 < od; ++b0) {
            for (int b1 = b0 + 1; b1 < od; ++b1) {
                BubbleResult res = WalkBranchPair(dbg,
                                                  outgoings[b0], outgoings[b1],
                                                  opts.max_bubble_depth);
                if (!res.found) continue;

                ++n_bubbles;
                double hi    = std::max(res.mult0, res.mult1);
                double lo    = std::min(res.mult0, res.mult1);
                double ratio = lo > 0.0 ? hi / lo : hi;

                if (ratio > 5.0) {
                    ++n_error;
                } else if (ratio < 2.0 &&
                           res.mult0 > valley_position &&
                           res.mult1 > valley_position) {
                    ++n_balanced;
                }
            }
        }
    }

    BubbleFeatures bf = {};
    bf.bubble_density =
        n_branching > 0 ? static_cast<double>(n_bubbles) / n_branching : 0.0;
    bf.error_bubble_fraction =
        n_bubbles > 0 ? static_cast<double>(n_error)    / n_bubbles : 0.0;
    bf.balanced_bubble_fraction =
        n_bubbles > 0 ? static_cast<double>(n_balanced) / n_bubbles : 0.0;
    bf.n_branching_sampled = n_branching;
    bf.n_bubbles_found     = n_bubbles;
    return bf;
}

// ---------------------------------------------------------------------------
// Top-level driver
// ---------------------------------------------------------------------------

TopoFeatures RunExtraction(SDBG& dbg, const ExtractionOptions& opts) {
    using Clock = std::chrono::steady_clock;
    auto ms = [](Clock::time_point a, Clock::time_point b) {
        return std::chrono::duration<double, std::milli>(b - a).count();
    };

    TopoFeatures tf = {};

    auto t0 = Clock::now();
    tf.hist = ExtractHistogramFeatures(dbg);
    auto t1 = Clock::now();

    tf.node = ExtractNodeFeatures(dbg, opts);
    auto t2 = Clock::now();

    {
        TipWalkResult tw = ExtractMeanTipLength(dbg, opts);
        tf.walk.mean_tip_length = tw.mean_tip_length;
        tf.walk.n_tips_walked   = tw.n_tips_walked;
    }
    auto t3 = Clock::now();

    {
        BubbleFeatures bf = ExtractBubbleFeatures(dbg, opts, tf.hist.valley_position);
        tf.walk.bubble_density           = bf.bubble_density;
        tf.walk.error_bubble_fraction    = bf.error_bubble_fraction;
        tf.walk.balanced_bubble_fraction = bf.balanced_bubble_fraction;
        tf.walk.n_branching_sampled      = bf.n_branching_sampled;
        tf.walk.n_bubbles_found          = bf.n_bubbles_found;
    }
    auto t4 = Clock::now();

    tf.timing.histogram_ms = ms(t0, t1);
    tf.timing.node_ms      = ms(t1, t2);
    tf.timing.walk_ms      = ms(t2, t3);
    tf.timing.bubble_ms    = ms(t3, t4);

    return tf;
}
