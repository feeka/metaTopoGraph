#include "topo_extractor.h"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <vector>

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

// Walk both branches from a 2-out-degree edge with bounded depth.
// A bubble is detected when both walk fronts share the same target node,
// which is identified by Forward(cur0) == Forward(cur1).
static BubbleResult FindBubbleFromBranch(SDBG& dbg, uint64_t edge,
                                         uint64_t max_depth) {
    uint64_t outgoings[8];
    int od = dbg.OutgoingEdges(edge, outgoings);
    if (od != 2) return {false, 0.0, 0.0};

    uint64_t cur0 = outgoings[0];
    uint64_t cur1 = outgoings[1];

    // Use the multiplicity of the first edge on each branch as representative.
    const double mult0 = dbg.EdgeMultiplicity(cur0);
    const double mult1 = dbg.EdgeMultiplicity(cur1);

    for (uint64_t depth = 0; depth < max_depth; ++depth) {
        // Both branches now point at the same target node → bubble closes.
        if (dbg.Forward(cur0) == dbg.Forward(cur1)) {
            return {true, mult0, mult1};
        }
        uint64_t next0 = dbg.UniqueNextEdge(cur0);
        uint64_t next1 = dbg.UniqueNextEdge(cur1);
        if (next0 == SDBG::kNullID || next1 == SDBG::kNullID) break;
        cur0 = next0;
        cur1 = next1;
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

    for (uint64_t i = 0; i < dbg.size(); ++i) {
        if (dbg.IsValidEdge(i)) {
            mul_t m = dbg.EdgeMultiplicity(i);
            uint32_t idx = std::min(static_cast<uint32_t>(m), HIST_CAP);
            hist[idx]++;
            total++;
        }
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

    // Smooth histogram (sigma=1.5, radius=3) for valley/mode finding.
    std::vector<double> smoothed = GaussianSmooth(hist, 1.5, 3);

    // Feature 1 — valley_position: local minimum between mult 2 and 50
    uint32_t valley_pos = 2;
    double   valley_val = smoothed[2];
    for (uint32_t i = 3; i <= std::min(50u, HIST_CAP); ++i) {
        if (smoothed[i] < valley_val) {
            valley_val = smoothed[i];
            valley_pos = i;
        }
    }
    hf.valley_position = valley_pos;

    // Feature 2 — valley_depth = count_at_valley / count_at_error_peak
    {
        double error_peak = smoothed[1] > 0.0 ? smoothed[1] : 1.0;
        hf.valley_depth = valley_val / error_peak;
    }

    // Features 6 & 7 — n_signal_modes and primary_mode_depth
    // Count local maxima right of valley; track argmax.
    int      n_modes          = 0;
    uint32_t primary_mode     = valley_pos;
    double   primary_mode_val = 0.0;
    const uint32_t right_end  = HIST_CAP - 1;

    for (uint32_t i = valley_pos + 1; i < right_end; ++i) {
        if (smoothed[i] > smoothed[i - 1] && smoothed[i] > smoothed[i + 1]) {
            ++n_modes;
            if (smoothed[i] > primary_mode_val) {
                primary_mode_val = smoothed[i];
                primary_mode     = i;
            }
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

    uint64_t n_valid = 0, n_tips = 0, n_branch = 0, n_linear = 0, n_high = 0;
    double   tip_mult_sum = 0.0, branch_mult_sum = 0.0;

    for (uint64_t i = 0; i < n; i += stride) {
        if (!dbg.IsValidEdge(i)) continue;
        ++n_valid;

        const int    indeg  = dbg.EdgeIndegree(i);
        const int    outdeg = dbg.EdgeOutdegree(i);
        const double m      = dbg.EdgeMultiplicity(i);

        if (outdeg == 0 || indeg == 0) {
            ++n_tips;
            tip_mult_sum += m;
        }
        if (indeg >= 2 || outdeg >= 2) {
            ++n_branch;
            branch_mult_sum += m;
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

    uint64_t n_tips      = 0;
    uint64_t total_steps = 0;

    for (uint64_t i = 0; i < n; i += stride) {
        if (!dbg.IsValidEdge(i)) continue;
        if (dbg.EdgeOutdegree(i) != 0) continue; // only tail tips
        if (dbg.EdgeIndegree(i)  == 0) continue; // isolated edge — skip

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

    return n_tips > 0 ? static_cast<double>(total_steps) / n_tips : 0.0;
}

// ---------------------------------------------------------------------------
// Phase 4: bubble detection (features 16-18)
// ---------------------------------------------------------------------------

BubbleFeatures ExtractBubbleFeatures(SDBG& dbg,
                                     const ExtractionOptions& opts,
                                     double valley_position) {
    const uint64_t n      = dbg.size();
    const uint64_t stride = std::max(uint64_t(1), n / opts.sample_size);

    uint64_t n_branching  = 0;
    uint64_t n_bubbles    = 0;
    uint64_t n_error      = 0;
    uint64_t n_balanced   = 0;

    for (uint64_t i = 0; i < n; i += stride) {
        if (!dbg.IsValidEdge(i)) continue;
        if (dbg.EdgeOutdegree(i) != 2) continue;

        ++n_branching;
        BubbleResult res = FindBubbleFromBranch(dbg, i, opts.max_bubble_depth);
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

    BubbleFeatures bf = {};
    bf.bubble_density =
        n_branching > 0 ? static_cast<double>(n_bubbles) / n_branching : 0.0;
    bf.error_bubble_fraction =
        n_bubbles > 0 ? static_cast<double>(n_error)    / n_bubbles : 0.0;
    bf.balanced_bubble_fraction =
        n_bubbles > 0 ? static_cast<double>(n_balanced) / n_bubbles : 0.0;
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

    tf.walk.mean_tip_length = ExtractMeanTipLength(dbg, opts);
    auto t3 = Clock::now();

    {
        BubbleFeatures bf = ExtractBubbleFeatures(dbg, opts, tf.hist.valley_position);
        tf.walk.bubble_density           = bf.bubble_density;
        tf.walk.error_bubble_fraction    = bf.error_bubble_fraction;
        tf.walk.balanced_bubble_fraction = bf.balanced_bubble_fraction;
    }
    auto t4 = Clock::now();

    tf.timing.histogram_ms = ms(t0, t1);
    tf.timing.node_ms      = ms(t1, t2);
    tf.timing.walk_ms      = ms(t2, t3);
    tf.timing.bubble_ms    = ms(t3, t4);

    return tf;
}
