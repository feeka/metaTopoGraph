"""
plot_features.py  —  visualise metaTopoGraph features.json output

Single-sample mode (one file):
    python plot_features.py features.json

Multi-sample mode (directory of *.json files):
    python plot_features.py runs/

Outputs (single sample):
  <name>_dashboard.png   — coverage model + node topology + bubble breakdown + radar
Multi-sample adds:
  features_pca.png       — PCA biplot with loadings
  features_correlation.png  — Pearson correlation matrix
  features_dendrogram.png   — hierarchical clustering
  features_scatter_matrix.png  — pairwise scatter of the 6 most variable features
"""

import json
import sys
import os
import glob
import math

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.stats import pearsonr


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def load_one(path: str) -> dict:
    with open(path) as f:
        d = json.load(f)
    flat = {}
    for group in ("histogram", "node", "walk"):
        for k, v in d[group].items():
            flat[k] = v
    flat["_name"] = os.path.splitext(os.path.basename(path))[0]
    flat["_timing"] = d["timing"]
    return flat


def load_dir(directory: str) -> list[dict]:
    paths = sorted(glob.glob(os.path.join(directory, "*.json")))
    if not paths:
        sys.exit(f"No *.json files found in {directory}")
    return [load_one(p) for p in paths]


FEATURE_NAMES = [
    "valley_position", "valley_depth", "mult_1_fraction",
    "mean_node_multiplicity", "multiplicity_cv", "n_signal_modes",
    "primary_mode_depth", "high_mult_tail_ratio",
    "tip_density", "branching_node_fraction", "linear_node_fraction",
    "high_degree_node_fraction", "mult_at_tips", "mult_at_branches",
    "mean_tip_length", "bubble_density",
    "error_bubble_fraction", "balanced_bubble_fraction",
]

GROUPS = {
    "Histogram":   FEATURE_NAMES[:8],
    "Node":        FEATURE_NAMES[8:14],
    "Walk/Bubble": FEATURE_NAMES[14:],
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def gaussian(x, mu, sigma, amp):
    return amp * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def reconstruct_coverage_model(s: dict):
    """
    Reconstruct a plausible coverage histogram curve from the 8 histogram
    features using a two-component Gaussian mixture:
      - error component: mean~1, amplitude inferred from mult_1_fraction
      - signal component: mean=primary_mode_depth, sigma=mean*cv
    Returns (x, y_error, y_signal, valley_x).
    """
    cov     = s["primary_mode_depth"]
    mean    = s["mean_node_multiplicity"]
    cv      = s["multiplicity_cv"]
    valley  = s["valley_position"]
    m1_frac = s["mult_1_fraction"]
    tail    = s["high_mult_tail_ratio"]

    x_max = int(max(cov * 3, valley * 4, 20))
    x = np.linspace(1, x_max, 500)

    sigma_signal = max(mean * cv, 1.0)
    amp_signal   = 1.0
    amp_error    = m1_frac * 10 if m1_frac > 0 else 0.0

    y_signal = gaussian(x, cov,    sigma_signal, amp_signal)
    y_error  = gaussian(x, 1.0,    0.5,          amp_error)
    y_tail   = tail * 0.3 * np.exp(-0.01 * (x - cov * 3))  # rough tail
    y_tail   = np.clip(y_tail, 0, None)

    return x, y_error, y_signal, y_tail, valley


def radar_chart(ax, categories, values, title, color="steelblue"):
    """Draw a radar / spider chart on a polar axis."""
    n = len(categories)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False).tolist()
    values_plot = values + [values[0]]
    angles      = angles  + [angles[0]]

    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.plot(angles, values_plot, color=color, linewidth=2)
    ax.fill(angles, values_plot, color=color, alpha=0.25)
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories, size=7)
    ax.set_title(title, size=9, pad=12)
    ax.set_ylim(0, 1)


def normalise_for_radar(s: dict, feature_list: list, ref_max: dict | None = None) -> list:
    """
    Min-max normalise features to [0,1].
    ref_max: optional dict of {feature: max_val} for cross-sample normalisation.
    """
    vals = [s.get(f, 0.0) for f in feature_list]
    out  = []
    for f, v in zip(feature_list, vals):
        mx = ref_max[f] if ref_max and f in ref_max else max(abs(v), 1e-9)
        out.append(min(abs(v) / mx, 1.0))
    return out


# ---------------------------------------------------------------------------
# Single-sample dashboard
# ---------------------------------------------------------------------------

def plot_single(sample: dict):
    name = sample["_name"]
    fig  = plt.figure(figsize=(18, 10))
    fig.suptitle(f"Graph topology dashboard — {name}", fontsize=14, y=1.01)
    gs   = gridspec.GridSpec(2, 4, figure=fig, wspace=0.4, hspace=0.5)

    # 1. Reconstructed coverage model
    ax1 = fig.add_subplot(gs[0, :2])
    x, y_err, y_sig, y_tail, valley = reconstruct_coverage_model(sample)
    ax1.fill_between(x, y_err,              alpha=0.5, color="#e15759", label="Error component")
    ax1.fill_between(x, y_sig + y_tail,     alpha=0.4, color="#4e79a7", label="Signal component")
    ax1.axvline(valley, color="black",  ls="--", lw=1.2, label=f"Valley @ {valley:.0f}")
    ax1.axvline(sample["primary_mode_depth"], color="green", ls=":", lw=1.5,
                label=f"Primary mode @ {sample['primary_mode_depth']:.0f}×")
    ax1.set_xlabel("Multiplicity")
    ax1.set_ylabel("Relative density")
    ax1.set_title("Reconstructed multiplicity distribution")
    ax1.legend(fontsize=8)
    ax1.set_xlim(left=0)

    # 2. Node topology composition — stacked horizontal bar
    ax2 = fig.add_subplot(gs[0, 2])
    tip    = sample["tip_density"]
    branch = sample["branching_node_fraction"]
    linear = sample["linear_node_fraction"]
    high   = sample["high_degree_node_fraction"]
    other  = max(0.0, 1.0 - tip - branch - linear - high)
    cats   = ["Linear", "Branching", "Tip", "High-degree", "Other"]
    vals   = [linear,    branch,      tip,   high,           other]
    colors = ["#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#bab0ac"]
    bottom = 0.0
    for cat, val, col in zip(cats, vals, colors):
        if val > 0:
            ax2.bar(0, val, bottom=bottom, color=col, width=0.5, label=f"{cat} {val:.2%}")
            if val > 0.02:
                ax2.text(0, bottom + val / 2, f"{val:.1%}", ha="center", va="center",
                         fontsize=8, color="white", fontweight="bold")
            bottom += val
    ax2.set_xlim(-0.5, 0.5)
    ax2.set_ylim(0, 1)
    ax2.set_xticks([])
    ax2.set_title("Node topology\ncomposition")
    ax2.legend(loc="upper left", bbox_to_anchor=(1, 1), fontsize=7)

    # 3. Bubble breakdown — pie
    ax3 = fig.add_subplot(gs[0, 3])
    bd    = sample["bubble_density"]
    err_f = sample["error_bubble_fraction"]
    bal_f = sample["balanced_bubble_fraction"]
    other_f = max(0.0, 1.0 - err_f - bal_f)
    no_bubble = max(0.0, 1.0 - bd) if bd <= 1.0 else 0.0
    if bd > 0:
        pie_vals = [err_f * bd, bal_f * bd, other_f * bd, no_bubble]
        pie_lbls = ["Error\nbubbles", "Balanced\nbubbles", "Other\nbubbles", "No bubble"]
        pie_cols = ["#e15759", "#59a14f", "#f28e2b", "#d3d3d3"]
        pie_vals_nz = [(v, l, c) for v, l, c in zip(pie_vals, pie_lbls, pie_cols) if v > 0]
        pv, pl, pc = zip(*pie_vals_nz)
        ax3.pie(pv, labels=pl, colors=pc, autopct="%1.1f%%", textprops={"fontsize": 8})
    else:
        ax3.pie([1], labels=["No bubbles\ndetected"], colors=["#d3d3d3"])
    ax3.set_title(f"Bubble composition\n(density={bd:.4f})")

    # 4. Radar: histogram features (bottom-left)
    ax4 = fig.add_subplot(gs[1, :2], polar=True)
    hist_feats  = GROUPS["Histogram"]
    hist_labels = [f.replace("_", "\n") for f in hist_feats]
    hist_vals   = normalise_for_radar(sample, hist_feats)
    radar_chart(ax4, hist_labels, hist_vals, "Histogram feature profile", color="#4e79a7")

    # 5. Radar: node + walk features (bottom-right)
    ax5 = fig.add_subplot(gs[1, 2:], polar=True)
    nw_feats  = GROUPS["Node"] + GROUPS["Walk/Bubble"]
    nw_labels = [f.replace("_", "\n") for f in nw_feats]
    nw_vals   = normalise_for_radar(sample, nw_feats)
    radar_chart(ax5, nw_labels, nw_vals, "Node/Walk feature profile", color="#f28e2b")

    out = f"{name}_dashboard.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")
    plt.show()


# ---------------------------------------------------------------------------
# Multi-sample plots
# ---------------------------------------------------------------------------

def _build_matrix(samples):
    return np.array([[s.get(f, 0.0) for f in FEATURE_NAMES] for s in samples],
                    dtype=float)


def plot_pca(samples, matrix):
    names = [s["_name"] for s in samples]
    mu    = matrix.mean(axis=0)
    std   = matrix.std(axis=0)
    std[std == 0] = 1.0
    Z = (matrix - mu) / std

    cov = np.cov(Z.T)
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]
    explained = eigvals / eigvals.sum() * 100

    scores = Z @ eigvecs[:, :2]

    fig, ax = plt.subplots(figsize=(9, 7))
    scatter = ax.scatter(scores[:, 0], scores[:, 1],
                         c=range(len(names)), cmap="tab20", s=80, zorder=3)
    for i, n in enumerate(names):
        ax.annotate(n, scores[i], fontsize=7, xytext=(4, 4),
                    textcoords="offset points")

    # Loading arrows (top 6 by magnitude)
    loadings = eigvecs[:, :2]
    magnitudes = np.linalg.norm(loadings, axis=1)
    top_idx = np.argsort(magnitudes)[-6:]
    scale   = scores.ptp(axis=0).max() * 0.45
    for i in top_idx:
        ax.annotate("", xy=(loadings[i, 0] * scale, loadings[i, 1] * scale),
                    xytext=(0, 0),
                    arrowprops=dict(arrowstyle="->", color="#e15759", lw=1.5))
        ax.text(loadings[i, 0] * scale * 1.1, loadings[i, 1] * scale * 1.1,
                FEATURE_NAMES[i], fontsize=7, color="#e15759")

    ax.axhline(0, color="grey", lw=0.5, ls="--")
    ax.axvline(0, color="grey", lw=0.5, ls="--")
    ax.set_xlabel(f"PC1 ({explained[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({explained[1]:.1f}%)")
    ax.set_title("PCA of graph topology features (with loading arrows)")
    plt.tight_layout()
    plt.savefig("features_pca.png", dpi=150)
    print("Saved features_pca.png")
    plt.show()


def plot_correlation(matrix):
    n = len(FEATURE_NAMES)
    corr = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if matrix[:, i].std() == 0 or matrix[:, j].std() == 0:
                corr[i, j] = 0.0
            else:
                corr[i, j], _ = pearsonr(matrix[:, i], matrix[:, j])

    fig, ax = plt.subplots(figsize=(12, 10))
    im = ax.imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1)
    plt.colorbar(im, ax=ax, fraction=0.03)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    short = [f.replace("_", "\n") for f in FEATURE_NAMES]
    ax.set_xticklabels(short, rotation=45, ha="right", fontsize=7)
    ax.set_yticklabels(short, fontsize=7)

    for i in range(n):
        for j in range(n):
            val = corr[i, j]
            if abs(val) > 0.4:
                ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                        fontsize=6, color="white" if abs(val) > 0.7 else "black")

    ax.set_title("Pearson correlation matrix of topology features")
    plt.tight_layout()
    plt.savefig("features_correlation.png", dpi=150)
    print("Saved features_correlation.png")
    plt.show()


def plot_dendrogram(samples, matrix):
    names = [s["_name"] for s in samples]
    mu    = matrix.mean(axis=0)
    std   = matrix.std(axis=0)
    std[std == 0] = 1.0
    Z_norm = (matrix - mu) / std

    linked = linkage(Z_norm, method="ward")

    fig, ax = plt.subplots(figsize=(max(8, len(names) * 0.6), 5))
    dendrogram(linked, labels=names, ax=ax, leaf_rotation=45, leaf_font_size=9,
               color_threshold=0.7 * max(linked[:, 2]))
    ax.set_title("Hierarchical clustering of samples by topology features (Ward)")
    ax.set_ylabel("Distance")
    plt.tight_layout()
    plt.savefig("features_dendrogram.png", dpi=150)
    print("Saved features_dendrogram.png")
    plt.show()


def plot_scatter_matrix(samples, matrix):
    """Pairwise scatter for the 6 most variable features."""
    names    = [s["_name"] for s in samples]
    variances = matrix.var(axis=0)
    top6     = np.argsort(variances)[-6:][::-1]
    top_names = [FEATURE_NAMES[i] for i in top6]
    sub       = matrix[:, top6]

    n    = len(top6)
    fig, axes = plt.subplots(n, n, figsize=(12, 12))
    fig.suptitle("Pairwise scatter — 6 most variable features", fontsize=12)
    colors = plt.cm.tab20(np.linspace(0, 1, len(names)))

    for r in range(n):
        for c in range(n):
            ax = axes[r][c]
            if r == c:
                ax.hist(sub[:, r], bins=max(5, len(samples) // 3),
                        color="steelblue", edgecolor="white")
                ax.set_title(top_names[r].replace("_", "\n"), fontsize=7)
            else:
                ax.scatter(sub[:, c], sub[:, r], c=colors[:len(names)], s=30, zorder=3)
                for i, nm in enumerate(names):
                    ax.annotate(nm, (sub[i, c], sub[i, r]), fontsize=5,
                                xytext=(2, 2), textcoords="offset points")
            if c == 0:
                ax.set_ylabel(top_names[r].replace("_", "\n"), fontsize=6)
            if r == n - 1:
                ax.set_xlabel(top_names[c].replace("_", "\n"), fontsize=6)
            ax.tick_params(labelsize=6)

    plt.tight_layout()
    plt.savefig("features_scatter_matrix.png", dpi=150)
    print("Saved features_scatter_matrix.png")
    plt.show()


def plot_multi(samples: list[dict]):
    matrix = _build_matrix(samples)

    # Also show a combined radar overlay for all samples
    fig, axes = plt.subplots(1, 2, figsize=(14, 6),
                             subplot_kw={"polar": True})
    fig.suptitle("Feature radar overlay — all samples", fontsize=13)
    cmap = plt.cm.tab20(np.linspace(0, 1, len(samples)))

    ref_max_hist = {f: max(s.get(f, 0.0) for s in samples) or 1.0
                    for f in GROUPS["Histogram"]}
    ref_max_nw   = {f: max(s.get(f, 0.0) for s in samples) or 1.0
                    for f in GROUPS["Node"] + GROUPS["Walk/Bubble"]}

    handles = []
    for s, col in zip(samples, cmap):
        hv = normalise_for_radar(s, GROUPS["Histogram"], ref_max_hist)
        radar_chart(axes[0],
                    [f.replace("_", "\n") for f in GROUPS["Histogram"]],
                    hv, "Histogram features", color=col)
        nv = normalise_for_radar(s, GROUPS["Node"] + GROUPS["Walk/Bubble"], ref_max_nw)
        radar_chart(axes[1],
                    [f.replace("_", "\n") for f in GROUPS["Node"] + GROUPS["Walk/Bubble"]],
                    nv, "Node / Walk features", color=col)
        handles.append(mpatches.Patch(color=col, label=s["_name"]))

    fig.legend(handles=handles, loc="lower center", ncol=min(6, len(samples)),
               fontsize=8, bbox_to_anchor=(0.5, -0.05))
    plt.tight_layout()
    plt.savefig("features_radar_overlay.png", dpi=150, bbox_inches="tight")
    print("Saved features_radar_overlay.png")
    plt.show()

    plot_pca(samples, matrix)
    plot_correlation(matrix)
    plot_dendrogram(samples, matrix)
    if len(samples) >= 4:
        plot_scatter_matrix(samples, matrix)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    target = sys.argv[1]

    if os.path.isdir(target):
        samples = load_dir(target)
        if len(samples) == 1:
            plot_single(samples[0])
        else:
            plot_multi(samples)
    elif os.path.isfile(target):
        plot_single(load_one(target))
    else:
        sys.exit(f"Not a file or directory: {target}")
