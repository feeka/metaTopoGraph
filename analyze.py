#!/usr/bin/env python3
"""
analyze.py — MetaTopoGraph feature analysis

Usage: python analyze.py --input <folder> --output <folder>

Reads every features_*.json in the input folder, computes descriptive
statistics, correlations, outliers, and binomial confidence intervals for
bubble fractions, then writes 5 CSVs and 5 plots to the output folder.
"""

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import beta as beta_dist
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

# ---------------------------------------------------------------------------
# Feature config
# ---------------------------------------------------------------------------

FEATURES = [
    # histogram (5)
    "valley_depth",
    "mult_1_fraction",
    "mean_node_multiplicity",
    "multiplicity_cv",
    "n_signal_modes",
    # node (6)
    "tip_density",
    "branching_node_fraction",
    "linear_node_fraction",
    "high_degree_node_fraction",
    "mult_at_tips",
    "mult_at_branches",
    # walk (4)
    "mean_tip_length",
    "bubble_density",
    "error_bubble_fraction",
    "balanced_bubble_fraction",
]

# Raw count columns kept for confidence interval computation only.
COUNT_COLS = ["n_tips_walked", "n_branching_sampled", "n_bubbles_found"]

# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def load_data(input_dir: Path) -> pd.DataFrame:
    rows = []
    for f in sorted(input_dir.glob("features_*.json")):
        biome = f.stem.replace("features_", "")
        with open(f) as fh:
            d = json.load(fh)
        row = {"biome": biome}
        for section in ("hist", "node", "walk"):
            if section in d:
                for k, v in d[section].items():
                    row[k] = v
        rows.append(row)
    if not rows:
        sys.exit(f"No features_*.json files found in {input_dir}")
    df = pd.DataFrame(rows).set_index("biome")
    return df

# ---------------------------------------------------------------------------
# 1. Summary statistics
# ---------------------------------------------------------------------------

def compute_summary(df: pd.DataFrame) -> pd.DataFrame:
    feat = df[FEATURES]
    summary = pd.DataFrame({
        "mean":   feat.mean(),
        "median": feat.median(),
        "std":    feat.std(),
        "min":    feat.min(),
        "max":    feat.max(),
        "p25":    feat.quantile(0.25),
        "p75":    feat.quantile(0.75),
        "cv":     feat.std() / feat.mean().replace(0, np.nan),
    })
    return summary

# ---------------------------------------------------------------------------
# 2. Correlations
# ---------------------------------------------------------------------------

def compute_correlations(df: pd.DataFrame):
    feat = df[FEATURES]
    pearson  = feat.corr(method="pearson")
    spearman = feat.corr(method="spearman")
    return pearson, spearman

# ---------------------------------------------------------------------------
# 3. Outliers (beyond 2 SD from the mean per feature)
# ---------------------------------------------------------------------------

def compute_outliers(df: pd.DataFrame) -> pd.DataFrame:
    feat = df[FEATURES]
    mean = feat.mean()
    std  = feat.std()
    records = []
    for biome in feat.index:
        for col in FEATURES:
            if std[col] == 0:
                continue
            z = (feat.loc[biome, col] - mean[col]) / std[col]
            if abs(z) > 2.0:
                records.append({
                    "biome":   biome,
                    "feature": col,
                    "value":   feat.loc[biome, col],
                    "z_score": round(z, 3),
                })
    return pd.DataFrame(records)

# ---------------------------------------------------------------------------
# 4. Sampling reliability (exact binomial confidence intervals)
# ---------------------------------------------------------------------------

def clopper_pearson(k: int, n: int, alpha: float = 0.05):
    """Clopper-Pearson exact interval via the beta-binomial relationship."""
    if n == 0:
        return 0.0, 1.0
    lo = beta_dist.ppf(alpha / 2,       k,     n - k + 1) if k > 0 else 0.0
    hi = beta_dist.ppf(1 - alpha / 2,   k + 1, n - k)     if k < n else 1.0
    return lo, hi

def compute_sampling_reliability(df: pd.DataFrame) -> pd.DataFrame:
    records = []
    for biome in df.index:
        n_branching = int(df.loc[biome, "n_branching_sampled"]) if "n_branching_sampled" in df.columns else 0
        n_bubbles   = int(df.loc[biome, "n_bubbles_found"])     if "n_bubbles_found"     in df.columns else 0

        # bubble_density: proportion of branching edges that produced a bubble
        bd_val  = df.loc[biome, "bubble_density"]
        k_bd    = n_bubbles
        lo, hi  = clopper_pearson(k_bd, n_branching)
        records.append({"biome": biome, "feature": "bubble_density",
                         "value": bd_val, "ci_lo": lo, "ci_hi": hi,
                         "n": n_branching, "k": k_bd})

        # error_bubble_fraction: proportion of bubbles that are error-type
        ef_val  = df.loc[biome, "error_bubble_fraction"]
        k_err   = round(ef_val * n_bubbles)
        lo, hi  = clopper_pearson(k_err, n_bubbles)
        records.append({"biome": biome, "feature": "error_bubble_fraction",
                         "value": ef_val, "ci_lo": lo, "ci_hi": hi,
                         "n": n_bubbles, "k": k_err})

        # balanced_bubble_fraction: proportion of bubbles that are balanced
        bf_val  = df.loc[biome, "balanced_bubble_fraction"]
        k_bal   = round(bf_val * n_bubbles)
        lo, hi  = clopper_pearson(k_bal, n_bubbles)
        records.append({"biome": biome, "feature": "balanced_bubble_fraction",
                         "value": bf_val, "ci_lo": lo, "ci_hi": hi,
                         "n": n_bubbles, "k": k_bal})

    return pd.DataFrame(records)

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def _hierarchical_order(corr_matrix: pd.DataFrame) -> list:
    """Return column order after hierarchical clustering on |correlation|."""
    dist = 1.0 - corr_matrix.abs().values
    np.fill_diagonal(dist, 0.0)
    dist_condensed = squareform(dist, checks=False)
    Z = linkage(dist_condensed, method="average")
    return corr_matrix.columns[leaves_list(Z)].tolist()

def plot_feature_distributions(df: pd.DataFrame, out: Path):
    """One histogram + one boxplot per feature, laid out in a grid."""
    n     = len(FEATURES)
    ncols = 3                          # feature-pairs per row
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols * 2, figsize=(18, nrows * 3))
    axes = axes.flatten()

    for idx, feat in enumerate(FEATURES):
        ax_hist = axes[idx * 2]
        ax_box  = axes[idx * 2 + 1]
        data = df[feat].dropna()

        ax_hist.hist(data, bins=min(10, len(data)), color="steelblue", edgecolor="white")
        ax_hist.set_title(feat, fontsize=8)
        ax_hist.set_ylabel("count", fontsize=7)
        ax_hist.tick_params(labelsize=7)

        sns.boxplot(y=data, ax=ax_box, color="steelblue")
        ax_box.set_xlabel(feat, fontsize=7)
        ax_box.tick_params(labelsize=7)

    for ax in axes[n * 2:]:
        ax.set_visible(False)

    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)

def plot_correlation_heatmap(pearson: pd.DataFrame, out: Path):
    """Pearson correlation heatmap with hierarchical column/row ordering."""
    order    = _hierarchical_order(pearson)
    reordered = pearson.loc[order, order]
    fig, ax  = plt.subplots(figsize=(10, 8))
    sns.heatmap(reordered, ax=ax, cmap="coolwarm", center=0,
                vmin=-1, vmax=1, annot=True, fmt=".2f",
                annot_kws={"size": 6}, linewidths=0.4)
    ax.set_title("Pearson correlation (hierarchical order)")
    ax.tick_params(axis="x", rotation=45, labelsize=7)
    ax.tick_params(axis="y", rotation=0,  labelsize=7)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)

def plot_biome_fingerprints(df: pd.DataFrame, out: Path):
    """Heatmap of z-scored feature values, one row per biome."""
    feat = df[FEATURES]
    z    = (feat - feat.mean()) / feat.std().replace(0, np.nan)
    fig, ax = plt.subplots(figsize=(12, max(5, len(df) * 0.55)))
    sns.heatmap(z, ax=ax, cmap="RdBu_r", center=0,
                linewidths=0.3, annot=True, fmt=".2f",
                annot_kws={"size": 6})
    ax.set_title("Biome feature fingerprints (z-scored)")
    ax.tick_params(axis="x", rotation=45, labelsize=7)
    ax.tick_params(axis="y", rotation=0,  labelsize=7)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)

def plot_parallel_coordinates(df: pd.DataFrame, out: Path):
    """One line per biome across all features (normalised 0–1)."""
    feat      = df[FEATURES].copy().astype(float)
    feat_min  = feat.min()
    feat_rng  = (feat.max() - feat_min).replace(0, np.nan)
    feat_norm = (feat - feat_min) / feat_rng

    n_biomes = len(df)
    cmap     = plt.get_cmap("tab20", n_biomes)
    colors   = [cmap(i) for i in range(n_biomes)]
    x        = list(range(len(FEATURES)))

    fig, ax = plt.subplots(figsize=(16, 5))
    for i, biome in enumerate(feat_norm.index):
        vals = feat_norm.loc[biome, FEATURES].values.astype(float)
        ax.plot(x, vals, color=colors[i], label=biome, linewidth=1.2, alpha=0.85)

    ax.set_xticks(x)
    ax.set_xticklabels(FEATURES, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("Normalised value [0–1]", fontsize=8)
    ax.set_title("Parallel coordinates — one line per biome")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)

def plot_feature_pairs(df: pd.DataFrame, out: Path):
    """Scatter-plot matrix of all feature pairs, coloured by biome."""
    feat = df[FEATURES].copy()
    feat["biome"] = feat.index
    g = sns.pairplot(feat, hue="biome", corner=True,
                     plot_kws={"s": 35, "alpha": 0.85},
                     diag_kind="hist",
                     height=1.1,
                     palette="tab20")
    g.figure.suptitle("Feature pairs", y=1.01, fontsize=10)
    g.figure.savefig(out, dpi=110, bbox_inches="tight")
    plt.close(g.figure)

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="MetaTopoGraph feature analysis — produces CSVs and plots."
    )
    parser.add_argument("--input",  required=True,
                        help="Folder containing features_*.json files")
    parser.add_argument("--output", required=True,
                        help="Folder to write CSVs and plots")
    args = parser.parse_args()

    input_dir  = Path(args.input)
    output_dir = Path(args.output)

    if not input_dir.is_dir():
        sys.exit(f"Input folder not found: {input_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading data from {input_dir} ...")
    df = load_data(input_dir)
    print(f"  {len(df)} samples: {list(df.index)}")

    # ---- CSVs ----
    print("Computing summary statistics ...")
    summary = compute_summary(df)
    summary.to_csv(output_dir / "feature_summary.csv")

    print("Computing correlations ...")
    pearson, spearman = compute_correlations(df)
    pearson.to_csv(output_dir  / "correlation_pearson.csv")
    spearman.to_csv(output_dir / "correlation_spearman.csv")

    print("Identifying outliers ...")
    outliers = compute_outliers(df)
    outliers.to_csv(output_dir / "outliers.csv", index=False)

    print("Computing sampling reliability ...")
    reliability = compute_sampling_reliability(df)
    reliability.to_csv(output_dir / "sampling_reliability.csv", index=False)

    # ---- Plots ----
    sns.set_theme(style="whitegrid")

    print("Plotting feature distributions ...")
    plot_feature_distributions(df, output_dir / "feature_distributions.png")

    print("Plotting correlation heatmap ...")
    plot_correlation_heatmap(pearson, output_dir / "correlation_heatmap.png")

    print("Plotting biome fingerprints ...")
    plot_biome_fingerprints(df, output_dir / "biome_fingerprints.png")

    print("Plotting parallel coordinates ...")
    plot_parallel_coordinates(df, output_dir / "parallel_coordinates.png")

    print("Plotting feature pairs ...")
    plot_feature_pairs(df, output_dir / "feature_pairs.png")

    print(f"\nDone. Outputs written to {output_dir}/")


if __name__ == "__main__":
    main()
