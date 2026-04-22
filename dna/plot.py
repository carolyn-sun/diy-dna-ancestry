"""
plot.py — Visualisation module (matplotlib)

Produces:
  1. PCA scatter plots   pca_PC1_PC2.png / pca_PC3_PC4.png
  2. ADMIXTURE bar charts admixture_K{k}.png
  3. CV error curve      cv_error.png (when multiple K values are present)
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")  # headless — no display window

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from rich.console import Console

console = Console()

# ──────────────────────────────────────────────────────────────────────────────
# Global style
# ──────────────────────────────────────────────────────────────────────────────

plt.rcParams.update({
    "figure.facecolor": "#0f1117",
    "axes.facecolor":   "#161b22",
    "axes.edgecolor":   "#30363d",
    "axes.labelcolor":  "#c9d1d9",
    "axes.titlecolor":  "#e6edf3",
    "xtick.color":      "#8b949e",
    "ytick.color":      "#8b949e",
    "text.color":       "#c9d1d9",
    "grid.color":       "#21262d",
    "grid.alpha":       0.6,
    "legend.facecolor": "#161b22",
    "legend.edgecolor": "#30363d",
    "font.family":      "DejaVu Sans",
    "font.size":        10,
    "axes.titlesize":   13,
    "figure.dpi":       150,
})

# Colour palette keyed by HGDP geographic region
REGION_PALETTE = {
    "Africa":              "#F4A261",
    "America":             "#E76F51",
    "Central_South_Asia":  "#2A9D8F",
    "East_Asia":           "#4CC9F0",
    "Europe":              "#7B68EE",
    "Middle_East":         "#F72585",
    "Oceania":             "#90E0EF",
    "Unknown":             "#8b949e",
    "USER":                "#FFFFFF",   # user sample: white star
}


# ──────────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────────

# Map SuperPop abbreviation codes → REGION_PALETTE keys
_SUPERPOP_TO_REGION: dict[str, str] = {
    "AFR": "Africa",
    "AMR": "America",
    "EAS": "East_Asia",
    "EUR": "Europe",
    "CSA": "Central_South_Asia",
    "SAS": "Central_South_Asia",
    "MID": "Middle_East",
    "ME":  "Middle_East",
    "MDE": "Middle_East",
    "OCE": "Oceania",
}


def _load_labels(ref_dir: str) -> pd.DataFrame | None:
    """Load the HGDP population label table.

    The Zenodo label file (hgdp_pop_labels.tsv) has columns:
        Sample | SuperPop | Project
    where SuperPop uses short codes (AFR, EAS, EUR, …).
    We normalise column names and expand codes to REGION_PALETTE keys.
    """
    label_path = Path(ref_dir) / "hgdp_pop_labels.tsv"
    if not label_path.exists():
        console.print(f"  [dim]Label file not found, skipping colours: {label_path}[/dim]")
        return None
    try:
        df = pd.read_csv(label_path, sep="\t")
        col_map = {}
        for col in df.columns:
            c = col.lower().replace(" ", "").replace("_", "")
            if "sample" in c or c in ("id", "s"):
                col_map[col] = "sample_id"
            elif "superpop" in c:
                col_map[col] = "region"   # SuperPop IS the region here
            elif "region" in c or "continent" in c:
                col_map[col] = "region"
            elif "pop" in c or "population" in c:
                col_map[col] = "population"
        df = df.rename(columns=col_map)
        # Expand SuperPop abbreviation → full region name expected by REGION_PALETTE
        if "region" in df.columns:
            df["region"] = df["region"].apply(
                lambda x: _SUPERPOP_TO_REGION.get(str(x).upper(), str(x))
            )
        return df
    except Exception as e:
        console.print(f"  [yellow]Failed to read label file: {e}[/yellow]")
        return None


def _assign_colors(fid_series: pd.Series, iid_series: pd.Series,
                   labels_df: Optional[pd.DataFrame]) -> list[str]:
    """Assign a hex colour to each sample."""
    colors = []
    for fid, iid in zip(fid_series, iid_series):
        if str(fid) == "USER":
            colors.append(REGION_PALETTE["USER"])
            continue
        if labels_df is not None and "sample_id" in labels_df.columns and "region" in labels_df.columns:
            row = labels_df[labels_df["sample_id"] == iid]
            if not row.empty:
                region = row.iloc[0]["region"]
                colors.append(REGION_PALETTE.get(str(region), REGION_PALETTE["Unknown"]))
                continue
        colors.append(REGION_PALETTE["Unknown"])
    return colors


def _save(fig: plt.Figure, path: str) -> None:
    fig.savefig(path, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    console.print(f"  [green]✓[/green] Saved: [cyan]{path}[/cyan]")



# ──────────────────────────────────────────────────────────────────────────────
# PCA-KNN ancestry proximity chart
# ──────────────────────────────────────────────────────────────────────────────

def _plot_pca_knn_ancestry(
    pca_results: dict, ref_dir: str, out_dir: str, n_neighbors: int = 30
) -> None:
    """
    Estimate ancestry from PCA position using K nearest neighbours.

    Finds the n_neighbors closest reference samples in the full PC space
    (all available PCs, usually 10) and reports their regional composition
    as a horizontal bar chart.  This directly reflects the user's position
    in the PCA scatter plots — more intuitive than ADMIXTURE-derived pie charts.
    """
    from collections import Counter

    # ── Load PCA data ─────────────────────────────────────────────────────────
    if "eigenvec_df" in pca_results:
        df = pca_results["eigenvec_df"]
    else:
        eigenvec_path = pca_results["eigenvec"]
        df = pd.read_csv(eigenvec_path, sep=r"\s+", header=None)
        n_pcs = df.shape[1] - 2
        df.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, n_pcs + 1)]

    pc_cols = [c for c in df.columns if c.startswith("PC")]

    # ── Separate user and reference samples ───────────────────────────────────
    user_mask = df["FID"] == "USER"
    ref_mask  = ~user_mask

    if not user_mask.any():
        console.print("  [yellow]PCA KNN: no USER sample found — skipping[/yellow]")
        return

    user_coords = df.loc[user_mask, pc_cols].values[0]   # (n_pcs,)
    ref_coords  = df.loc[ref_mask,  pc_cols].values       # (n_ref, n_pcs)
    ref_iids    = df.loc[ref_mask, "IID"].values

    # ── Load region labels ────────────────────────────────────────────────────
    labels_df = _load_labels(ref_dir)
    if (labels_df is None
            or "sample_id" not in labels_df.columns
            or "region"    not in labels_df.columns):
        console.print("  [yellow]PCA KNN: no region labels found — skipping[/yellow]")
        return

    iid_to_region = dict(zip(labels_df["sample_id"], labels_df["region"]))

    # ── KNN in full PC space (Euclidean) ──────────────────────────────────────
    dists = np.linalg.norm(ref_coords - user_coords, axis=1)
    k     = min(n_neighbors, len(dists))
    nn_idx = np.argsort(dists)[:k]

    nn_regions = [iid_to_region.get(ref_iids[i], "Unknown") for i in nn_idx]
    nn_dists   = dists[nn_idx]

    region_counts = Counter(nn_regions)
    total = sum(region_counts.values())

    # ── Order regions by proximity fraction ───────────────────────────────────
    region_order = [
        "East_Asia", "Europe", "Middle_East", "Central_South_Asia",
        "America", "Africa", "Oceania", "Unknown",
    ]
    present = [(r, region_counts[r] / total)
               for r in region_order if r in region_counts]
    present += [(r, v / total) for r, v in region_counts.items()
                if r not in region_order and v > 0]
    present = sorted(present, key=lambda x: x[1])   # ascending → longest bar on top

    regions = [r for r, _ in present]
    values  = [v for _, v in present]
    colors  = [REGION_PALETTE.get(r, REGION_PALETTE["Unknown"]) for r in regions]

    # ── Draw ──────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(9, max(3, len(regions) * 0.7 + 1.5)))
    fig.patch.set_facecolor("#0f1117")
    ax.set_facecolor("#0f1117")

    bars = ax.barh(
        [r.replace("_", " ") for r in regions], values,
        color=colors, edgecolor="#0f1117", linewidth=1.2, height=0.65,
    )
    for bar, val in zip(bars, values):
        ax.text(
            bar.get_width() + 0.005,
            bar.get_y() + bar.get_height() / 2,
            f"{val * 100:.1f}%",
            va="center", ha="left", fontsize=11, color="#e6edf3", fontweight="bold",
        )

    ax.set_xlim(0, max(values) * 1.28)
    ax.xaxis.set_visible(False)
    ax.tick_params(colors="#8b949e", labelsize=11)
    ax.set_yticklabels([r.replace("_", " ") for r in regions],
                       color="#e6edf3", fontsize=11)
    for spine in ax.spines.values():
        spine.set_edgecolor("#30363d")

    ax.set_title(
        f"Ancestry by PCA Proximity  ({k} nearest neighbours, {len(pc_cols)} PCs)",
        fontsize=13, pad=14, color="#e6edf3",
    )
    fig.text(
        0.5, 0.01,
        f"Euclidean distance in {len(pc_cols)}-D PC space  ·  "
        f"nearest: {nn_dists[0]:.4f}  ·  furthest of {k}: {nn_dists[-1]:.4f}",
        ha="center", fontsize=8, color="#8b949e",
    )
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    _save(fig, str(Path(out_dir) / "pca_ancestry_knn.png"))


# ──────────────────────────────────────────────────────────────────────────────
# PCA scatter plot
# ──────────────────────────────────────────────────────────────────────────────

def _plot_pca(pca_results: dict, ref_dir: str, out_dir: str, pc_x: int = 1, pc_y: int = 2) -> None:
    """Draw a PCA scatter plot for two specified components."""
    if "eigenvec_df" in pca_results:
        df = pca_results["eigenvec_df"]
    else:
        eigenvec_path = pca_results["eigenvec"]
        df = pd.read_csv(eigenvec_path, sep=r"\s+", header=None)
        n_pcs = df.shape[1] - 2
        df.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, n_pcs + 1)]

    eigenval_df = pca_results.get("eigenval_df")
    if eigenval_df is None and "eigenval" in pca_results:
        eigenval_df = pd.read_csv(pca_results["eigenval"], header=None, names=["eigenval"])
        total = eigenval_df["eigenval"].sum()
        eigenval_df["var_pct"] = eigenval_df["eigenval"] / total * 100

    labels_df = _load_labels(ref_dir)
    colors = _assign_colors(df["FID"], df["IID"], labels_df)

    pc_x_col, pc_y_col = f"PC{pc_x}", f"PC{pc_y}"
    x_label, y_label = f"PC{pc_x}", f"PC{pc_y}"
    if eigenval_df is not None:
        try:
            x_pct = eigenval_df.iloc[pc_x - 1]["var_pct"]
            y_pct = eigenval_df.iloc[pc_y - 1]["var_pct"]
            x_label = f"PC{pc_x} ({x_pct:.1f}%)"
            y_label = f"PC{pc_y} ({y_pct:.1f}%)"
        except (IndexError, KeyError):
            pass

    fig, ax = plt.subplots(figsize=(10, 8))

    mask_ref  = df["FID"] != "USER"
    mask_user = df["FID"] == "USER"
    ref_colors = [c for c, m in zip(colors, mask_ref) if m]

    ax.scatter(df.loc[mask_ref, pc_x_col], df.loc[mask_ref, pc_y_col],
               c=ref_colors, s=18, alpha=0.7, linewidths=0, zorder=2)

    if mask_user.any():
        ax.scatter(df.loc[mask_user, pc_x_col], df.loc[mask_user, pc_y_col],
                   c=REGION_PALETTE["USER"], s=220, marker="*", zorder=5,
                   linewidths=0.8, edgecolors="#000000", label="Your sample")

    # Legend patches
    seen_regions: set[str] = set()
    if labels_df is not None and "sample_id" in labels_df.columns and "region" in labels_df.columns:
        region_iid_map = dict(zip(labels_df["sample_id"], labels_df["region"]))
        for iid, fid in zip(df["IID"], df["FID"]):
            if str(fid) != "USER":
                seen_regions.add(str(region_iid_map.get(iid, "Unknown")))
    else:
        seen_regions = set(REGION_PALETTE.keys()) - {"USER"}

    patches = [mpatches.Patch(color=REGION_PALETTE.get(r, REGION_PALETTE["Unknown"]), label=r)
               for r in sorted(seen_regions) if r in REGION_PALETTE]
    if mask_user.any():
        patches.append(mpatches.Patch(color=REGION_PALETTE["USER"], label="Your sample ★"))
    ax.legend(handles=patches, loc="upper right", fontsize=8, framealpha=0.8, ncol=2)

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(f"PCA — HGDP reference populations + your sample (PC{pc_x} vs PC{pc_y})")
    ax.grid(True, ls="--", lw=0.5)
    fig.tight_layout()
    _save(fig, str(Path(out_dir) / f"pca_PC{pc_x}_PC{pc_y}.png"))


# ──────────────────────────────────────────────────────────────────────────────
# ADMIXTURE bar chart
# ──────────────────────────────────────────────────────────────────────────────

_ADMIX_COLORS = [
    "#4CC9F0", "#F72585", "#7B68EE", "#F4A261", "#2A9D8F",
    "#E76F51", "#90E0EF", "#FAB5E0", "#C77DFF", "#FCC8A5",
    "#06D6A0", "#FFD166", "#EF476F", "#118AB2", "#073B4C",
]


def _plot_admixture(k: int, q_file: str, fam_file: str, ref_dir: str, out_dir: str) -> None:
    """Draw an ADMIXTURE bar chart for a single K value."""
    q_df = pd.read_csv(q_file, sep=r"\s+", header=None)
    q_df.columns = [f"K{i+1}" for i in range(k)]

    fam_df = pd.read_csv(fam_file, sep=r"\s+", header=None,
                         names=["FID", "IID", "PAT", "MAT", "SEX", "PHEN"])
    q_df = pd.concat([fam_df[["FID", "IID"]], q_df], axis=1)

    labels_df = _load_labels(ref_dir)
    region_order = [
        "Africa", "America", "Central_South_Asia",
        "East_Asia", "Europe", "Middle_East", "Oceania", "Unknown", "USER",
    ]
    if labels_df is not None and "sample_id" in labels_df.columns and "region" in labels_df.columns:
        q_df = q_df.merge(
            labels_df[["sample_id", "region"]].rename(columns={"sample_id": "IID"}),
            on="IID", how="left",
        )
        q_df["region"] = q_df["region"].fillna(
            q_df["FID"].apply(lambda x: "USER" if x == "USER" else "Unknown")
        )
        q_df["_order"] = q_df["region"].map({r: i for i, r in enumerate(region_order)}).fillna(99)
        q_df = q_df.sort_values(["_order", "region", "IID"]).reset_index(drop=True)
    else:
        q_df["region"] = q_df["FID"].apply(lambda x: "USER" if x == "USER" else "Unknown")

    # Build color list — cycle through extras if K > len(_ADMIX_COLORS)
    if k <= len(_ADMIX_COLORS):
        colors = _ADMIX_COLORS[:k]
    else:
        import matplotlib.cm as _cm
        tab20 = [_cm.tab20(i / 20) for i in range(20)]
        extras = ["#{:02x}{:02x}{:02x}".format(int(r*255), int(g*255), int(b*255))
                  for r, g, b, _ in tab20]
        colors = (_ADMIX_COLORS + extras)[:k]
    k_cols  = [f"K{i+1}" for i in range(k)]

    fig, ax = plt.subplots(figsize=(max(14, len(q_df) * 0.04), 5))
    bottom = np.zeros(len(q_df))

    for i, col in enumerate(k_cols):
        ax.bar(range(len(q_df)), q_df[col].values, bottom=bottom,
               color=colors[i], width=1.0, linewidth=0, label=f"Component {i+1}")
        bottom += q_df[col].values

    # Mark user sample
    for idx in q_df[q_df["FID"] == "USER"].index.tolist():
        ax.axvline(x=idx, color="white", lw=1.5, alpha=0.9)
        ax.text(idx, 1.02, "★ you", ha="center", va="bottom",
                fontsize=7, color="white", fontweight="bold")

    # Region separators
    if "region" in q_df.columns:
        bounds: list[tuple[int, int, str]] = []
        prev, start = q_df.iloc[0]["region"], 0
        for i, row in q_df.iterrows():
            if row["region"] != prev:
                bounds.append((start, i - 1, prev))
                start, prev = i, row["region"]
        bounds.append((start, len(q_df) - 1, prev))

        for s, e, region in bounds:
            ax.axvline(x=s - 0.5, color="#30363d", lw=1.2)
            ax.text((s + e) / 2, -0.08, region.replace("_", "\n"),
                    ha="center", va="top", fontsize=6.5, color="#8b949e",
                    transform=ax.get_xaxis_transform(), clip_on=False)

    ax.set_xlim(-0.5, len(q_df) - 0.5)
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["0%", "25%", "50%", "75%", "100%"])
    ax.set_xticks([])
    ax.set_title(f"ADMIXTURE unsupervised clustering (K={k}) — {len(q_df)} samples")
    ax.set_ylabel("Ancestry proportion")
    ax.legend(loc="upper left", fontsize=8, ncol=min(k, 5),
              framealpha=0.8, bbox_to_anchor=(0, 1.15))
    fig.tight_layout()
    _save(fig, str(Path(out_dir) / f"admixture_K{k}.png"))


# ──────────────────────────────────────────────────────────────────────────────
# Ancestry pie chart (geographic breakdown for the user sample)
# ──────────────────────────────────────────────────────────────────────────────

def _plot_ancestry_pie(k: int, q_file: str, fam_file: str,
                       ref_dir: str, out_dir: str) -> None:
    """
    Produce a pie chart showing the user's ancestry broken down by geographic
    region.

    Method
    ------
    For each NMF / ADMIXTURE component i we compute its *centroid* — the
    average proportion of that component across all reference samples belonging
    to each geographic region.  Each component is then assigned to its dominant
    region.  The user's Q values are mapped through this assignment to produce
    a region-level percentage.

    If a component cannot be uniquely assigned (tie), it is split proportionally.
    """
    # ── Load Q and FAM ────────────────────────────────────────────────────────
    q_raw = pd.read_csv(q_file, sep=r"\s+", header=None)
    q_raw.columns = [f"K{i+1}" for i in range(k)]
    fam_df = pd.read_csv(fam_file, sep=r"\s+", header=None,
                         names=["FID", "IID", "PAT", "MAT", "SEX", "PHEN"])
    if len(fam_df) != len(q_raw):
        fam_df = fam_df.iloc[:len(q_raw)].reset_index(drop=True)
    q_df = pd.concat([fam_df[["FID", "IID"]], q_raw], axis=1)

    # ── Load region labels ────────────────────────────────────────────────────
    labels_df = _load_labels(ref_dir)
    if labels_df is not None and "sample_id" in labels_df.columns and "region" in labels_df.columns:
        q_df = q_df.merge(
            labels_df[["sample_id", "region"]].rename(columns={"sample_id": "IID"}),
            on="IID", how="left",
        )
        q_df["region"] = q_df["region"].fillna(
            q_df["FID"].apply(lambda x: "USER" if x == "USER" else "Unknown")
        )
    else:
        q_df["region"] = q_df["FID"].apply(lambda x: "USER" if x == "USER" else "Unknown")

    # ── Extract user row ──────────────────────────────────────────────────────
    user_mask = q_df["FID"] == "USER"
    if not user_mask.any():
        console.print("  [yellow]Pie chart: user sample not found in Q file — skipping[/yellow]")
        return
    user_q = q_df.loc[user_mask, [f"K{i+1}" for i in range(k)]].values[0]  # shape (k,)

    # ── Compute component centroids per region (reference samples only) ───────
    ref_df = q_df[~user_mask].copy()
    k_cols = [f"K{i+1}" for i in range(k)]
    region_order = [
        "Africa", "America", "Central_South_Asia",
        "East_Asia", "Europe", "Middle_East", "Oceania",
    ]
    present_regions = [r for r in region_order if r in ref_df["region"].values]

    # centroid_matrix[r, i] = mean Q_i for reference samples in region r
    centroid = np.zeros((len(present_regions), k))
    for ri, region in enumerate(present_regions):
        mask = ref_df["region"] == region
        if mask.sum() > 0:
            centroid[ri] = ref_df.loc[mask, k_cols].values.mean(axis=0)

    # ── Soft assignment (weighted) ────────────────────────────────────────────
    # Each component i contributes user_q[i] ancestry, split across regions
    # proportionally to centroid[:, i].  This prevents a component that
    # represents both Africa and Europe from being incorrectly 100% Africa.
    #
    #   weight[r, i] = centroid[r, i] / sum_r(centroid[r, i])
    #   ancestry[r] += user_q[i] * weight[r, i]  for all i
    #
    col_sums = centroid.sum(axis=0, keepdims=True).clip(1e-12)
    weights = centroid / col_sums               # shape (n_regions, k)
    ancestry_vec = weights @ user_q             # shape (n_regions,)

    ancestry: dict[str, float] = {
        r: float(v) for r, v in zip(present_regions, ancestry_vec)
    }

    # Drop negligible regions (< 1%)
    ancestry = {r: v for r, v in ancestry.items() if v > 0.01}
    if not ancestry:
        console.print("  [yellow]Pie chart: all ancestry proportions are zero — skipping[/yellow]")
        return

    # ── Draw pie chart ────────────────────────────────────────────────────────
    labels  = list(ancestry.keys())
    values  = np.array(list(ancestry.values()))
    colors  = [REGION_PALETTE.get(r, REGION_PALETTE["Unknown"]) for r in labels]

    fig, ax = plt.subplots(figsize=(7, 7))
    fig.patch.set_facecolor("#0f1117")
    ax.set_facecolor("#0f1117")

    wedges, texts, autotexts = ax.pie(
        values,
        labels=None,
        colors=colors,
        autopct=lambda p: f"{p:.1f}%" if p > 2 else "",
        startangle=140,
        wedgeprops={"linewidth": 1.5, "edgecolor": "#0f1117"},
        pctdistance=0.75,
    )
    for at in autotexts:
        at.set_fontsize(10)
        at.set_color("white")
        at.set_fontweight("bold")

    # Legend with percentages
    legend_labels = [f"{r.replace('_', ' ')}  {v*100:.1f}%" for r, v in zip(labels, values)]
    ax.legend(
        wedges, legend_labels,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.12),
        ncol=2,
        fontsize=10,
        framealpha=0.15,
        edgecolor="#30363d",
    )

    ax.set_title(
        f"Estimated Ancestry Composition (K={k})",
        fontsize=14, pad=18, color="#e6edf3",
    )
    # Subtitle note
    fig.text(0.5, 0.01,
             "Components assigned to dominant reference population per region",
             ha="center", fontsize=8, color="#8b949e")

    fig.tight_layout()
    _save(fig, str(Path(out_dir) / f"ancestry_pie_K{k}.png"))


# ──────────────────────────────────────────────────────────────────────────────
# CV error curve
# ──────────────────────────────────────────────────────────────────────────────

def _plot_cv_error(admix_results: dict[int, dict], out_dir: str) -> None:
    """Plot CV error vs. K."""
    ks      = sorted(admix_results.keys())
    cv_errs = [admix_results[k].get("cv_error") for k in ks]

    if all(v is None for v in cv_errs):
        console.print("  [dim]No valid CV error data — skipping curve[/dim]")
        return

    valid_ks  = [k for k, v in zip(ks, cv_errs) if v is not None]
    valid_cvs = [v for v in cv_errs if v is not None]

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(valid_ks, valid_cvs, "o-", color="#4CC9F0", lw=2, ms=8, zorder=3)

    best_k = valid_ks[int(np.argmin(valid_cvs))]
    best_v = min(valid_cvs)
    ax.scatter([best_k], [best_v], c="#F72585", s=120, zorder=4)
    ax.annotate(f"K={best_k}\nbest", xy=(best_k, best_v),
                xytext=(best_k + 0.3, best_v + (max(valid_cvs) - min(valid_cvs)) * 0.15),
                color="#F72585", fontsize=9)

    ax.set_xlabel("K (number of clusters)")
    ax.set_ylabel("CV Error")
    ax.set_title("ADMIXTURE Cross-Validation Error by K")
    ax.set_xticks(valid_ks)
    ax.grid(True, ls="--", lw=0.5)
    fig.tight_layout()
    _save(fig, str(Path(out_dir) / "cv_error.png"))


# ──────────────────────────────────────────────────────────────────────────────
# Public entry point
# ──────────────────────────────────────────────────────────────────────────────

def make_all_plots(
    pca_results: dict,
    admix_results: dict[int, dict],
    ref_dir: str,
    user_vcf: Optional[str],
    out_dir: str,
) -> None:
    """
    Generate all plots.

    Args:
        pca_results:   Dict from pca.run_pca()
        admix_results: Dict from admixture.run_admixture() (keyed by K)
        ref_dir:       HGDP reference panel directory (contains label file)
        user_vcf:      User VCF path (may be None; used only for plot titles)
        out_dir:       Output directory
    """
    os.makedirs(out_dir, exist_ok=True)

    console.print("\n  [bold]Plot: PCA (PC1 vs PC2)[/bold]")
    _plot_pca(pca_results, ref_dir=ref_dir, out_dir=out_dir, pc_x=1, pc_y=2)

    console.print("  [bold]Plot: PCA (PC3 vs PC4)[/bold]")
    _plot_pca(pca_results, ref_dir=ref_dir, out_dir=out_dir, pc_x=3, pc_y=4)

    console.print("  [bold]Plot: Ancestry by PCA proximity (KNN)[/bold]")
    _plot_pca_knn_ancestry(pca_results, ref_dir=ref_dir, out_dir=out_dir, n_neighbors=30)

    # Locate FAM file — prefer the one stored in result (matches Q row count)
    merged_fam = None
    if admix_results:
        first = next(iter(admix_results.values()))
        # _run_admixture_k / _run_nmf_fallback embed the correct fam_file path
        if "fam_file" in first and Path(first["fam_file"]).exists():
            merged_fam = first["fam_file"]
        else:
            stem = Path(first["q_file"]).stem.rsplit(".", 2)[0]
            for candidate in [
                str(Path(first["q_file"]).parent / (stem + ".fam")),
                str(Path(first["q_file"]).parent.parent / "work" / "admix_sub.fam"),
                str(Path(first["q_file"]).parent.parent / "work" / "admix_input.fam"),
                str(Path(first["q_file"]).parent.parent / "work" / "merged.fam"),
            ]:
                if Path(candidate).exists():
                    merged_fam = candidate
                    break

    for k, result in sorted(admix_results.items()):
        console.print(f"  [bold]Plot: ADMIXTURE K={k}[/bold]")
        if merged_fam:
            _plot_admixture(k=k, q_file=result["q_file"], fam_file=merged_fam,
                            ref_dir=ref_dir, out_dir=out_dir)
            console.print(f"  [bold]Plot: Ancestry pie K={k}[/bold]")
            _plot_ancestry_pie(k=k, q_file=result["q_file"], fam_file=merged_fam,
                               ref_dir=ref_dir, out_dir=out_dir)
        else:
            console.print(f"  [yellow]merged.fam not found — skipping K={k} bar chart[/yellow]")

    if len(admix_results) > 1:
        console.print("  [bold]Plot: CV Error curve[/bold]")
        _plot_cv_error(admix_results, out_dir=out_dir)

    console.print(f"\n  [bold green]All plots saved to: {os.path.abspath(out_dir)}[/bold green]")
