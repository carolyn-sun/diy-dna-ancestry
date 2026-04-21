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

def _load_labels(ref_dir: str) -> pd.DataFrame | None:
    """Load the HGDP population label table."""
    label_path = Path(ref_dir) / "hgdp_pop_labels.tsv"
    if not label_path.exists():
        console.print(f"  [dim]Label file not found, skipping population colours: {label_path}[/dim]")
        return None
    try:
        df = pd.read_csv(label_path, sep="\t")
        col_map = {}
        for col in df.columns:
            c = col.lower()
            if "sample" in c or c == "id":
                col_map[col] = "sample_id"
            elif "pop" in c or "population" in c:
                col_map[col] = "population"
            elif "region" in c or "continent" in c:
                col_map[col] = "region"
        return df.rename(columns=col_map)
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

    colors  = _ADMIX_COLORS[:k]
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

    # Locate merged FAM file
    merged_fam = None
    if admix_results:
        first = next(iter(admix_results.values()))
        stem = Path(first["q_file"]).stem.rsplit(".", 2)[0]  # e.g. "merged"
        for candidate in [
            str(Path(first["q_file"]).parent / (stem + ".fam")),
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
        else:
            console.print(f"  [yellow]merged.fam not found — skipping K={k} bar chart[/yellow]")

    if len(admix_results) > 1:
        console.print("  [bold]Plot: CV Error curve[/bold]")
        _plot_cv_error(admix_results, out_dir=out_dir)

    console.print(f"\n  [bold green]All plots saved to: {os.path.abspath(out_dir)}[/bold green]")
