"""
plot.py — 可视化模块（matplotlib）

生成图表：
  1. PCA 散点图      pca.png / pca_pc3_pc4.png
  2. ADMIXTURE 条形图 admixture_K{k}.png
  3. CV error 曲线   cv_error.png（若有多个 K）
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")  # 无头环境（headless）不弹窗

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from rich.console import Console

console = Console()

# ──────────────────────────────────────────────────────────────────────────────
# 全局样式
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

# 七大洲 / 地理区域对应色板（与 HGDP 标准人群对应）
REGION_PALETTE = {
    "Africa":        "#F4A261",
    "America":       "#E76F51",
    "Central_South_Asia": "#2A9D8F",
    "East_Asia":     "#4CC9F0",
    "Europe":        "#7B68EE",
    "Middle_East":   "#F72585",
    "Oceania":       "#90E0EF",
    "Unknown":       "#8b949e",
    "USER":          "#FFFFFF",   # 用户样本：白色 + 星形
}


# ──────────────────────────────────────────────────────────────────────────────
# 工具函数
# ──────────────────────────────────────────────────────────────────────────────

def _load_labels(ref_dir: str) -> pd.DataFrame | None:
    """加载 HGDP 人群标签表。"""
    label_path = Path(ref_dir) / "hgdp_pop_labels.tsv"
    if not label_path.exists():
        console.print(f"  [dim]标签文件不存在，跳过人群着色：{label_path}[/dim]")
        return None
    try:
        df = pd.read_csv(label_path, sep="\t")
        # 尝试标准化列名
        col_map = {}
        for col in df.columns:
            c = col.lower()
            if "sample" in c or "id" in c:
                col_map[col] = "sample_id"
            elif "pop" in c or "population" in c:
                col_map[col] = "population"
            elif "region" in c or "continent" in c:
                col_map[col] = "region"
        df = df.rename(columns=col_map)
        return df
    except Exception as e:
        console.print(f"  [yellow]标签文件读取失败：{e}[/yellow]")
        return None


def _assign_colors(fid_series: pd.Series, iid_series: pd.Series,
                   labels_df: Optional[pd.DataFrame]) -> list[str]:
    """为每个样本分配颜色。"""
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
    console.print(f"  [green]✓[/green] 已保存：[cyan]{path}[/cyan]")


# ──────────────────────────────────────────────────────────────────────────────
# PCA 散点图
# ──────────────────────────────────────────────────────────────────────────────

def _plot_pca(
    pca_results: dict,
    ref_dir: str,
    out_dir: str,
    pc_x: int = 1,
    pc_y: int = 2,
) -> None:
    """绘制 PCA 散点图。"""
    # 加载数据
    if "eigenvec_df" in pca_results:
        df = pca_results["eigenvec_df"]
    else:
        eigenvec_path = pca_results["eigenvec"]
        df = pd.read_csv(eigenvec_path, sep=r"\s+", header=None)
        n_pcs = df.shape[1] - 2
        df.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, n_pcs + 1)]

    # 加载 variance explained
    eigenval_df = None
    if "eigenval_df" in pca_results:
        eigenval_df = pca_results["eigenval_df"]
    elif "eigenval" in pca_results:
        eigenval_df = pd.read_csv(pca_results["eigenval"], header=None, names=["eigenval"])
        total = eigenval_df["eigenval"].sum()
        eigenval_df["var_pct"] = eigenval_df["eigenval"] / total * 100

    # 标签
    labels_df = _load_labels(ref_dir)
    colors = _assign_colors(df["FID"], df["IID"], labels_df)

    pc_x_col = f"PC{pc_x}"
    pc_y_col = f"PC{pc_y}"

    # 轴标签（含方差贡献）
    x_label = f"PC{pc_x}"
    y_label = f"PC{pc_y}"
    if eigenval_df is not None:
        try:
            x_pct = eigenval_df.iloc[pc_x - 1]["var_pct"]
            y_pct = eigenval_df.iloc[pc_y - 1]["var_pct"]
            x_label = f"PC{pc_x} ({x_pct:.1f}%)"
            y_label = f"PC{pc_y} ({y_pct:.1f}%)"
        except (IndexError, KeyError):
            pass

    fig, ax = plt.subplots(figsize=(10, 8))

    # 参考人群点
    mask_ref  = df["FID"] != "USER"
    mask_user = df["FID"] == "USER"

    ref_colors = [c for c, m in zip(colors, mask_ref) if m]
    ax.scatter(
        df.loc[mask_ref, pc_x_col],
        df.loc[mask_ref, pc_y_col],
        c=ref_colors,
        s=18, alpha=0.7, linewidths=0,
        label="_nolegend_",
        zorder=2,
    )

    # 用户样本（星星标记）
    if mask_user.any():
        ax.scatter(
            df.loc[mask_user, pc_x_col],
            df.loc[mask_user, pc_y_col],
            c=REGION_PALETTE["USER"],
            s=220, marker="*", zorder=5,
            linewidths=0.8, edgecolors="#000000",
            label="你的样本",
        )

    # 图例：地理区域
    seen_regions: set[str] = set()
    if labels_df is not None and "sample_id" in labels_df.columns and "region" in labels_df.columns:
        region_iid_map = dict(zip(labels_df["sample_id"], labels_df["region"]))
        for iid, fid in zip(df["IID"], df["FID"]):
            if str(fid) != "USER":
                region = region_iid_map.get(iid, "Unknown")
                seen_regions.add(str(region))
    else:
        seen_regions = set(REGION_PALETTE.keys()) - {"USER"}

    legend_patches = [
        mpatches.Patch(color=REGION_PALETTE.get(r, REGION_PALETTE["Unknown"]), label=r)
        for r in sorted(seen_regions) if r in REGION_PALETTE
    ]
    if mask_user.any():
        legend_patches.append(
            mpatches.Patch(color=REGION_PALETTE["USER"], label="你的样本 ★")
        )
    ax.legend(
        handles=legend_patches,
        loc="upper right", fontsize=8,
        framealpha=0.8, ncol=2,
    )

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(f"PCA — HGDP 参考人群 + 你的样本（PC{pc_x} vs PC{pc_y}）")
    ax.grid(True, ls="--", lw=0.5)
    fig.tight_layout()

    out_path = str(Path(out_dir) / f"pca_PC{pc_x}_PC{pc_y}.png")
    _save(fig, out_path)


# ──────────────────────────────────────────────────────────────────────────────
# ADMIXTURE 条形图
# ──────────────────────────────────────────────────────────────────────────────

# 每个 K 值使用的颜色序列（最多支持 K=10）
_ADMIX_PALETTES = [
    ["#4CC9F0", "#F72585", "#7B68EE", "#F4A261", "#2A9D8F",
     "#E76F51", "#90E0EF", "#FAB5E0", "#C77DFF", "#FCC8A5"],
]


def _plot_admixture(
    k: int,
    q_file: str,
    fam_file: str,
    ref_dir: str,
    out_dir: str,
) -> None:
    """绘制单个 K 值的 ADMIXTURE 条形图。"""
    # 加载 Q 矩阵
    q_df = pd.read_csv(q_file, sep=r"\s+", header=None)
    q_df.columns = [f"K{i+1}" for i in range(k)]

    # 加载 FAM（获取样本 FID/IID）
    fam_df = pd.read_csv(
        fam_file, sep=r"\s+", header=None,
        names=["FID", "IID", "PAT", "MAT", "SEX", "PHEN"],
    )
    q_df = pd.concat([fam_df[["FID", "IID"]], q_df], axis=1)

    # 加载标签，按地理区域排序（USER 放最后）
    labels_df = _load_labels(ref_dir)
    if labels_df is not None and "sample_id" in labels_df.columns and "region" in labels_df.columns:
        q_df = q_df.merge(
            labels_df[["sample_id", "region", "population"]].rename(
                columns={"sample_id": "IID"}
            ),
            on="IID", how="left",
        )
        q_df["region"] = q_df["region"].fillna(
            q_df["FID"].apply(lambda x: "USER" if x == "USER" else "Unknown")
        )
        region_order = [
            "Africa", "America", "Central_South_Asia",
            "East_Asia", "Europe", "Middle_East", "Oceania", "Unknown", "USER",
        ]
        q_df["region_order"] = q_df["region"].map(
            {r: i for i, r in enumerate(region_order)}
        ).fillna(99)
        q_df = q_df.sort_values(["region_order", "region", "IID"]).reset_index(drop=True)
    else:
        q_df["region"] = q_df["FID"].apply(lambda x: "USER" if x == "USER" else "Unknown")

    colors = _ADMIX_PALETTES[0][:k]
    k_cols = [f"K{i+1}" for i in range(k)]

    fig, ax = plt.subplots(figsize=(max(14, len(q_df) * 0.04), 5))

    bottom = np.zeros(len(q_df))
    for i, col in enumerate(k_cols):
        ax.bar(
            range(len(q_df)),
            q_df[col].values,
            bottom=bottom,
            color=colors[i],
            width=1.0,
            linewidth=0,
            label=f"祖源成分 {i+1}",
        )
        bottom += q_df[col].values

    # 标记用户样本
    user_indices = q_df[q_df["FID"] == "USER"].index.tolist()
    for idx in user_indices:
        ax.axvline(x=idx, color="white", lw=1.5, alpha=0.9)
        ax.text(
            idx, 1.02, "★你", ha="center", va="bottom",
            fontsize=7, color="white", fontweight="bold",
        )

    # 地理区域分割线 & 标签
    if "region" in q_df.columns:
        region_bounds: list[tuple[int, int, str]] = []
        prev_region, start = q_df.iloc[0]["region"], 0
        for i, row in q_df.iterrows():
            if row["region"] != prev_region:
                region_bounds.append((start, i - 1, prev_region))
                start = i
                prev_region = row["region"]
        region_bounds.append((start, len(q_df) - 1, prev_region))

        for s, e, region in region_bounds:
            ax.axvline(x=s - 0.5, color="#30363d", lw=1.2)
            mid = (s + e) / 2
            ax.text(
                mid, -0.08, region.replace("_", "\n"),
                ha="center", va="top", fontsize=6.5, color="#8b949e",
                transform=ax.get_xaxis_transform(), clip_on=False,
            )

    ax.set_xlim(-0.5, len(q_df) - 0.5)
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["0%", "25%", "50%", "75%", "100%"])
    ax.set_xticks([])
    ax.set_title(f"ADMIXTURE 无监督聚类（K={k}）— {len(q_df)} 个样本")
    ax.set_ylabel("祖源比例")
    ax.legend(
        loc="upper left", fontsize=8, ncol=min(k, 5),
        framealpha=0.8, bbox_to_anchor=(0, 1.15),
    )

    fig.tight_layout()
    out_path = str(Path(out_dir) / f"admixture_K{k}.png")
    _save(fig, out_path)


# ──────────────────────────────────────────────────────────────────────────────
# CV error 曲线
# ──────────────────────────────────────────────────────────────────────────────

def _plot_cv_error(admix_results: dict[int, dict], out_dir: str) -> None:
    """绘制各 K 值 CV error 曲线。"""
    ks       = sorted(admix_results.keys())
    cv_errs  = [admix_results[k].get("cv_error") for k in ks]

    if all(v is None for v in cv_errs):
        console.print("  [dim]没有有效的 CV error 数据，跳过曲线绘制[/dim]")
        return

    fig, ax = plt.subplots(figsize=(7, 4))
    valid_ks  = [k for k, v in zip(ks, cv_errs) if v is not None]
    valid_cvs = [v for v in cv_errs if v is not None]

    ax.plot(valid_ks, valid_cvs, "o-", color="#4CC9F0", lw=2, ms=8, zorder=3)

    best_k = valid_ks[int(np.argmin(valid_cvs))]
    best_v = min(valid_cvs)
    ax.scatter([best_k], [best_v], c="#F72585", s=120, zorder=4)
    ax.annotate(
        f"K={best_k}\n最优",
        xy=(best_k, best_v),
        xytext=(best_k + 0.3, best_v + (max(valid_cvs) - min(valid_cvs)) * 0.15),
        color="#F72585", fontsize=9,
    )

    ax.set_xlabel("K（聚类数）")
    ax.set_ylabel("CV Error（交叉验证误差）")
    ax.set_title("ADMIXTURE Cross-Validation Error by K")
    ax.set_xticks(valid_ks)
    ax.grid(True, ls="--", lw=0.5)
    fig.tight_layout()

    out_path = str(Path(out_dir) / "cv_error.png")
    _save(fig, out_path)


# ──────────────────────────────────────────────────────────────────────────────
# 公开入口
# ──────────────────────────────────────────────────────────────────────────────

def make_all_plots(
    pca_results: dict,
    admix_results: dict[int, dict],
    ref_dir: str,
    user_vcf: Optional[str],
    out_dir: str,
) -> None:
    """
    生成所有图表。

    Args:
        pca_results:   来自 pca.run_pca() 的结果字典
        admix_results: 来自 admixture.run_admixture() 的结果字典（键为 K 值）
        ref_dir:       HGDP 参考面板目录（含标签文件）
        user_vcf:      用户 VCF 路径（可为 None，仅用于图表标题）
        out_dir:       输出目录
    """
    os.makedirs(out_dir, exist_ok=True)

    console.print("\n  [bold]绘图：PCA (PC1 vs PC2)[/bold]")
    _plot_pca(pca_results, ref_dir=ref_dir, out_dir=out_dir, pc_x=1, pc_y=2)

    console.print("  [bold]绘图：PCA (PC3 vs PC4)[/bold]")
    _plot_pca(pca_results, ref_dir=ref_dir, out_dir=out_dir, pc_x=3, pc_y=4)

    # 找到合并后的 FAM 文件
    merged_fam = None
    if admix_results:
        first = next(iter(admix_results.values()))
        q_file = first["q_file"]
        # FAM 文件与 BED 同前缀
        bed_prefix = str(Path(q_file).parent / Path(q_file).stem.rsplit(".", 2)[0])
        fam_candidate = bed_prefix + ".fam"
        # 查找 work 目录下的 merged.fam
        fam_candidate2 = str(Path(q_file).parent.parent / "work" / "merged.fam")
        for fc in [fam_candidate, fam_candidate2]:
            if Path(fc).exists():
                merged_fam = fc
                break

    for k, result in sorted(admix_results.items()):
        console.print(f"  [bold]绘图：ADMIXTURE K={k}[/bold]")
        if merged_fam:
            _plot_admixture(
                k=k,
                q_file=result["q_file"],
                fam_file=merged_fam,
                ref_dir=ref_dir,
                out_dir=out_dir,
            )
        else:
            console.print(f"  [yellow]未找到 merged.fam，跳过 K={k} 条形图[/yellow]")

    if len(admix_results) > 1:
        console.print("  [bold]绘图：CV Error 曲线[/bold]")
        _plot_cv_error(admix_results, out_dir=out_dir)

    console.print(f"\n  [bold green]✓ 所有图表已保存至：{os.path.abspath(out_dir)}[/bold green]")
