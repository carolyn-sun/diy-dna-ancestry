"""
pca.py — PLINK PCA 主成分分析模块

命令：
  plink --bfile merged --pca 10 --out results/pca/pca

输出：
  pca.eigenvec   ← 各样本的主成分坐标（PC1~PC10）
  pca.eigenval   ← 各主成分的特征值（方差贡献）
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pandas as pd
from rich.console import Console

console = Console()


def _run_plink(args: list[str], step: str) -> None:
    cmd = ["plink"] + args
    console.print(f"  [dim]$ {' '.join(cmd)}[/dim]")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        console.print(f"[red]PLINK {step} 失败（退出码 {result.returncode}）：[/red]")
        console.print(result.stderr or result.stdout)
        raise RuntimeError(f"PLINK 步骤失败：{step}")


def run_pca(
    bed: str,
    out_dir: str,
    n_pcs: int = 10,
    threads: int = 4,
) -> dict:
    """
    跑 PLINK PCA。

    Args:
        bed:    合并后的 BED 前缀
        out_dir: 输出目录
        n_pcs:  计算的主成分数量（默认 10）
        threads: 线程数

    Returns:
        dict，含 eigenvec / eigenval 文件路径，和解析好的 DataFrame
    """
    os.makedirs(out_dir, exist_ok=True)
    out_prefix = str(Path(out_dir) / "pca")

    console.print(f"  [bold]计算前 {n_pcs} 个主成分[/bold]")
    _run_plink([
        "--bfile", bed,
        "--pca", str(n_pcs),
        "--out", out_prefix,
        "--threads", str(threads),
        "--allow-no-sex",
    ], step="PCA")

    eigenvec_path = out_prefix + ".eigenvec"
    eigenval_path = out_prefix + ".eigenval"

    if not Path(eigenvec_path).exists():
        raise FileNotFoundError(f"PCA 结果文件未生成：{eigenvec_path}")

    # 解析 .eigenvec
    pc_cols = [f"PC{i}" for i in range(1, n_pcs + 1)]
    eigenvec_df = pd.read_csv(
        eigenvec_path, sep=r"\s+", header=None,
        names=["FID", "IID"] + pc_cols,
    )

    # 解析 .eigenval（方差贡献率）
    eigenval_df = pd.read_csv(eigenval_path, header=None, names=["eigenval"])
    total_var = eigenval_df["eigenval"].sum()
    eigenval_df["var_pct"] = eigenval_df["eigenval"] / total_var * 100

    # 打印方差贡献表
    console.print(f"\n  {'PC':<6} {'方差贡献':<12} {'累计':<10}")
    console.print(f"  {'─'*30}")
    cumsum = 0.0
    for i, row in eigenval_df.head(5).iterrows():
        cumsum += row["var_pct"]
        console.print(
            f"  PC{i+1:<4} {row['var_pct']:>7.2f}%     {cumsum:>7.2f}%"
        )
    console.print(f"  {'─'*30}")

    n_user = (eigenvec_df["FID"] == "USER").sum()
    console.print(
        f"\n  [green]✓[/green] PCA 完成："
        f"{len(eigenvec_df)} 个样本（含 {n_user} 个用户样本），"
        f"{n_pcs} 个主成分"
    )

    return {
        "eigenvec": eigenvec_path,
        "eigenval": eigenval_path,
        "eigenvec_df": eigenvec_df,
        "eigenval_df": eigenval_df,
    }
