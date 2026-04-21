"""
pca.py — PLINK PCA (principal component analysis)

Command:
  plink --bfile merged --pca 10 --out results/pca/pca

Outputs:
  pca.eigenvec  — per-sample PC coordinates (PC1–PC10)
  pca.eigenval  — eigenvalues (variance explained per PC)
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
        console.print(f"[red]PLINK failed at '{step}' (exit {result.returncode}):[/red]")
        console.print(result.stderr or result.stdout)
        raise RuntimeError(f"PLINK step failed: {step}")


def run_pca(
    bed: str,
    out_dir: str,
    n_pcs: int = 10,
    threads: int = 4,
) -> dict:
    """
    Run PLINK PCA on a merged dataset.

    Args:
        bed:     Merged BED prefix
        out_dir: Output directory
        n_pcs:   Number of principal components to compute (default 10)
        threads: Number of threads

    Returns:
        Dict containing eigenvec/eigenval file paths and parsed DataFrames
    """
    os.makedirs(out_dir, exist_ok=True)
    out_prefix = str(Path(out_dir) / "pca")

    console.print(f"  [bold]Computing top {n_pcs} principal components[/bold]")
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
        raise FileNotFoundError(f"PCA output not found: {eigenvec_path}")

    # Parse .eigenvec
    pc_cols = [f"PC{i}" for i in range(1, n_pcs + 1)]
    eigenvec_df = pd.read_csv(
        eigenvec_path, sep=r"\s+", header=None,
        names=["FID", "IID"] + pc_cols,
    )

    # Parse .eigenval (variance explained)
    eigenval_df = pd.read_csv(eigenval_path, header=None, names=["eigenval"])
    total_var = eigenval_df["eigenval"].sum()
    eigenval_df["var_pct"] = eigenval_df["eigenval"] / total_var * 100

    # Print variance-explained summary
    console.print(f"\n  {'PC':<6} {'Var%':<12} {'Cumulative':<10}")
    console.print(f"  {'─' * 32}")
    cumsum = 0.0
    for i, row in eigenval_df.head(5).iterrows():
        cumsum += row["var_pct"]
        console.print(f"  PC{i+1:<4} {row['var_pct']:>7.2f}%     {cumsum:>7.2f}%")
    console.print(f"  {'─' * 32}")

    n_user = (eigenvec_df["FID"] == "USER").sum()
    console.print(
        f"\n  [green]✓[/green] PCA done: "
        f"{len(eigenvec_df)} samples ({n_user} user sample), "
        f"{n_pcs} components"
    )

    return {
        "eigenvec": eigenvec_path,
        "eigenval": eigenval_path,
        "eigenvec_df": eigenvec_df,
        "eigenval_df": eigenval_df,
    }
