"""
merge.py — Merge user data with the HGDP reference panel

Pipeline:
  1. Find common SNPs between user BIM and HGDP BIM
  2. Extract common SNPs from both datasets
  3. plink --bmerge to combine
  4. Handle strand-flip SNPs if needed

The user sample FID is tagged as "USER" for easy identification in plots.
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
    if result.returncode not in (0, 3):  # PLINK exit 3 = warning but success
        console.print(f"[red]PLINK failed at '{step}' (exit {result.returncode}):[/red]")
        console.print(result.stderr or result.stdout)
        raise RuntimeError(f"PLINK step failed: {step}")


def _load_bim(bim_path: str) -> pd.DataFrame:
    """Load a .bim file into a DataFrame."""
    return pd.read_csv(
        bim_path, sep="\t", header=None,
        names=["chrom", "snp", "cm", "pos", "a1", "a2"],
    )


def _common_snps(user_bim: str, ref_bim: str) -> set[str]:
    """Return the set of SNP IDs shared by both BIM files."""
    user_df = _load_bim(user_bim)
    ref_df  = _load_bim(ref_bim)
    return set(user_df["snp"]) & set(ref_df["snp"])


def merge_with_hgdp(
    user_bed: str,
    ref_dir: str,
    out_dir: str,
) -> str:
    """
    Merge the user's LD-pruned dataset with the HGDP reference panel.

    Args:
        user_bed: User BED prefix (returned by qc.py)
        ref_dir:  HGDP reference panel directory (contains hgdp_pruned.bed/.bim/.fam)
        out_dir:  Working directory

    Returns:
        Path prefix of the merged BED dataset
    """
    os.makedirs(out_dir, exist_ok=True)

    ref_prefix    = str(Path(ref_dir) / "hgdp_pruned")
    merged_prefix = str(Path(out_dir) / "merged")

    # Verify reference panel files exist
    for ext in (".bed", ".bim", ".fam"):
        if not Path(ref_prefix + ext).exists():
            raise FileNotFoundError(
                f"HGDP reference panel file not found: {ref_prefix + ext}\n"
                "Run 'dna download' first."
            )

    # ── Step 2a: Find common SNPs ─────────────────────────────────────────────
    console.print("  [bold]2a[/bold] Computing common SNPs")
    common = _common_snps(user_bed + ".bim", ref_prefix + ".bim")
    console.print(f"    User × HGDP common SNPs: {len(common):,}")

    if len(common) < 1000:
        console.print(
            f"  [yellow]Warning: only {len(common)} common SNPs "
            "(< 1000). Results may be unreliable. "
            "Check that both datasets use the same genome build (hg19 vs hg38).[/yellow]"
        )

    common_list_path = str(Path(out_dir) / "common_snps.txt")
    with open(common_list_path, "w") as f:
        f.write("\n".join(sorted(common)))

    # ── Step 2b: Extract common SNPs ─────────────────────────────────────────
    console.print("  [bold]2b[/bold] Extracting common SNPs from both datasets")

    user_common = str(Path(out_dir) / "user_common")
    _run_plink([
        "--bfile", user_bed,
        "--extract", common_list_path,
        "--make-bed", "--out", user_common,
    ], step="extract user common SNPs")

    ref_common = str(Path(out_dir) / "ref_common")
    _run_plink([
        "--bfile", ref_prefix,
        "--extract", common_list_path,
        "--make-bed", "--out", ref_common,
    ], step="extract HGDP common SNPs")

    # ── Step 2c: Tag user sample ──────────────────────────────────────────────
    console.print("  [bold]2c[/bold] Tagging user sample (FID → USER)")
    fam_path = user_common + ".fam"
    fam_df = pd.read_csv(
        fam_path, sep=r"\s+", header=None,
        names=["fid", "iid", "pat", "mat", "sex", "phen"],
    )
    fam_df["fid"] = "USER"
    fam_df.to_csv(fam_path, sep="\t", header=False, index=False)

    # ── Step 2d: bmerge ───────────────────────────────────────────────────────
    console.print("  [bold]2d[/bold] bmerge: user data + HGDP")
    _run_plink([
        "--bfile", ref_common,
        "--bmerge", user_common + ".bed", user_common + ".bim", user_common + ".fam",
        "--make-bed", "--out", merged_prefix,
        "--allow-no-sex",
    ], step="bmerge")

    # Handle strand flips
    missnp_path = merged_prefix + "-merge.missnp"
    if Path(missnp_path).exists():
        console.print("  [yellow]Strand-flip SNPs detected — correcting...[/yellow]")

        user_flip = str(Path(out_dir) / "user_flipped")
        _run_plink([
            "--bfile", user_common,
            "--flip", missnp_path,
            "--make-bed", "--out", user_flip,
        ], step="strand flip")

        _run_plink([
            "--bfile", ref_common,
            "--bmerge", user_flip + ".bed", user_flip + ".bim", user_flip + ".fam",
            "--make-bed", "--out", merged_prefix,
            "--allow-no-sex",
        ], step="re-merge after flip")

    from dna.qc import _count_samples, _count_variants
    n_sam = _count_samples(merged_prefix)
    n_snp = _count_variants(merged_prefix)
    console.print(
        f"  [green]✓[/green] Merged: {n_sam} samples "
        f"({n_sam - 1} HGDP reference + 1 user), {n_snp:,} SNPs"
    )
    return merged_prefix
