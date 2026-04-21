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


def _align_snp_ids(user_bim: str, ref_bim: str) -> set[str]:
    """Return common SNP IDs, with automatic fallback to chr:pos matching.

    Some VCFs use rsIDs while the reference panel uses 'chr:pos' IDs (or vice
    versa).  When direct name matching yields fewer than 100 hits we fall back
    to matching on (chromosome, base-pair position).  Matched user BIM entries
    are then *rewritten in place* to use the reference panel's SNP IDs so that
    downstream --extract calls work correctly.
    """
    user_df = _load_bim(user_bim)
    ref_df  = _load_bim(ref_bim)

    # ── Try direct SNP-ID match first ────────────────────────────────────────
    common_ids = set(user_df["snp"]) & set(ref_df["snp"])
    if len(common_ids) >= 100:
        return common_ids

    # ── Fallback: match by (chrom, pos) ──────────────────────────────────────
    console.print(
        f"    [yellow]Only {len(common_ids)} SNPs matched by ID; "
        "falling back to chromosome:position matching "
        "(SNP ID format mismatch between VCF and reference panel)[/yellow]"
    )

    # Build a (chrom, pos) → ref_snp_id lookup
    # Normalise chromosome labels: strip leading 'chr' for comparison
    ref_df["_chrom_norm"] = ref_df["chrom"].astype(str).str.lstrip("chr")
    ref_df["_key"] = ref_df["_chrom_norm"] + ":" + ref_df["pos"].astype(str)
    key_to_ref_snp: dict[str, str] = dict(zip(ref_df["_key"], ref_df["snp"]))

    user_df["_chrom_norm"] = user_df["chrom"].astype(str).str.lstrip("chr")
    user_df["_key"] = user_df["_chrom_norm"] + ":" + user_df["pos"].astype(str)
    user_df["_new_snp"] = user_df["_key"].map(key_to_ref_snp)

    matched = user_df[user_df["_new_snp"].notna()].copy()
    if matched.empty:
        return set()

    # Rewrite the user BIM so SNP IDs match the reference panel
    user_df.loc[matched.index, "snp"] = matched["_new_snp"]
    user_df[["chrom", "snp", "cm", "pos", "a1", "a2"]].to_csv(
        user_bim, sep="\t", header=False, index=False
    )
    console.print(
        f"    SNP IDs remapped via chr:pos — {len(matched):,} variants aligned"
    )
    return set(matched["_new_snp"])


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
    common = _align_snp_ids(user_bed + ".bim", ref_prefix + ".bim")
    console.print(f"    User × HGDP common SNPs: {len(common):,}")

    if len(common) == 0:
        raise RuntimeError(
            "No common SNPs found between your VCF and the HGDP reference panel.\n"
            "Possible causes:\n"
            "  • Genome build mismatch — your VCF may be hg38 while HGDP panel is hg19 (or vice versa)\n"
            "  • The VCF was not produced by a standard genotyping array\n"
            "Check your VCF header for 'reference' or 'genome_build' metadata."
        )
    if len(common) < 1000:
        console.print(
            f"  [yellow]Warning: only {len(common):,} common SNPs "
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

    # ── Step 2d: bmerge (up to 3 passes) ─────────────────────────────────────
    # Pass 1: try merging as-is
    # Pass 2: if missnp found, flip those SNPs in the user data and retry
    # Pass 3: if missnp still remain after flipping, exclude them and merge again
    console.print("  [bold]2d[/bold] bmerge: user data + HGDP")
    _run_plink([
        "--bfile", ref_common,
        "--bmerge", user_common + ".bed", user_common + ".bim", user_common + ".fam",
        "--make-bed", "--out", merged_prefix,
        "--allow-no-sex",
    ], step="bmerge pass 1")

    missnp_path = merged_prefix + "-merge.missnp"
    if Path(missnp_path).exists():
        n_missnp = sum(1 for _ in open(missnp_path))
        console.print(
            f"  [yellow]Strand-flip SNPs detected ({n_missnp:,}) — flipping and retrying...[/yellow]"
        )

        user_flip = str(Path(out_dir) / "user_flipped")
        _run_plink([
            "--bfile", user_common,
            "--flip", missnp_path,
            "--make-bed", "--out", user_flip,
            "--allow-no-sex",
        ], step="strand flip")

        # Remove stale missnp file before pass 2
        Path(missnp_path).unlink(missing_ok=True)

        _run_plink([
            "--bfile", ref_common,
            "--bmerge", user_flip + ".bed", user_flip + ".bim", user_flip + ".fam",
            "--make-bed", "--out", merged_prefix,
            "--allow-no-sex",
        ], step="bmerge pass 2")

        # If missnp still exists after flipping, those are ambiguous SNPs
        # (A/T or C/G strand-ambiguous) that cannot be resolved — exclude them
        if Path(missnp_path).exists():
            n_remaining = sum(1 for _ in open(missnp_path))
            console.print(
                f"  [yellow]{n_remaining:,} ambiguous SNPs still unresolved — excluding them[/yellow]"
            )

            user_clean = str(Path(out_dir) / "user_clean")
            _run_plink([
                "--bfile", user_flip,
                "--exclude", missnp_path,
                "--make-bed", "--out", user_clean,
                "--allow-no-sex",
            ], step="exclude ambiguous SNPs")

            ref_clean = str(Path(out_dir) / "ref_clean")
            _run_plink([
                "--bfile", ref_common,
                "--exclude", missnp_path,
                "--make-bed", "--out", ref_clean,
                "--allow-no-sex",
            ], step="exclude ambiguous SNPs from ref")

            Path(missnp_path).unlink(missing_ok=True)

            _run_plink([
                "--bfile", ref_clean,
                "--bmerge", user_clean + ".bed", user_clean + ".bim", user_clean + ".fam",
                "--make-bed", "--out", merged_prefix,
                "--allow-no-sex",
            ], step="bmerge pass 3")

    from dna.qc import _count_samples, _count_variants
    n_sam = _count_samples(merged_prefix)
    n_snp = _count_variants(merged_prefix)

    if n_sam == 0 or n_snp == 0:
        raise RuntimeError(
            f"Merge produced an empty dataset ({n_sam} samples, {n_snp} SNPs).\n"
            "This usually means all SNPs were excluded during strand-flip resolution.\n"
            "Possible causes:\n"
            "  • Genome build mismatch (hg19 vs hg38) — positions don't truly overlap\n"
            "  • The reference panel HGDP BIM uses a different allele coding\n"
            "Check results/work/merged.log for details."
        )

    console.print(
        f"  [green]✓[/green] Merged: {n_sam} samples "
        f"({n_sam - 1} HGDP reference + 1 user), {n_snp:,} SNPs"
    )
    return merged_prefix
