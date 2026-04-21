"""
qc.py — VCF quality control & LD pruning

Pipeline:
  1. VCF → PLINK BED  (plink --vcf --make-bed)
  2. QC filtering      (--geno / --maf / --hwe)
  3. LD pruning        (--indep-pairwise 50 10 0.2)
  4. Extract pruned SNPs (--extract .prune.in)

Returns the LD-pruned BED prefix (without the .bed extension).
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

from rich.console import Console

console = Console()


# ──────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ──────────────────────────────────────────────────────────────────────────────

def _run_plink(args: list[str], step: str) -> None:
    """Run a PLINK command; raise RuntimeError on failure.

    PLINK exit codes:
      0 — success
      3 — warnings issued but output still written (treated as success)
      other — genuine failure
    """
    cmd = ["plink"] + args
    console.print(f"  [dim]$ {' '.join(cmd)}[/dim]")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode not in (0, 3):
        console.print(f"[red]PLINK failed at '{step}' (exit {result.returncode}):[/red]")
        console.print(result.stderr or result.stdout)
        raise RuntimeError(f"PLINK step failed: {step}")


def _count_variants(bed_prefix: str) -> int:
    """Count SNPs from the corresponding .bim file."""
    bim = Path(bed_prefix + ".bim")
    if not bim.exists():
        return 0
    with open(bim) as f:
        return sum(1 for _ in f)


def _count_samples(bed_prefix: str) -> int:
    """Count samples from the corresponding .fam file."""
    fam = Path(bed_prefix + ".fam")
    if not fam.exists():
        return 0
    with open(fam) as f:
        return sum(1 for _ in f)


# ──────────────────────────────────────────────────────────────────────────────
# Public entry point
# ──────────────────────────────────────────────────────────────────────────────

def run_qc(
    vcf_path: str,
    out_dir: str,
    geno: float = 0.05,
    maf: float = 0.01,
    hwe: float = 1e-6,
    threads: int = 4,
) -> str:
    """
    Run QC and LD pruning on an input VCF.

    Args:
        vcf_path: Path to the input VCF file
        out_dir:  Working directory for intermediate files
        geno:     Genotype missingness threshold (0–1)
        maf:      Minimum allele frequency
        hwe:      Hardy-Weinberg equilibrium p-value cutoff
        threads:  Number of PLINK threads

    Returns:
        Path prefix of the LD-pruned BED dataset (without .bed extension)
    """
    os.makedirs(out_dir, exist_ok=True)
    prefix = str(Path(out_dir) / "user")

    # ── Step 1a: VCF → BED ───────────────────────────────────────────────────
    console.print("  [bold]1a[/bold] VCF → PLINK BED")
    _run_plink([
        "--vcf", vcf_path,
        "--make-bed",
        "--out", prefix + "_raw",
        "--allow-extra-chr",   # allow non-standard chromosome names (e.g. 23andMe)
        "--double-id",         # use VCF sample ID as both FID and IID
        "--threads", str(threads),
    ], step="VCF → BED")

    n_snp_raw = _count_variants(prefix + "_raw")
    n_sam_raw = _count_samples(prefix + "_raw")
    console.print(f"    Raw: {n_sam_raw} samples, {n_snp_raw:,} SNPs")

    # ── Step 1b: QC filtering ─────────────────────────────────────────────────
    console.print(f"  [bold]1b[/bold] QC filtering (geno={geno}, maf={maf}, hwe={hwe:.0e})")
    _run_plink([
        "--bfile", prefix + "_raw",
        "--geno", str(geno),
        "--maf",  str(maf),
        "--hwe",  str(hwe),
        "--autosome",          # autosomes only (chr 1–22)
        "--snps-only",         # exclude INDELs
        "--make-bed",
        "--out", prefix + "_qc",
        "--threads", str(threads),
    ], step="QC filtering")

    n_snp_qc = _count_variants(prefix + "_qc")
    console.print(f"    After QC: {n_snp_qc:,} SNPs retained")
    if n_snp_qc == 0:
        raise RuntimeError(
            "No SNPs survived QC filtering. Possible causes:\n"
            "  • The VCF is empty or uses unsupported format\n"
            "  • The geno/maf/hwe thresholds are too strict\n"
            "  • --autosome filtering removed all variants (non-standard chr names)\n"
            "Try loosening the thresholds with --geno 0.1 --maf 0.005 --hwe 1e-10"
        )

    # ── Step 1c: LD pruning ───────────────────────────────────────────────────
    console.print("  [bold]1c[/bold] LD pruning (window=50, step=10, r²<0.2)")
    _run_plink([
        "--bfile", prefix + "_qc",
        "--indep-pairwise", "50", "10", "0.2",
        "--out", prefix + "_ld",
        "--threads", str(threads),
    ], step="LD pruning")

    prune_in = prefix + "_ld.prune.in"
    if not Path(prune_in).exists():
        # Read PLINK's log for useful diagnostics
        log_path = prefix + "_ld.log"
        log_tail = ""
        if Path(log_path).exists():
            with open(log_path) as lf:
                lines = lf.readlines()
            log_tail = "\n  Last lines of PLINK log:\n" + "".join(
                "    " + l for l in lines[-20:]
            )
        raise FileNotFoundError(
            f"LD pruning output not found: {prune_in}\n"
            "This usually means PLINK's --indep-pairwise step produced no output.\n"
            "Possible causes:\n"
            "  • Too few SNPs passed QC (check geno/maf/hwe thresholds)\n"
            "  • PLINK exited with an unrecognised error code"
            + log_tail
        )
    with open(prune_in) as f:
        n_prune = sum(1 for _ in f)
    console.print(f"    After LD pruning: {n_prune:,} independent SNPs")

    # ── Step 1d: Extract pruned SNPs ─────────────────────────────────────────
    console.print("  [bold]1d[/bold] Extract pruned SNP set")
    _run_plink([
        "--bfile", prefix + "_qc",
        "--extract", prune_in,
        "--make-bed",
        "--out", prefix + "_pruned",
        "--threads", str(threads),
    ], step="extract pruned SNPs")

    final_prefix = prefix + "_pruned"
    console.print(
        f"  [green]✓[/green] QC complete: "
        f"{_count_samples(final_prefix)} samples, "
        f"{_count_variants(final_prefix):,} independent SNPs"
    )
    return final_prefix
