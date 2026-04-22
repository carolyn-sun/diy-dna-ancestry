"""
cli.py — diy-dna-ancestry main command entry point

Subcommands:
  dna init      Environment check
  dna download  Download HGDP reference panel
  dna run       Full analysis pipeline
  dna plot      Re-plot from existing results
"""

import sys
import click
from rich.console import Console

console = Console()


# ──────────────────────────────────────────────────────────────────────────────
# Main group
# ──────────────────────────────────────────────────────────────────────────────

@click.group()
@click.version_option(package_name="diy-dna-ancestry")
def main():
    """🧬 diy-dna-ancestry: personal DNA ancestry analysis tool

    Powered by the HGDP reference panel, PLINK, ADMIXTURE, and PCA.

    \b
    Quick start:
      dna init            # check environment
      dna download        # download reference panel
      dna run --vcf your.vcf --k 3,5
    """
    pass


# ──────────────────────────────────────────────────────────────────────────────
# dna init
# ──────────────────────────────────────────────────────────────────────────────

@main.command()
def init():
    """Check runtime environment (Python, PLINK, ADMIXTURE)."""
    from dna.init_env import run_check
    ok = run_check()
    sys.exit(0 if ok else 1)


# ──────────────────────────────────────────────────────────────────────────────
# dna download
# ──────────────────────────────────────────────────────────────────────────────

@main.command()
@click.option(
    "--out-dir", default="data/hgdp", show_default=True,
    help="Directory to download the reference panel into"
)
@click.option(
    "--force", is_flag=True, default=False,
    help="Force re-download even if files already exist"
)
def download(out_dir: str, force: bool):
    """Download the HGDP reference panel (~50k LD-pruned SNPs)."""
    from dna.download import download_hgdp
    download_hgdp(out_dir=out_dir, force=force)


# ──────────────────────────────────────────────────────────────────────────────
# dna run
# ──────────────────────────────────────────────────────────────────────────────

@main.command()
@click.option(
    "--vcf", required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Input VCF file path"
)
@click.option(
    "--k", "k_values", default="3,5", show_default=True,
    help="ADMIXTURE K values, comma-separated (e.g. 3,5 or 2,3,4,5)"
)
@click.option("--threads", default=4, show_default=True, help="Number of parallel threads")
@click.option("--out", "out_dir", default="results", show_default=True, help="Output directory")
@click.option("--ref-dir", default="data/hgdp", show_default=True, help="HGDP reference panel directory")
@click.option("--geno", default=0.05, show_default=True, help="Genotype missingness threshold")
@click.option("--maf",  default=0.01, show_default=True, help="Minimum allele frequency")
@click.option("--hwe",  default=1e-6, show_default=True, help="Hardy-Weinberg p-value cutoff")
@click.option("--skip-plot", is_flag=True, default=False, help="Skip the plotting step")
@click.option(
    "--nmf-fallback", "nmf_fallback", is_flag=True, default=False,
    help=(
        "Enable Python NMF as a fallback when the ADMIXTURE binary crashes "
        "(e.g. SIGSEGV on incompatible hardware). "
        "Results are approximate — prefer native ADMIXTURE when possible."
    ),
)
@click.option(
    "--admixture-bin", "admixture_bin", default="admixture", show_default=True,
    help=(
        "Path (or name) of the ADMIXTURE executable. "
        "Use this to point at a 32-bit build on WSL, e.g. "
        "--admixture-bin ~/bin/admixture32."
    ),
)
def run(vcf: str, k_values: str, threads: int, out_dir: str, ref_dir: str,
        geno: float, maf: float, hwe: float, skip_plot: bool,
        nmf_fallback: bool, admixture_bin: str):
    """Full pipeline: QC → merge → ADMIXTURE → PCA → plots."""
    from dna.init_env import run_check
    from dna.qc import run_qc
    from dna.merge import merge_with_hgdp
    from dna.admixture import run_admixture
    from dna.pca import run_pca
    from dna.plot import make_all_plots
    import os

    try:
        ks = [int(k.strip()) for k in k_values.split(",")]
    except ValueError:
        console.print("[red]--k must be comma-separated integers, e.g. 3,5[/red]")
        sys.exit(1)

    os.makedirs(out_dir, exist_ok=True)
    work_dir = os.path.join(out_dir, "work")
    os.makedirs(work_dir, exist_ok=True)

    console.rule("[bold cyan]diy-dna-ancestry pipeline[/bold cyan]")

    console.print("\n[bold]0 / 5  Environment check[/bold]")
    if not run_check(verbose=False):
        console.print("[red]Environment check failed — run 'dna init' to diagnose[/red]")
        sys.exit(1)

    console.print("\n[bold]1 / 5  VCF QC & LD pruning[/bold]")
    user_bed = run_qc(
        vcf_path=vcf, out_dir=work_dir,
        geno=geno, maf=maf, hwe=hwe, threads=threads,
    )

    console.print("\n[bold]2 / 5  Merge with HGDP reference panel[/bold]")
    merged_bed = merge_with_hgdp(user_bed=user_bed, ref_dir=ref_dir, out_dir=work_dir)

    console.print(f"\n[bold]3 / 5  ADMIXTURE (K = {ks})[/bold]")
    if nmf_fallback:
        console.print("  [yellow]NMF fallback mode enabled (--nmf-fallback)[/yellow]")
    if admixture_bin != "admixture":
        console.print(f"  [dim]ADMIXTURE binary: {admixture_bin}[/dim]")
    admix_results = run_admixture(
        bed=merged_bed, ks=ks,
        out_dir=os.path.join(out_dir, "admixture"),
        threads=threads,
        allow_nmf_fallback=nmf_fallback,
        admixture_bin=admixture_bin,
    )

    console.print("\n[bold]4 / 5  PCA[/bold]")
    pca_results = run_pca(
        bed=merged_bed,
        out_dir=os.path.join(out_dir, "pca"),
        threads=threads,
    )

    if not skip_plot:
        console.print("\n[bold]5 / 5  Generate plots[/bold]")
        make_all_plots(
            pca_results=pca_results, admix_results=admix_results,
            ref_dir=ref_dir, user_vcf=vcf, out_dir=out_dir,
        )
    else:
        console.print("\n[dim]5 / 5  Plotting skipped (--skip-plot)[/dim]")

    console.rule("[bold green]Done[/bold green]")
    console.print(f"Results saved to: [cyan]{os.path.abspath(out_dir)}[/cyan]")


# ──────────────────────────────────────────────────────────────────────────────
# dna plot
# ──────────────────────────────────────────────────────────────────────────────

@main.command()
@click.option(
    "--results", "results_dir", default="results", show_default=True,
    type=click.Path(exists=True, file_okay=False),
    help="Results directory produced by 'dna run'"
)
@click.option("--ref-dir", default="data/hgdp", show_default=True, help="HGDP reference panel directory")
def plot(results_dir: str, ref_dir: str):
    """Re-plot from existing results (no need to re-run the analysis)."""
    from dna.plot import make_all_plots
    import os, glob

    q_files = glob.glob(os.path.join(results_dir, "admixture", "*.Q"))
    admix_results = {}
    for qf in q_files:
        stem = os.path.basename(qf)
        try:
            k = int(stem.split(".")[-2])
            admix_results[k] = qf
        except (ValueError, IndexError):
            pass

    pca_eigenvec = os.path.join(results_dir, "pca", "pca.eigenvec")
    if not os.path.exists(pca_eigenvec):
        console.print(f"[red]PCA results not found: {pca_eigenvec}[/red]")
        sys.exit(1)

    make_all_plots(
        pca_results={"eigenvec": pca_eigenvec, "eigenval": pca_eigenvec.replace(".eigenvec", ".eigenval")},
        admix_results=admix_results,
        ref_dir=ref_dir,
        user_vcf=None,
        out_dir=results_dir,
    )
    console.print(f"[green]Plots saved to: {os.path.abspath(results_dir)}[/green]")


if __name__ == "__main__":
    main()
