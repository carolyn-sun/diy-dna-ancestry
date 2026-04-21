"""
cli.py — diy-dna-ancestry 主命令入口

子命令：
  dna init      环境检查
  dna download  下载 HGDP 参考面板
  dna run       完整分析流水线
  dna plot      单独绘图
"""

import sys
import click
from rich.console import Console
from rich import print as rprint

console = Console()


# ──────────────────────────────────────────────────────────────────────────────
# 主组
# ──────────────────────────────────────────────────────────────────────────────

@click.group()
@click.version_option(package_name="diy-dna-ancestry")
def main():
    """🧬 diy-dna-ancestry：个人 DNA 祖源分析工具

    基于 HGDP 参考面板 + PLINK + ADMIXTURE + PCA。

    \b
    快速开始：
      dna init            # 检查环境
      dna download        # 下载参考面板
      dna run --vcf your.vcf --k 3,5
    """
    pass


# ──────────────────────────────────────────────────────────────────────────────
# dna init
# ──────────────────────────────────────────────────────────────────────────────

@main.command()
def init():
    """检查运行环境（Python、PLINK、ADMIXTURE）。"""
    from dna.init_env import run_check
    ok = run_check()
    sys.exit(0 if ok else 1)


# ──────────────────────────────────────────────────────────────────────────────
# dna download
# ──────────────────────────────────────────────────────────────────────────────

@main.command()
@click.option(
    "--out-dir", default="data/hgdp", show_default=True,
    help="参考面板下载目录"
)
@click.option(
    "--force", is_flag=True, default=False,
    help="强制重新下载（即使文件已存在）"
)
def download(out_dir: str, force: bool):
    """下载 HGDP 参考面板（~50k SNP LD 剪枝版本）。"""
    from dna.download import download_hgdp
    download_hgdp(out_dir=out_dir, force=force)


# ──────────────────────────────────────────────────────────────────────────────
# dna run
# ──────────────────────────────────────────────────────────────────────────────

@main.command()
@click.option(
    "--vcf", required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="输入 VCF 文件路径"
)
@click.option(
    "--k", "k_values", default="3,5", show_default=True,
    help="ADMIXTURE K 值列表，逗号分隔（如 3,5 或 2,3,4,5）"
)
@click.option(
    "--threads", default=4, show_default=True,
    help="并行线程数"
)
@click.option(
    "--out", "out_dir", default="results", show_default=True,
    help="输出目录"
)
@click.option(
    "--ref-dir", default="data/hgdp", show_default=True,
    help="HGDP 参考面板目录（由 dna download 创建）"
)
@click.option("--geno", default=0.05, show_default=True, help="基因型缺失率阈值")
@click.option("--maf",  default=0.01, show_default=True, help="最小等位基因频率")
@click.option("--hwe",  default=1e-6, show_default=True, help="Hardy-Weinberg p 值阈值")
@click.option(
    "--skip-plot", is_flag=True, default=False,
    help="跳过绘图步骤"
)
def run(vcf: str, k_values: str, threads: int, out_dir: str, ref_dir: str,
        geno: float, maf: float, hwe: float, skip_plot: bool):
    """完整分析流水线：QC → 合并 → ADMIXTURE → PCA → 绘图。"""
    from dna.init_env import run_check
    from dna.qc import run_qc
    from dna.merge import merge_with_hgdp
    from dna.admixture import run_admixture
    from dna.pca import run_pca
    from dna.plot import make_all_plots
    import os

    # 解析 K 值
    try:
        ks = [int(k.strip()) for k in k_values.split(",")]
    except ValueError:
        console.print("[red]--k 参数格式错误，请使用逗号分隔的整数，如 3,5[/red]")
        sys.exit(1)

    os.makedirs(out_dir, exist_ok=True)
    work_dir = os.path.join(out_dir, "work")
    os.makedirs(work_dir, exist_ok=True)

    console.rule("[bold cyan]diy-dna-ancestry 祖源分析流水线[/bold cyan]")

    # 0. 环境检查
    console.print("\n[bold]0 / 5  环境检查[/bold]")
    if not run_check(verbose=False):
        console.print("[red]环境检查失败，请先运行 dna init 排查问题[/red]")
        sys.exit(1)

    # 1. QC & LD 剪枝
    console.print("\n[bold]1 / 5  VCF 质控 & LD 剪枝[/bold]")
    user_bed = run_qc(
        vcf_path=vcf,
        out_dir=work_dir,
        geno=geno, maf=maf, hwe=hwe,
        threads=threads,
    )

    # 2. 合并 HGDP
    console.print("\n[bold]2 / 5  与 HGDP 参考面板合并[/bold]")
    merged_bed = merge_with_hgdp(
        user_bed=user_bed,
        ref_dir=ref_dir,
        out_dir=work_dir,
    )

    # 3. ADMIXTURE
    console.print(f"\n[bold]3 / 5  ADMIXTURE（K = {ks}）[/bold]")
    admix_results = run_admixture(
        bed=merged_bed,
        ks=ks,
        out_dir=os.path.join(out_dir, "admixture"),
        threads=threads,
    )

    # 4. PCA
    console.print("\n[bold]4 / 5  PCA 主成分分析[/bold]")
    pca_results = run_pca(
        bed=merged_bed,
        out_dir=os.path.join(out_dir, "pca"),
        threads=threads,
    )

    # 5. 绘图
    if not skip_plot:
        console.print("\n[bold]5 / 5  生成可视化图表[/bold]")
        make_all_plots(
            pca_results=pca_results,
            admix_results=admix_results,
            ref_dir=ref_dir,
            user_vcf=vcf,
            out_dir=out_dir,
        )
    else:
        console.print("\n[dim]5 / 5  绘图步骤已跳过（--skip-plot）[/dim]")

    console.rule("[bold green]分析完成[/bold green]")
    console.print(f"结果保存在：[cyan]{os.path.abspath(out_dir)}[/cyan]")


# ──────────────────────────────────────────────────────────────────────────────
# dna plot
# ──────────────────────────────────────────────────────────────────────────────

@main.command()
@click.option(
    "--results", "results_dir", default="results", show_default=True,
    type=click.Path(exists=True, file_okay=False),
    help="已有结果目录（由 dna run 生成）"
)
@click.option(
    "--ref-dir", default="data/hgdp", show_default=True,
    help="HGDP 参考面板目录"
)
def plot(results_dir: str, ref_dir: str):
    """从已有结果重新绘图（无需重跑分析）。"""
    from dna.plot import make_all_plots
    import os, glob

    # 自动发现 admixture 结果
    q_files = glob.glob(os.path.join(results_dir, "admixture", "*.Q"))
    admix_results = {}
    for qf in q_files:
        stem = os.path.basename(qf)  # e.g. merged.3.Q
        try:
            k = int(stem.split(".")[-2])
            admix_results[k] = qf
        except (ValueError, IndexError):
            pass

    pca_eigenvec = os.path.join(results_dir, "pca", "pca.eigenvec")
    pca_eigenval = os.path.join(results_dir, "pca", "pca.eigenval")
    if not os.path.exists(pca_eigenvec):
        console.print(f"[red]找不到 PCA 结果：{pca_eigenvec}[/red]")
        sys.exit(1)

    make_all_plots(
        pca_results={"eigenvec": pca_eigenvec, "eigenval": pca_eigenval},
        admix_results=admix_results,
        ref_dir=ref_dir,
        user_vcf=None,
        out_dir=results_dir,
    )
    console.print(f"[green]绘图完成，保存在：{os.path.abspath(results_dir)}[/green]")


if __name__ == "__main__":
    main()
