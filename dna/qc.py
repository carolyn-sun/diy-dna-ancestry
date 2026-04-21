"""
qc.py — VCF 质控 & LD 剪枝模块

流程：
  1. VCF → PLINK BED（plink --vcf --make-bed）
  2. 质控过滤（--geno / --maf / --hwe）
  3. LD 剪枝（--indep-pairwise 50 10 0.2）
  4. 提取独立 SNPs（--extract .prune.in）

返回：经过剪枝的 BED 前缀路径（不含 .bed 后缀）
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

from rich.console import Console

console = Console()


# ──────────────────────────────────────────────────────────────────────────────
# 内部工具函数
# ──────────────────────────────────────────────────────────────────────────────

def _run_plink(args: list[str], step: str) -> None:
    """执行 PLINK 命令，失败则抛出 RuntimeError。"""
    cmd = ["plink"] + args
    console.print(f"  [dim]$ {' '.join(cmd)}[/dim]")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        console.print(f"[red]PLINK {step} 失败（退出码 {result.returncode}）：[/red]")
        console.print(result.stderr or result.stdout)
        raise RuntimeError(f"PLINK 步骤失败：{step}")


def _count_variants(bed_prefix: str) -> int:
    """从 .bim 文件统计 SNP 数量。"""
    bim = Path(bed_prefix + ".bim")
    if not bim.exists():
        return 0
    with open(bim) as f:
        return sum(1 for _ in f)


def _count_samples(bed_prefix: str) -> int:
    """从 .fam 文件统计样本数量。"""
    fam = Path(bed_prefix + ".fam")
    if not fam.exists():
        return 0
    with open(fam) as f:
        return sum(1 for _ in f)


# ──────────────────────────────────────────────────────────────────────────────
# 公开入口
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
    对输入 VCF 执行质控和 LD 剪枝。

    Args:
        vcf_path: 输入 VCF 文件路径
        out_dir:  工作目录（中间文件存放处）
        geno:     基因型缺失率阈值（0~1）
        maf:      最小等位基因频率
        hwe:      Hardy-Weinberg 均衡检验 p 值阈值
        threads:  PLINK 使用的线程数

    Returns:
        LD 剪枝后的 BED 前缀路径字符串（不含 .bed 后缀）
    """
    os.makedirs(out_dir, exist_ok=True)
    prefix = str(Path(out_dir) / "user")

    # ── 步骤 1：VCF → BED ────────────────────────────────────────────────────
    console.print("  [bold]1a[/bold] VCF → PLINK BED 转换")
    _run_plink([
        "--vcf", vcf_path,
        "--make-bed",
        "--out", prefix + "_raw",
        "--allow-extra-chr",          # 允许非标准染色体名（如 23andMe 格式）
        "--double-id",                # 使用 VCF 样本 ID 作为 FID/IID
        "--threads", str(threads),
    ], step="VCF → BED")

    n_snp_raw = _count_variants(prefix + "_raw")
    n_sam_raw = _count_samples(prefix + "_raw")
    console.print(f"    原始数据：{n_sam_raw} 个样本，{n_snp_raw:,} 个 SNPs")

    # ── 步骤 2：质控过滤 ──────────────────────────────────────────────────────
    console.print(f"  [bold]1b[/bold] 质控过滤（geno={geno}, maf={maf}, hwe={hwe:.0e}）")
    _run_plink([
        "--bfile", prefix + "_raw",
        "--geno", str(geno),
        "--maf",  str(maf),
        "--hwe",  str(hwe),
        "--autosome",                 # 仅保留常染色体（1–22）
        "--snps-only",                # 排除 INDEL
        "--make-bed",
        "--out", prefix + "_qc",
        "--threads", str(threads),
    ], step="QC 过滤")

    n_snp_qc = _count_variants(prefix + "_qc")
    console.print(f"    质控后保留：{n_snp_qc:,} 个 SNPs")

    # ── 步骤 3：LD 剪枝（生成 .prune.in 文件）────────────────────────────────
    console.print("  [bold]1c[/bold] LD 剪枝（窗口 50 SNPs，步长 10，r² < 0.2）")
    _run_plink([
        "--bfile", prefix + "_qc",
        "--indep-pairwise", "50", "10", "0.2",
        "--out", prefix + "_ld",
        "--threads", str(threads),
    ], step="LD 剪枝")

    prune_in = prefix + "_ld.prune.in"
    with open(prune_in) as f:
        n_prune = sum(1 for _ in f)
    console.print(f"    LD 剪枝后保留：{n_prune:,} 个独立 SNPs")

    # ── 步骤 4：提取独立 SNPs ─────────────────────────────────────────────────
    console.print("  [bold]1d[/bold] 提取独立 SNPs")
    _run_plink([
        "--bfile", prefix + "_qc",
        "--extract", prune_in,
        "--make-bed",
        "--out", prefix + "_pruned",
        "--threads", str(threads),
    ], step="提取独立 SNPs")

    final_prefix = prefix + "_pruned"
    console.print(
        f"  [green]✓[/green] 用户数据质控完成："
        f"{_count_samples(final_prefix)} 个样本，"
        f"{_count_variants(final_prefix):,} 个独立 SNPs"
    )
    return final_prefix
