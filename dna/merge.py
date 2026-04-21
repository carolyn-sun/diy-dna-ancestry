"""
merge.py — 将用户数据与 HGDP 参考面板合并

流程：
  1. 找出用户数据与 HGDP 面板的共同 SNPs（通过 BIM 文件交集）
  2. 分别提取交集 SNPs
  3. plink --bmerge 合并
  4. 处理 strand flip（链翻转 SNPs）

注意：合并前会将用户样本 ID 标记为 "USER"，方便后续绘图识别。
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
    if result.returncode not in (0, 3):  # PLINK 有时以 3 表示警告但成功
        console.print(f"[red]PLINK {step} 失败（退出码 {result.returncode}）：[/red]")
        console.print(result.stderr or result.stdout)
        raise RuntimeError(f"PLINK 步骤失败：{step}")


def _load_bim(bim_path: str) -> pd.DataFrame:
    """加载 .bim 文件，返回 DataFrame。"""
    return pd.read_csv(
        bim_path, sep="\t", header=None,
        names=["chrom", "snp", "cm", "pos", "a1", "a2"],
    )


def _common_snps(user_bim: str, ref_bim: str) -> set[str]:
    """返回两个 BIM 文件中共同的 SNP ID 集合。"""
    user_df = _load_bim(user_bim)
    ref_df  = _load_bim(ref_bim)
    common = set(user_df["snp"]) & set(ref_df["snp"])
    return common


def merge_with_hgdp(
    user_bed: str,
    ref_dir: str,
    out_dir: str,
) -> str:
    """
    将用户 LD 剪枝数据与 HGDP 参考面板合并。

    Args:
        user_bed: 用户 BED 前缀（由 qc.py 返回）
        ref_dir:  HGDP 参考面板目录（含 hgdp_pruned.bed/.bim/.fam）
        out_dir:  工作目录

    Returns:
        合并后的 BED 前缀路径字符串
    """
    os.makedirs(out_dir, exist_ok=True)

    ref_prefix = str(Path(ref_dir) / "hgdp_pruned")
    merged_prefix = str(Path(out_dir) / "merged")

    # ── 检查参考面板文件 ─────────────────────────────────────────────────────
    for ext in (".bed", ".bim", ".fam"):
        if not Path(ref_prefix + ext).exists():
            raise FileNotFoundError(
                f"找不到 HGDP 参考面板文件：{ref_prefix + ext}\n"
                "请先运行：dna download"
            )

    # ── 步骤 1：找共同 SNPs ──────────────────────────────────────────────────
    console.print("  [bold]2a[/bold] 计算共同 SNPs")
    common = _common_snps(user_bed + ".bim", ref_prefix + ".bim")
    console.print(f"    用户 × HGDP 共同 SNPs：{len(common):,} 个")

    if len(common) < 1000:
        console.print(
            f"  [yellow]⚠ 共同 SNPs 较少（{len(common)} < 1000），"
            "结果可能不可靠。建议检查基因组版本是否一致（hg19 vs hg38）[/yellow]"
        )

    # 写出共同 SNP 列表
    common_list_path = str(Path(out_dir) / "common_snps.txt")
    with open(common_list_path, "w") as f:
        f.write("\n".join(sorted(common)))

    # ── 步骤 2：分别提取共同 SNPs ────────────────────────────────────────────
    console.print("  [bold]2b[/bold] 提取共同 SNPs（用户数据 & HGDP）")

    user_common = str(Path(out_dir) / "user_common")
    _run_plink([
        "--bfile", user_bed,
        "--extract", common_list_path,
        "--make-bed",
        "--out", user_common,
    ], step="用户数据提取共同 SNPs")

    ref_common = str(Path(out_dir) / "ref_common")
    _run_plink([
        "--bfile", ref_prefix,
        "--extract", common_list_path,
        "--make-bed",
        "--out", ref_common,
    ], step="HGDP 提取共同 SNPs")

    # ── 步骤 3：标记用户样本 ─────────────────────────────────────────────────
    console.print("  [bold]2c[/bold] 标记用户样本（FID → USER）")
    fam_path = user_common + ".fam"
    fam_df = pd.read_csv(
        fam_path, sep=r"\s+", header=None,
        names=["fid", "iid", "pat", "mat", "sex", "phen"],
    )
    fam_df["fid"] = "USER"
    fam_df.to_csv(fam_path, sep="\t", header=False, index=False)

    # ── 步骤 4：合并 ─────────────────────────────────────────────────────────
    console.print("  [bold]2d[/bold] bmerge 合并用户数据与 HGDP")
    _run_plink([
        "--bfile", ref_common,
        "--bmerge", user_common + ".bed", user_common + ".bim", user_common + ".fam",
        "--make-bed",
        "--out", merged_prefix,
        "--allow-no-sex",
    ], step="bmerge")

    # ── 处理 strand flip ────────────────────────────────────────────────────
    missnp_path = merged_prefix + "-merge.missnp"
    if Path(missnp_path).exists():
        console.print(f"  [yellow]检测到链翻转 SNPs，正在修正...[/yellow]")

        # 对用户数据 flip
        user_flip = str(Path(out_dir) / "user_flipped")
        _run_plink([
            "--bfile", user_common,
            "--flip", missnp_path,
            "--make-bed",
            "--out", user_flip,
        ], step="strand flip")

        # 重新合并
        _run_plink([
            "--bfile", ref_common,
            "--bmerge", user_flip + ".bed", user_flip + ".bim", user_flip + ".fam",
            "--make-bed",
            "--out", merged_prefix,
            "--allow-no-sex",
        ], step="re-merge after flip")

    from dna.qc import _count_samples, _count_variants
    n_sam = _count_samples(merged_prefix)
    n_snp = _count_variants(merged_prefix)
    console.print(
        f"  [green]✓[/green] 合并完成：{n_sam} 个样本"
        f"（含 {n_sam - 1} 个 HGDP 参考 + 用户样本），"
        f"{n_snp:,} 个 SNPs"
    )
    return merged_prefix
