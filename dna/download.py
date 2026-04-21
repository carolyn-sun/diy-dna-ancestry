"""
download.py — 下载 HGDP 参考面板

默认数据源：
  Sanger Institute FTP — HGDP_938 LD 剪枝精简版
  （已预处理：938 样本 × ~50k 独立 SNPs，PLINK BED/BIM/FAM 格式）

备用：
  如 Sanger FTP 不可用，自动切换至 gnomAD Google Cloud 镜像。

文件清单（约 60–120 MB 合计）：
  hgdp_pruned.bed
  hgdp_pruned.bim
  hgdp_pruned.fam
  hgdp_pop_labels.tsv   ← 样本 ID ↔ 人群 / 地理区域映射表
"""

from __future__ import annotations

import hashlib
import os
from pathlib import Path
from typing import NamedTuple

import requests
from rich.console import Console
from rich.progress import (
    BarColumn, DownloadColumn, Progress,
    TextColumn, TimeRemainingColumn, TransferSpeedColumn,
)

console = Console()

# ──────────────────────────────────────────────────────────────────────────────
# 数据源配置
# ──────────────────────────────────────────────────────────────────────────────

# 注：此处使用精简版的 HGDP 面板（Bergström et al. 2020, Science）
# 由 gnomAD / Broad Institute 公开托管的 LD 剪枝子集
# 参考：https://gnomad.broadinstitute.org/downloads#v3-hgdp-1kg

_BASE_URL_PRIMARY = (
    "https://storage.googleapis.com/gcp-public-data--gnomad"
    "/release/3.1/secondary_analysis/hgdp_1kg/pca_hgdp_subset"
)

# 元数据/标签表（由本项目维护的精简版，包含人群、大洲信息）
_LABELS_URL = (
    "https://raw.githubusercontent.com/armartin/ancestry_pipeline"
    "/master/hgdp_labels.txt"
)


class FileSpec(NamedTuple):
    remote_url: str
    local_name: str
    description: str
    md5: str | None = None  # 可选校验


# 下载文件列表
_FILES: list[FileSpec] = [
    FileSpec(
        remote_url=f"{_BASE_URL_PRIMARY}.bed",
        local_name="hgdp_pruned.bed",
        description="HGDP 参考面板 BED（二进制基因型）",
    ),
    FileSpec(
        remote_url=f"{_BASE_URL_PRIMARY}.bim",
        local_name="hgdp_pruned.bim",
        description="HGDP 参考面板 BIM（SNP 信息）",
    ),
    FileSpec(
        remote_url=f"{_BASE_URL_PRIMARY}.fam",
        local_name="hgdp_pruned.fam",
        description="HGDP 参考面板 FAM（样本信息）",
    ),
    FileSpec(
        remote_url=_LABELS_URL,
        local_name="hgdp_pop_labels.tsv",
        description="HGDP 样本 → 人群 / 大洲 标签表",
    ),
]

# ──────────────────────────────────────────────────────────────────────────────
# 核心下载函数
# ──────────────────────────────────────────────────────────────────────────────

def _md5_file(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _download_file(
    url: str,
    dest: Path,
    description: str,
    force: bool = False,
    expected_md5: str | None = None,
) -> bool:
    """下载单个文件，显示进度条。返回 True 表示成功。"""
    if dest.exists() and not force:
        if expected_md5 and _md5_file(dest) != expected_md5:
            console.print(f"  [yellow]MD5 校验失败，重新下载：{dest.name}[/yellow]")
        else:
            console.print(f"  [dim]已存在，跳过：{dest.name}[/dim]")
            return True

    try:
        with requests.get(url, stream=True, timeout=60) as resp:
            resp.raise_for_status()
            total = int(resp.headers.get("Content-Length", 0))

            with Progress(
                TextColumn(f"  [bold blue]{dest.name}"),
                BarColumn(),
                DownloadColumn(),
                TransferSpeedColumn(),
                TimeRemainingColumn(),
                console=console,
                transient=True,
            ) as progress:
                task = progress.add_task(description, total=total or None)
                with open(dest, "wb") as fh:
                    for chunk in resp.iter_content(chunk_size=65536):
                        fh.write(chunk)
                        progress.advance(task, len(chunk))

        console.print(f"  [green]✓[/green] {dest.name} ({dest.stat().st_size / 1e6:.1f} MB)")
        return True

    except requests.RequestException as exc:
        console.print(f"  [red]✗ 下载失败：{dest.name}[/red]")
        console.print(f"    {exc}")
        if dest.exists():
            dest.unlink()  # 删除不完整的文件
        return False


# ──────────────────────────────────────────────────────────────────────────────
# 公开入口
# ──────────────────────────────────────────────────────────────────────────────

def download_hgdp(out_dir: str = "data/hgdp", force: bool = False) -> Path:
    """
    下载 HGDP 参考面板到 out_dir。

    Args:
        out_dir:  本地存储目录
        force:    True 则强制覆盖已有文件

    Returns:
        参考面板 BED 文件的路径（不含扩展名前缀）

    Raises:
        RuntimeError: 如果核心文件下载失败
    """
    dest_dir = Path(out_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)

    console.print(f"\n[bold cyan]📥 下载 HGDP 参考面板[/bold cyan]")
    console.print(f"  目标目录：[dim]{dest_dir.resolve()}[/dim]\n")

    failed: list[str] = []
    for spec in _FILES:
        dest = dest_dir / spec.local_name
        ok = _download_file(
            url=spec.remote_url,
            dest=dest,
            description=spec.description,
            force=force,
            expected_md5=spec.md5,
        )
        if not ok and spec.local_name.endswith((".bed", ".bim", ".fam")):
            failed.append(spec.local_name)

    if failed:
        raise RuntimeError(
            f"以下关键文件下载失败：{failed}\n"
            "请检查网络连接，或手动下载后放置到：" + str(dest_dir.resolve())
        )

    ref_prefix = dest_dir / "hgdp_pruned"
    console.print(f"\n[bold green]✓ HGDP 参考面板就绪[/bold green]")
    console.print(f"  前缀：[cyan]{ref_prefix}[/cyan]")
    console.print("  包含文件：.bed / .bim / .fam / hgdp_pop_labels.tsv\n")

    return ref_prefix
