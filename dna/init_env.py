"""
init_env.py — 环境检查模块

检测项：
  - Python 版本（需要 3.11+）
  - conda 是否激活
  - PLINK 是否可用
  - ADMIXTURE 是否可用
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from dataclasses import dataclass, field

from rich.console import Console
from rich.table import Table
from rich import box

console = Console()


# ──────────────────────────────────────────────────────────────────────────────
# 数据结构
# ──────────────────────────────────────────────────────────────────────────────

@dataclass
class ToolStatus:
    name: str
    ok: bool
    version: str = ""
    path: str = ""
    note: str = ""


# ──────────────────────────────────────────────────────────────────────────────
# 单项检测
# ──────────────────────────────────────────────────────────────────────────────

def _check_python() -> ToolStatus:
    vi = sys.version_info
    ver = f"{vi.major}.{vi.minor}.{vi.micro}"
    ok = (vi.major, vi.minor) >= (3, 11)
    return ToolStatus(
        name="Python",
        ok=ok,
        version=ver,
        path=sys.executable,
        note="" if ok else "需要 Python 3.11+，请更新",
    )


def _check_conda() -> ToolStatus:
    """检查当前是否在 conda 环境中运行。"""
    import os
    conda_env = os.environ.get("CONDA_DEFAULT_ENV", "")
    conda_prefix = os.environ.get("CONDA_PREFIX", "")
    ok = bool(conda_env and conda_prefix)
    return ToolStatus(
        name="conda 环境",
        ok=ok,
        version=conda_env or "(未激活)",
        path=conda_prefix,
        note="" if ok else "请先运行：conda activate dna-ancestry",
    )


def _check_tool(name: str, version_flag: str = "--version") -> ToolStatus:
    """检查外部二进制工具。"""
    path = shutil.which(name)
    if path is None:
        return ToolStatus(
            name=name.upper(),
            ok=False,
            note=f"{name} 未在 PATH 中找到，请确认 conda 环境已激活",
        )
    try:
        result = subprocess.run(
            [name, version_flag],
            capture_output=True, text=True, timeout=10,
        )
        raw = (result.stdout or result.stderr or "").strip().splitlines()
        ver = raw[0] if raw else "unknown"
    except (subprocess.TimeoutExpired, OSError):
        ver = "unknown"

    return ToolStatus(name=name.upper(), ok=True, version=ver, path=path)


# ──────────────────────────────────────────────────────────────────────────────
# 汇总检查
# ──────────────────────────────────────────────────────────────────────────────

def run_check(verbose: bool = True) -> bool:
    """运行所有环境检查，返回 True 表示全部通过。"""
    statuses: list[ToolStatus] = [
        _check_python(),
        _check_conda(),
        _check_tool("plink"),
        _check_tool("admixture"),
    ]

    if verbose:
        _print_table(statuses)

    all_ok = all(s.ok for s in statuses)

    if verbose:
        if all_ok:
            console.print("\n[bold green]✓ 所有工具就绪，可以开始分析！[/bold green]")
        else:
            console.print("\n[bold yellow]⚠  部分工具未就绪，请参考上表中的提示。[/bold yellow]")
            console.print(
                "\n安装指南："
                "\n  [cyan]conda activate dna-ancestry[/cyan]"
                "\n  [cyan]conda install -c bioconda -c conda-forge plink admixture[/cyan]"
            )

    return all_ok


def _print_table(statuses: list[ToolStatus]) -> None:
    table = Table(
        title="🧬 diy-dna-ancestry 环境检查",
        box=box.ROUNDED,
        show_header=True,
        header_style="bold cyan",
        border_style="dim",
    )
    table.add_column("工具",    style="bold", min_width=14)
    table.add_column("状态",    justify="center", min_width=6)
    table.add_column("版本",    min_width=20)
    table.add_column("路径",    min_width=30)
    table.add_column("备注",    style="dim")

    for s in statuses:
        icon = "[green]✓[/green]" if s.ok else "[red]✗[/red]"
        table.add_row(
            s.name,
            icon,
            s.version or "—",
            s.path or "—",
            s.note or "",
        )

    console.print()
    console.print(table)
