"""
admixture.py — ADMIXTURE 无监督聚类模块

流程：
  对每个 K 值运行 admixture --cv merged.bed K
  解析 .Q 文件（个体祖源比例）和 CV error
  输出结果字典供 plot.py 使用
"""

from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path

from rich.console import Console
from rich.table import Table
from rich import box

console = Console()


# ──────────────────────────────────────────────────────────────────────────────
# 核心函数
# ──────────────────────────────────────────────────────────────────────────────

def _run_admixture_k(
    bed: str,
    k: int,
    out_dir: str,
    threads: int,
) -> dict:
    """
    跑单个 K 值的 ADMIXTURE。

    Returns dict with keys: k, q_file, cv_error, log_file
    """
    os.makedirs(out_dir, exist_ok=True)

    bed_path = Path(bed + ".bed")
    log_path = Path(out_dir) / f"admixture_K{k}.log"
    q_dest   = Path(out_dir) / f"merged.{k}.Q"
    p_dest   = Path(out_dir) / f"merged.{k}.P"

    # ADMIXTURE 会在当前目录输出 .Q 和 .P，需要 chdir
    cmd = [
        "admixture",
        "--cv",                  # 交叉验证（cross-validation）
        "-j", str(threads),      # 线程数
        str(bed_path.resolve()), # BED 文件完整路径
        str(k),
    ]

    console.print(f"  [dim]$ {' '.join(cmd)}[/dim]")
    console.print(f"  [bold]运行 K={k}[/bold]（日志：{log_path.name}）")

    with open(log_path, "w") as log_fh:
        proc = subprocess.run(
            cmd,
            stdout=log_fh,
            stderr=subprocess.STDOUT,
            cwd=out_dir,          # 输出文件将落在 out_dir
        )

    if proc.returncode != 0:
        raise RuntimeError(
            f"ADMIXTURE K={k} 失败（退出码 {proc.returncode}）\n"
            f"请查看日志：{log_path}"
        )

    # ADMIXTURE 输出的文件名基于输入文件名（去掉路径）
    stem = bed_path.stem  # e.g. "merged"
    raw_q = Path(out_dir) / f"{stem}.{k}.Q"
    raw_p = Path(out_dir) / f"{stem}.{k}.P"

    if not raw_q.exists():
        raise FileNotFoundError(
            f"ADMIXTURE 未生成 .Q 文件：{raw_q}\n"
            f"请查看日志：{log_path}"
        )

    # 解析 CV error
    cv_error = _parse_cv_error(log_path)
    if cv_error is not None:
        console.print(f"  [green]✓[/green] K={k} 完成，CV error = {cv_error:.6f}")
    else:
        console.print(f"  [green]✓[/green] K={k} 完成（CV error 未解析）")

    return {
        "k": k,
        "q_file": str(raw_q),
        "p_file": str(raw_p) if raw_p.exists() else None,
        "cv_error": cv_error,
        "log_file": str(log_path),
    }


def _parse_cv_error(log_path: Path) -> float | None:
    """从 ADMIXTURE 日志解析 CV error 值。"""
    pattern = re.compile(r"CV error \(K=\d+\):\s*([\d.]+)")
    try:
        with open(log_path) as f:
            for line in f:
                m = pattern.search(line)
                if m:
                    return float(m.group(1))
    except OSError:
        pass
    return None


def _print_cv_table(results: list[dict]) -> None:
    """打印 CV error 汇总表。"""
    table = Table(
        title="ADMIXTURE Cross-Validation Error",
        box=box.SIMPLE_HEAVY,
        show_header=True,
        header_style="bold cyan",
    )
    table.add_column("K", justify="center", style="bold")
    table.add_column("CV Error", justify="right")
    table.add_column("推荐", justify="center")

    valid = [r for r in results if r["cv_error"] is not None]
    if not valid:
        return

    best_k = min(valid, key=lambda r: r["cv_error"])["k"]

    for r in results:
        cv_str = f"{r['cv_error']:.6f}" if r["cv_error"] is not None else "—"
        rec    = "★ 最低" if r["k"] == best_k else ""
        table.add_row(str(r["k"]), cv_str, rec)

    console.print()
    console.print(table)


# ──────────────────────────────────────────────────────────────────────────────
# 公开入口
# ──────────────────────────────────────────────────────────────────────────────

def run_admixture(
    bed: str,
    ks: list[int],
    out_dir: str,
    threads: int = 4,
) -> dict[int, dict]:
    """
    对多个 K 值运行 ADMIXTURE。

    Args:
        bed:     合并后的 BED 前缀（不含 .bed 后缀）
        ks:      K 值列表，如 [3, 5]
        out_dir: 输出目录
        threads: 线程数

    Returns:
        dict，键为 K 值，值为包含 q_file / cv_error 等的字典
    """
    os.makedirs(out_dir, exist_ok=True)

    all_results: list[dict] = []
    for k in sorted(ks):
        console.print(f"\n  [cyan]── ADMIXTURE K={k} ────────────────────────[/cyan]")
        result = _run_admixture_k(bed=bed, k=k, out_dir=out_dir, threads=threads)
        all_results.append(result)

    _print_cv_table(all_results)

    return {r["k"]: r for r in all_results}
