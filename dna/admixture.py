"""
admixture.py — ADMIXTURE unsupervised clustering

Runs `admixture --cv merged.bed K` for each requested K value,
parses the .Q file (individual ancestry proportions) and CV error,
and returns a results dict for use by plot.py.
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
# Core functions
# ──────────────────────────────────────────────────────────────────────────────

def _run_admixture_k(bed: str, k: int, out_dir: str, threads: int) -> dict:
    """
    Run ADMIXTURE for a single K value.

    Returns a dict with keys: k, q_file, cv_error, log_file
    """
    os.makedirs(out_dir, exist_ok=True)

    bed_path = Path(bed + ".bed")
    log_path = Path(out_dir) / f"admixture_K{k}.log"

    cmd = [
        "admixture",
        "--cv",
        "-j", str(threads),
        str(bed_path.resolve()),
        str(k),
    ]

    console.print(f"  [dim]$ {' '.join(cmd)}[/dim]")
    console.print(f"  [bold]Running K={k}[/bold]  (log: {log_path.name})")

    with open(log_path, "w") as log_fh:
        proc = subprocess.run(
            cmd,
            stdout=log_fh,
            stderr=subprocess.STDOUT,
            cwd=out_dir,  # ADMIXTURE writes .Q/.P into cwd
        )

    if proc.returncode != 0:
        raise RuntimeError(
            f"ADMIXTURE K={k} failed (exit {proc.returncode})\n"
            f"See log: {log_path}"
        )

    stem  = bed_path.stem  # e.g. "merged"
    raw_q = Path(out_dir) / f"{stem}.{k}.Q"
    raw_p = Path(out_dir) / f"{stem}.{k}.P"

    if not raw_q.exists():
        raise FileNotFoundError(
            f"ADMIXTURE did not produce a .Q file: {raw_q}\n"
            f"See log: {log_path}"
        )

    cv_error = _parse_cv_error(log_path)
    if cv_error is not None:
        console.print(f"  [green]✓[/green] K={k} done, CV error = {cv_error:.6f}")
    else:
        console.print(f"  [green]✓[/green] K={k} done (CV error not parsed)")

    return {
        "k": k,
        "q_file": str(raw_q),
        "p_file": str(raw_p) if raw_p.exists() else None,
        "cv_error": cv_error,
        "log_file": str(log_path),
    }


def _parse_cv_error(log_path: Path) -> float | None:
    """Parse the CV error value from an ADMIXTURE log file."""
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
    """Print a CV error summary table."""
    table = Table(
        title="ADMIXTURE Cross-Validation Error",
        box=box.SIMPLE_HEAVY,
        show_header=True,
        header_style="bold cyan",
    )
    table.add_column("K",        justify="center", style="bold")
    table.add_column("CV Error", justify="right")
    table.add_column("Best",     justify="center")

    valid = [r for r in results if r["cv_error"] is not None]
    if not valid:
        return

    best_k = min(valid, key=lambda r: r["cv_error"])["k"]

    for r in results:
        cv_str = f"{r['cv_error']:.6f}" if r["cv_error"] is not None else "—"
        best   = "★ lowest" if r["k"] == best_k else ""
        table.add_row(str(r["k"]), cv_str, best)

    console.print()
    console.print(table)


# ──────────────────────────────────────────────────────────────────────────────
# Public entry point
# ──────────────────────────────────────────────────────────────────────────────

def run_admixture(
    bed: str,
    ks: list[int],
    out_dir: str,
    threads: int = 4,
) -> dict[int, dict]:
    """
    Run ADMIXTURE for multiple K values.

    Args:
        bed:     Merged BED prefix (without .bed extension)
        ks:      List of K values, e.g. [3, 5]
        out_dir: Output directory
        threads: Number of threads

    Returns:
        Dict mapping K → result dict (q_file, cv_error, log_file, ...)
    """
    os.makedirs(out_dir, exist_ok=True)

    all_results: list[dict] = []
    for k in sorted(ks):
        console.print(f"\n  [cyan]── ADMIXTURE K={k} ──────────────────────────────[/cyan]")
        result = _run_admixture_k(bed=bed, k=k, out_dir=out_dir, threads=threads)
        all_results.append(result)

    _print_cv_table(all_results)

    return {r["k"]: r for r in all_results}
