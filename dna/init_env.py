"""
init_env.py — Environment check module

Checks:
  - Python version (>= 3.11)
  - Whether a conda environment is active
  - Whether PLINK is available
  - Whether ADMIXTURE is available
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from dataclasses import dataclass

from rich.console import Console
from rich.table import Table
from rich import box

console = Console()


# ──────────────────────────────────────────────────────────────────────────────
# Data structures
# ──────────────────────────────────────────────────────────────────────────────

@dataclass
class ToolStatus:
    name: str
    ok: bool
    version: str = ""
    path: str = ""
    note: str = ""


# ──────────────────────────────────────────────────────────────────────────────
# Individual checks
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
        note="" if ok else "Python 3.11+ required — please upgrade",
    )


def _check_conda() -> ToolStatus:
    """Check whether we are running inside a conda environment."""
    import os
    conda_env = os.environ.get("CONDA_DEFAULT_ENV", "")
    conda_prefix = os.environ.get("CONDA_PREFIX", "")
    ok = bool(conda_env and conda_prefix)
    return ToolStatus(
        name="conda env",
        ok=ok,
        version=conda_env or "(not active)",
        path=conda_prefix,
        note="" if ok else "Run: conda activate dna-ancestry",
    )


def _check_tool(name: str, version_flag: str = "--version") -> ToolStatus:
    """Check whether an external binary is available."""
    path = shutil.which(name)
    if path is None:
        return ToolStatus(
            name=name.upper(),
            ok=False,
            note=f"{name} not found in PATH — make sure the conda env is active",
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
# Aggregate check
# ──────────────────────────────────────────────────────────────────────────────

def run_check(verbose: bool = True) -> bool:
    """Run all environment checks. Returns True if all pass."""
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
            console.print("\n[bold green]All tools ready — you can start the analysis![/bold green]")
        else:
            console.print("\n[bold yellow]Some tools are missing. See the notes above.[/bold yellow]")
            console.print(
                "\nInstallation hint:\n"
                "  [cyan]conda activate dna-ancestry[/cyan]\n"
                "  [cyan]conda install -c bioconda -c conda-forge plink admixture[/cyan]"
            )

    return all_ok


def _print_table(statuses: list[ToolStatus]) -> None:
    table = Table(
        title="🧬 diy-dna-ancestry environment check",
        box=box.ROUNDED,
        show_header=True,
        header_style="bold cyan",
        border_style="dim",
    )
    table.add_column("Tool",    style="bold", min_width=14)
    table.add_column("Status",  justify="center", min_width=6)
    table.add_column("Version", min_width=20)
    table.add_column("Path",    min_width=30)
    table.add_column("Note",    style="dim")

    for s in statuses:
        icon = "[green]OK[/green]" if s.ok else "[red]MISSING[/red]"
        table.add_row(
            s.name,
            icon,
            s.version or "—",
            s.path or "—",
            s.note or "",
        )

    console.print()
    console.print(table)
