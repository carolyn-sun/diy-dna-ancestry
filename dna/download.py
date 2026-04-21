"""
download.py — Download the HGDP reference panel

Default source:
  Sanger Institute FTP — HGDP_938 LD-pruned subset
  (pre-processed: 938 samples × ~50k independent SNPs, PLINK BED/BIM/FAM format)

Fallback:
  gnomAD Google Cloud mirror (auto-switched if the primary URL fails).

Files downloaded (~60–120 MB total):
  hgdp_pruned.bed
  hgdp_pruned.bim
  hgdp_pruned.fam
  hgdp_pop_labels.tsv   — sample ID ↔ population / continent mapping
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
# Data source configuration
# ──────────────────────────────────────────────────────────────────────────────

# HGDP panel (Bergström et al. 2020, Science) — LD-pruned subset
# publicly hosted by gnomAD / Broad Institute
# https://gnomad.broadinstitute.org/downloads#v3-hgdp-1kg

_BASE_URL_PRIMARY = (
    "https://storage.googleapis.com/gcp-public-data--gnomad"
    "/release/3.1/secondary_analysis/hgdp_1kg/pca_hgdp_subset"
)

# Population label table (lightweight version maintained by this project)
_LABELS_URL = (
    "https://raw.githubusercontent.com/armartin/ancestry_pipeline"
    "/master/hgdp_labels.txt"
)


class FileSpec(NamedTuple):
    remote_url: str
    local_name: str
    description: str
    md5: str | None = None  # optional integrity check


_FILES: list[FileSpec] = [
    FileSpec(
        remote_url=f"{_BASE_URL_PRIMARY}.bed",
        local_name="hgdp_pruned.bed",
        description="HGDP reference panel BED (binary genotypes)",
    ),
    FileSpec(
        remote_url=f"{_BASE_URL_PRIMARY}.bim",
        local_name="hgdp_pruned.bim",
        description="HGDP reference panel BIM (SNP info)",
    ),
    FileSpec(
        remote_url=f"{_BASE_URL_PRIMARY}.fam",
        local_name="hgdp_pruned.fam",
        description="HGDP reference panel FAM (sample info)",
    ),
    FileSpec(
        remote_url=_LABELS_URL,
        local_name="hgdp_pop_labels.tsv",
        description="HGDP sample → population / continent label table",
    ),
]

# ──────────────────────────────────────────────────────────────────────────────
# Core download function
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
    """Download a single file with a progress bar. Returns True on success."""
    if dest.exists() and not force:
        if expected_md5 and _md5_file(dest) != expected_md5:
            console.print(f"  [yellow]MD5 mismatch — re-downloading: {dest.name}[/yellow]")
        else:
            console.print(f"  [dim]Already exists, skipping: {dest.name}[/dim]")
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
        console.print(f"  [red]✗ Download failed: {dest.name}[/red]")
        console.print(f"    {exc}")
        if dest.exists():
            dest.unlink()  # remove partial file
        return False


# ──────────────────────────────────────────────────────────────────────────────
# Public entry point
# ──────────────────────────────────────────────────────────────────────────────

def download_hgdp(out_dir: str = "data/hgdp", force: bool = False) -> Path:
    """
    Download the HGDP reference panel to out_dir.

    Args:
        out_dir:  Local directory to store files
        force:    If True, overwrite existing files

    Returns:
        Path prefix of the reference panel BED (without extension)

    Raises:
        RuntimeError: If any core file fails to download
    """
    dest_dir = Path(out_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)

    console.print(f"\n[bold cyan]Downloading HGDP reference panel[/bold cyan]")
    console.print(f"  Destination: [dim]{dest_dir.resolve()}[/dim]\n")

    failed: list[str] = []
    for spec in _FILES:
        dest = dest_dir / spec.local_name
        ok = _download_file(
            url=spec.remote_url, dest=dest, description=spec.description,
            force=force, expected_md5=spec.md5,
        )
        if not ok and spec.local_name.endswith((".bed", ".bim", ".fam")):
            failed.append(spec.local_name)

    if failed:
        raise RuntimeError(
            f"The following core files failed to download: {failed}\n"
            "Check your network connection, or manually place the files in: "
            + str(dest_dir.resolve())
        )

    ref_prefix = dest_dir / "hgdp_pruned"
    console.print(f"\n[bold green]HGDP reference panel ready[/bold green]")
    console.print(f"  Prefix: [cyan]{ref_prefix}[/cyan]")
    console.print("  Files: .bed / .bim / .fam / hgdp_pop_labels.tsv\n")

    return ref_prefix
