"""
download.py — Download the HGDP + 1KG reference panel

Source:
  Zenodo record 10.5281/zenodo.14286454
  "BED/BIM/FAM files for HGDP + 1KG data from gnomAD v3.1.2"
  - Unrelated samples only
  - Variants filtered: AF > 1%, HWE p < 1e-12
  - ~3,000 samples, ~1M SNPs (will be further LD-pruned by qc.py)

Files:
  HGDP+1KG_SNPData.tar.gz             — BED/BIM/FAM archive
  hgdp_1kg_sample_info...tsv          — sample → population / superpopulation labels
"""

from __future__ import annotations

import hashlib
import os
import tarfile
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
# Source URLs  (Zenodo record 14286454, verified 2025-04)
# ──────────────────────────────────────────────────────────────────────────────

_ZENODO_BASE = "https://zenodo.org/records/14286454/files"

_SNP_ARCHIVE_URL = "https://zenodo.org/records/14286454/files/HGDP+1KG_SNPData.tar.gz?download=1"
_LABELS_URL = (
    "https://zenodo.org/records/14286454/files/"
    "hgdp_1kg_sample_info.unrelateds.pca_outliers_removed.with_project.tsv"
    "?download=1"
)

_SNP_ARCHIVE_NAME  = "HGDP+1KG_SNPData.tar.gz"
_LABELS_LOCAL_NAME = "hgdp_pop_labels.tsv"


# ──────────────────────────────────────────────────────────────────────────────
# Low-level helpers
# ──────────────────────────────────────────────────────────────────────────────

def _md5_file(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _download_file(url: str, dest: Path, label: str, force: bool = False) -> bool:
    """Stream-download a file with a rich progress bar. Returns True on success."""
    if dest.exists() and not force:
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
                task = progress.add_task(label, total=total or None)
                with open(dest, "wb") as fh:
                    for chunk in resp.iter_content(chunk_size=65536):
                        fh.write(chunk)
                        progress.advance(task, len(chunk))

        size_mb = dest.stat().st_size / 1e6
        console.print(f"  [green]✓[/green] {dest.name}  ({size_mb:.1f} MB)")
        return True

    except requests.RequestException as exc:
        console.print(f"  [red]✗ Download failed: {dest.name}[/red]  {exc}")
        if dest.exists():
            dest.unlink()
        return False


def _extract_tar(archive: Path, dest_dir: Path) -> list[Path]:
    """Extract a .tar.gz archive; return paths of extracted files."""
    console.print(f"  Extracting {archive.name} ...")
    extracted: list[Path] = []
    with tarfile.open(archive, "r:gz") as tf:
        for member in tf.getmembers():
            tf.extract(member, path=dest_dir)
            extracted.append(dest_dir / member.name)
    console.print(f"  [green]✓[/green] Extracted {len(extracted)} file(s)")
    return extracted


def _find_bed_prefix(dest_dir: Path) -> Path | None:
    """Locate the .bed file and return its prefix (without extension)."""
    for bed in dest_dir.rglob("*.bed"):
        return bed.with_suffix("")   # strip .bed
    return None


def _rename_to_canonical(bed_prefix: Path, dest_dir: Path) -> Path:
    """
    Rename BED/BIM/FAM to a fixed name 'hgdp_pruned.*' so downstream code
    can always reference the same prefix.
    """
    canon = dest_dir / "hgdp_pruned"
    if canon.with_suffix(".bed").exists():
        return canon  # already renamed
    for ext in (".bed", ".bim", ".fam"):
        src = bed_prefix.with_suffix(ext)
        dst = canon.with_suffix(ext)
        if src.exists() and not dst.exists():
            src.rename(dst)
    return canon


# ──────────────────────────────────────────────────────────────────────────────
# Public entry point
# ──────────────────────────────────────────────────────────────────────────────

def download_hgdp(out_dir: str = "data/hgdp", force: bool = False) -> Path:
    """
    Download the HGDP + 1KG reference panel (Zenodo 10.5281/zenodo.14286454).

    Steps:
      1. Download HGDP+1KG_SNPData.tar.gz
      2. Extract BED/BIM/FAM and rename to hgdp_pruned.*
      3. Download sample-info / population labels TSV

    Args:
        out_dir: Local directory to store files
        force:   If True, re-download even if files exist

    Returns:
        Path prefix of the reference panel BED (without extension),
        i.e. out_dir/hgdp_pruned

    Raises:
        RuntimeError: If a critical file fails to download or extract
    """
    dest_dir = Path(out_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)

    console.print("\n[bold cyan]Downloading HGDP + 1KG reference panel[/bold cyan]")
    console.print(f"  Source: Zenodo 10.5281/zenodo.14286454")
    console.print(f"  Destination: [dim]{dest_dir.resolve()}[/dim]\n")

    canon_prefix = dest_dir / "hgdp_pruned"

    # ── 1. BED/BIM/FAM archive ────────────────────────────────────────────────
    archive_path = dest_dir / _SNP_ARCHIVE_NAME

    bed_ready = (
        canon_prefix.with_suffix(".bed").exists() and
        canon_prefix.with_suffix(".bim").exists() and
        canon_prefix.with_suffix(".fam").exists()
    )

    if bed_ready and not force:
        console.print("  [dim]BED/BIM/FAM already present, skipping download.[/dim]")
    else:
        ok = _download_file(
            url=_SNP_ARCHIVE_URL,
            dest=archive_path,
            label="HGDP+1KG SNP data",
            force=force,
        )
        if not ok:
            raise RuntimeError(
                "Failed to download SNP archive from Zenodo.\n"
                "Check your connection or visit: https://doi.org/10.5281/zenodo.14286454"
            )

        extracted = _extract_tar(archive_path, dest_dir)

        bed_prefix = _find_bed_prefix(dest_dir)
        if bed_prefix is None:
            raise RuntimeError(
                f"No .bed file found after extracting {archive_path.name}. "
                f"Extracted files: {[p.name for p in extracted]}"
            )

        canon_prefix = _rename_to_canonical(bed_prefix, dest_dir)

        # Remove archive to save disk space
        archive_path.unlink(missing_ok=True)
        console.print("  Removed archive file to save disk space.")

    # ── 2. Population labels ──────────────────────────────────────────────────
    labels_path = dest_dir / _LABELS_LOCAL_NAME
    ok = _download_file(
        url=_LABELS_URL,
        dest=labels_path,
        label="Population labels",
        force=force,
    )
    if not ok:
        console.print(
            "  [yellow]Warning: population labels failed to download. "
            "Population colouring in plots will be disabled.[/yellow]"
        )

    # ── Done ──────────────────────────────────────────────────────────────────
    console.print(f"\n[bold green]Reference panel ready[/bold green]")
    console.print(f"  Prefix : [cyan]{canon_prefix}[/cyan]")
    console.print( "  Files  : hgdp_pruned.bed / .bim / .fam / hgdp_pop_labels.tsv\n")

    return canon_prefix
