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

    ADMIXTURE 1.3.0 requires:
      - Chromosomes in .bim must be plain integers 1-22 (no 'chr' prefix, no 0)
      - .bim/.fam must be in the same directory as .bed
      - .Q/.P are written to the working directory

    We therefore:
      1. Re-export the BED through PLINK (--chr 1-22) to guarantee numeric chrs
      2. Run ADMIXTURE from that directory using the filename only
      3. Move outputs to out_dir

    Returns a dict with keys: k, q_file, cv_error, log_file
    """
    import shutil

    os.makedirs(out_dir, exist_ok=True)

    bed_path = Path(bed + ".bed").resolve()
    bed_dir  = bed_path.parent
    out_path = Path(out_dir).resolve()
    log_path = out_path / f"admixture_K{k}.log"

    # ── Step: re-export through PLINK to guarantee ADMIXTURE-compatible format ──
    # ADMIXTURE 1.3.0 silently rejects BED files whose .bim has non-numeric
    # chromosome codes (chr1, 0, MT, …) and exits 255 with "Usage".
    clean_prefix = str(bed_dir / "admix_input")

    # Print plink version for diagnostics
    plink_ver = subprocess.run(["plink", "--version"], capture_output=True, text=True)
    console.print(f"  [dim]plink version: {plink_ver.stdout.strip() or plink_ver.stderr.strip()}[/dim]")

    plink_reexport = [
        "plink",
        "--bfile", bed_path.stem,
        "--chr", "1-22",          # keep only autosomes with numeric codes
        "--make-bed",
        "--out", "admix_input",
        "--allow-no-sex",
    ]
    console.print(f"  [dim]Preparing ADMIXTURE input (PLINK re-export)...[/dim]")
    prep = subprocess.run(
        plink_reexport,
        capture_output=True,
        text=True,
        cwd=str(bed_dir),
    )
    if prep.returncode not in (0, 3):
        raise RuntimeError(
            f"PLINK re-export for ADMIXTURE failed (exit {prep.returncode}):\n"
            + (prep.stderr or prep.stdout)
        )

    # Verify admix_input.bed exists and has correct PLINK 1.9 magic bytes (6c1b01)
    admix_bed = bed_dir / "admix_input.bed"
    admix_bim = bed_dir / "admix_input.bim"
    admix_fam = bed_dir / "admix_input.fam"
    for f in (admix_bed, admix_bim, admix_fam):
        if not f.exists() or f.stat().st_size == 0:
            raise RuntimeError(
                f"PLINK re-export did not produce: {f}\n"
                f"PLINK stdout:\n{prep.stdout}\nPLINK stderr:\n{prep.stderr}"
            )
    with open(admix_bed, "rb") as bf:
        magic = bf.read(3)
    console.print(f"  [dim]BED magic bytes: {magic.hex()} (need 6c1b01 for ADMIXTURE)[/dim]")
    if magic != b"\x6c\x1b\x01":
        raise RuntimeError(
            f"BED file has unexpected magic bytes: {magic.hex()}\n"
            "ADMIXTURE 1.3.0 requires PLINK 1.9 SNP-major BED (magic: 6c1b01).\n"
            "This may mean a PLINK 2 BED was produced. Check your plink version."
        )

    # Print first 3 lines of .bim for diagnosis
    if admix_bim.exists():
        with open(admix_bim) as bf:
            preview_lines = [bf.readline().rstrip() for _ in range(3)]
        console.print(
            f"  [dim].bim preview: " + " | ".join(preview_lines) + "[/dim]"
        )
    n_snp = sum(1 for _ in open(admix_bim))
    n_sam = sum(1 for _ in open(admix_fam))
    console.print(f"  [dim]admix_input: {n_sam} samples, {n_snp:,} SNPs[/dim]")

    stem = "admix_input"

    # ── Run ADMIXTURE ────────────────────────────────────────────────────────
    # ADMIXTURE 1.3.0 argument order: file and K must come before flags.
    # We try increasingly minimal invocations as fallback.
    bed_abs = str(admix_bed.resolve())
    candidate_cmds = [
        # Preferred: file K --cv -j N  (positional args first, flags after)
        ["admixture", bed_abs, str(k), "--cv", "-j", str(threads)],
        # Fallback 1: no threading flag
        ["admixture", bed_abs, str(k), "--cv"],
        # Fallback 2: bare minimum — no flags at all (CV error won't be available)
        ["admixture", bed_abs, str(k)],
    ]

    proc = None
    cmd  = None
    for attempt_cmd in candidate_cmds:
        console.print(f"  [dim]$ {' '.join(attempt_cmd)}[/dim]")
        console.print(f"  [bold]Running K={k}[/bold]  (log: {log_path.name})")
        with open(log_path, "w") as log_fh:
            proc = subprocess.run(
                attempt_cmd,
                stdout=log_fh,
                stderr=subprocess.STDOUT,
                cwd=str(bed_dir),
            )
        if proc.returncode == 0:
            cmd = attempt_cmd
            break
        # Check log for "Usage" — if so, try next invocation style
        try:
            with open(log_path) as lf:
                log_txt = lf.read()
        except OSError:
            log_txt = ""
        if "Usage:" not in log_txt:
            # A real error (not a file/arg parsing issue) — stop retrying
            cmd = attempt_cmd
            break
        console.print(
            f"  [yellow]Attempt failed (exit {proc.returncode}), trying simpler invocation...[/yellow]"
        )

    if proc.returncode != 0:
        tail = ""
        try:
            with open(log_path) as lf:
                lines = lf.readlines()
            tail = "\n" + "".join(lines[-30:])
        except OSError:
            pass
        raise RuntimeError(
            f"ADMIXTURE K={k} failed (exit {proc.returncode})\n"
            f"Log ({log_path}):{tail}"
        )

    # Move .Q and .P from bed_dir to out_dir
    for ext in (f".{k}.Q", f".{k}.P"):
        src = bed_dir / f"{stem}{ext}"
        if src.exists():
            dst = out_path / f"merged{ext}"   # rename to 'merged' for plot.py
            shutil.move(str(src), str(dst))

    raw_q = out_path / f"merged.{k}.Q"
    raw_p = out_path / f"merged.{k}.P"

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
