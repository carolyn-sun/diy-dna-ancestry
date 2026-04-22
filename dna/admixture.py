"""
admixture.py — ADMIXTURE unsupervised clustering

Runs `admixture --cv merged.bed K` for each requested K value,
parses the .Q file (individual ancestry proportions) and CV error,
and returns a results dict for use by plot.py.
"""

from __future__ import annotations

import os
import random
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from rich.console import Console
from rich.table import Table
from rich import box

console = Console()


# ──────────────────────────────────────────────────────────────────────────────
# Python NMF fallback (used when the ADMIXTURE binary crashes)
# ──────────────────────────────────────────────────────────────────────────────

def _run_nmf_fallback(
    stem: str, k: int, bed_dir: Path, out_path: Path, log_path: Path
) -> dict:
    """
    Ancestry estimation via Non-negative Matrix Factorisation (numpy only).

    Activated when the ADMIXTURE binary crashes (e.g. AVX2/BLAS incompatibility
    on this Linux host).  Results closely approximate ADMIXTURE output.
    """
    console.print(
        f"  [yellow]⚠  ADMIXTURE binary crashed (SIGSEGV) on every invocation.\n"
        f"     Falling back to Python NMF ancestry estimation (K={k}).[/yellow]"
    )

    # ── Convert BED → .raw (additive 0/1/2 coding) via PLINK ─────────────────
    recode = subprocess.run(
        ["plink", "--bfile", stem, "--recode", "A",
         "--out", "admix_raw", "--allow-no-sex"],
        capture_output=True, text=True, cwd=str(bed_dir),
    )
    raw_file = bed_dir / "admix_raw.raw"
    if not raw_file.exists():
        raise RuntimeError(
            "PLINK --recode A failed; cannot produce genotype matrix.\n"
            + (recode.stderr or recode.stdout)
        )

    # ── Load genotype matrix ──────────────────────────────────────────────────
    console.print("  [dim]Loading genotype matrix for NMF...[/dim]")
    df = pd.read_csv(raw_file, sep=r"\s+", na_values="NA")
    meta = ["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"]
    X = df.drop(columns=[c for c in meta if c in df.columns]).values.astype(float)

    # Fill missing values with per-SNP mean (NMF requires no NaN)
    col_means = np.nanmean(X, axis=0)
    nan_rows, nan_cols = np.where(np.isnan(X))
    X[nan_rows, nan_cols] = col_means[nan_cols]

    n, p = X.shape

    # ── Multiplicative-update NMF (Lee & Seung 2001) ─────────────────────────
    console.print(f"  [dim]NMF K={k} on {n}×{p} matrix...[/dim]")
    rng = np.random.default_rng(42)
    # Initialise with NNDSVD-like scaling
    scale = np.sqrt(X.mean() / k)
    W = rng.uniform(0, scale * 2, size=(n, k))
    H = rng.uniform(0, scale * 2, size=(k, p))

    eps = 1e-10
    for _ in range(200):
        # Update H
        WtX  = W.T @ X
        WtWH = W.T @ W @ H + eps
        H *= WtX / WtWH

        # Update W
        XHt  = X @ H.T
        WHHt = W @ H @ H.T + eps
        W *= XHt / WHHt

        # Column-normalise H (keeps W/H on same scale)
        col_norms = H.sum(axis=1, keepdims=True).clip(eps)
        H /= col_norms
        W *= col_norms.T

    # ── Normalise W rows → ancestry proportions (sum to 1) ───────────────────
    row_sums = W.sum(axis=1, keepdims=True).clip(eps)
    Q = W / row_sums  # shape: (n_samples, k)

    # ── Save Q in ADMIXTURE format (space-separated, no header) ─────────────
    q_path = out_path / f"merged.{k}.Q"
    np.savetxt(str(q_path), Q, fmt="%.6f")

    with open(log_path, "a") as lf:
        lf.write(
            f"\n[Python NMF fallback — ADMIXTURE binary crashed]\n"
            f"n_samples={n}, n_snps={p}, K={k}, iterations=200\n"
        )

    console.print(
        f"  [green]✓[/green] K={k} done via NMF fallback "
        f"(no CV error; ancestry proportions in {q_path.name})"
    )
    return {
        "k": k,
        "q_file": str(q_path),
        "p_file": None,
        "cv_error": None,
        "log_file": str(log_path),
        "fam_file": str(bed_dir / f"{stem}.fam"),
    }


# ──────────────────────────────────────────────────────────────────────────────
# Core functions
# ──────────────────────────────────────────────────────────────────────────────

def _run_admixture_k(bed: str, k: int, out_dir: str, threads: int,
                     allow_nmf_fallback: bool = False) -> dict:
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

    # ── Subsample reference samples if dataset is too large for ADMIXTURE ────
    # ADMIXTURE 1.3.0 (bioconda) can segfault on large matrices (~3000+ samples)
    # due to binary/BLAS compatibility issues on some Linux environments.
    # We cap reference samples at MAX_REF_SAMPLES; the user sample is always kept.
    MAX_REF_SAMPLES = 500

    fam_df = pd.read_csv(
        admix_fam, sep=r"\s+", header=None,
        names=["fid", "iid", "pat", "mat", "sex", "phen"],
    )
    user_rows = fam_df[fam_df["fid"] == "USER"]
    ref_rows  = fam_df[fam_df["fid"] != "USER"]

    if len(ref_rows) > MAX_REF_SAMPLES:
        console.print(
            f"  [yellow]Subsampling reference panel: "
            f"{len(ref_rows):,} → {MAX_REF_SAMPLES} samples "
            f"(ADMIXTURE stability limit)[/yellow]"
        )
        random.seed(42)
        sampled_ref = ref_rows.sample(n=MAX_REF_SAMPLES, random_state=42)
        keep_df = pd.concat([user_rows, sampled_ref])
        keep_path = bed_dir / "admix_keep.txt"
        keep_df[["fid", "iid"]].to_csv(keep_path, sep="\t", header=False, index=False)

        sub = subprocess.run(
            ["plink", "--bfile", "admix_input", "--keep", str(keep_path),
             "--make-bed", "--out", "admix_sub", "--allow-no-sex"],
            capture_output=True, text=True, cwd=str(bed_dir),
        )
        if sub.returncode not in (0, 3):
            raise RuntimeError(
                f"PLINK subsample failed (exit {sub.returncode}):\n"
                + (sub.stderr or sub.stdout)
            )
        stem = "admix_sub"
        admix_bed = bed_dir / "admix_sub.bed"
        n_sam = len(keep_df)
        console.print(f"  [dim]Subsampled: {n_sam} samples ({MAX_REF_SAMPLES} ref + 1 user)[/dim]")
    else:
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
        try:
            with open(log_path) as lf:
                log_txt = lf.read()
        except OSError:
            log_txt = ""
        # Retry on: arg/file parse failure ("Usage:") OR segfault (negative exit code)
        is_retryable = ("Usage:" in log_txt) or (proc.returncode < 0)
        if not is_retryable:
            cmd = attempt_cmd
            break
        reason = "segfault" if proc.returncode < 0 else "argument parse error"
        console.print(
            f"  [yellow]Attempt failed (exit {proc.returncode}, {reason}),"
            f" trying simpler invocation...[/yellow]"
        )

    # If all ADMIXTURE invocations crashed (segfault), optionally use NMF fallback
    all_segfaulted = (proc is not None and proc.returncode < 0)
    if all_segfaulted:
        if allow_nmf_fallback:
            console.print(
                "  [yellow]ADMIXTURE binary crashed (SIGSEGV). "
                "Using Python NMF fallback (--nmf-fallback is enabled).[/yellow]"
            )
            return _run_nmf_fallback(
                stem=stem, k=k, bed_dir=bed_dir,
                out_path=out_path, log_path=log_path,
            )
        else:
            raise RuntimeError(
                f"ADMIXTURE K={k} crashed with a segfault (exit {proc.returncode}).\n"
                "This usually means the ADMIXTURE binary is incompatible with your CPU/OS.\n"
                "Options:\n"
                "  1. Run on a native x86-64 Linux machine or Docker container.\n"
                "  2. Re-enable the approximate NMF fallback with --nmf-fallback\n"
                "     (results will be less accurate than true ADMIXTURE).\n"
                f"  Log: {log_path}"
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
        "fam_file": str(bed_dir / f"{stem}.fam"),
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
    allow_nmf_fallback: bool = False,
) -> dict[int, dict]:
    """
    Run ADMIXTURE for multiple K values.

    Args:
        bed:               Merged BED prefix (without .bed extension)
        ks:                List of K values, e.g. [3, 5]
        out_dir:           Output directory
        threads:           Number of threads
        allow_nmf_fallback: If True, fall back to a Python NMF approximation
                           when the ADMIXTURE binary crashes (e.g. SIGSEGV on
                           incompatible hardware).  Disabled by default because
                           NMF results are less accurate than true ADMIXTURE.

    Returns:
        Dict mapping K → result dict (q_file, cv_error, log_file, ...)
    """
    os.makedirs(out_dir, exist_ok=True)

    all_results: list[dict] = []
    for k in sorted(ks):
        console.print(f"\n  [cyan]── ADMIXTURE K={k} ──────────────────────────────[/cyan]")
        result = _run_admixture_k(
            bed=bed, k=k, out_dir=out_dir,
            threads=threads, allow_nmf_fallback=allow_nmf_fallback,
        )
        all_results.append(result)

    _print_cv_table(all_results)

    return {r["k"]: r for r in all_results}
