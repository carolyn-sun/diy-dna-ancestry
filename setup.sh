#!/usr/bin/env bash
# =============================================================================
# setup.sh — diy-dna-ancestry environment initialisation (conda)
# Usage: bash setup.sh
# =============================================================================

set -euo pipefail

# ── Colour helpers ────────────────────────────────────────────────────────────
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

info()    { echo -e "${CYAN}[INFO]${NC} $*"; }
success() { echo -e "${GREEN}[OK]${NC}   $*"; }
warn()    { echo -e "${YELLOW}[WARN]${NC} $*"; }
error()   { echo -e "${RED}[ERR]${NC}  $*"; }
header()  { echo -e "\n${BOLD}${CYAN}══════════════════════════════════════════${NC}"; \
            echo -e "${BOLD}  $*${NC}"; \
            echo -e "${BOLD}${CYAN}══════════════════════════════════════════${NC}"; }

# ── Step 1: Check conda ───────────────────────────────────────────────────────
header "Step 1 / 3: Check conda installation"

if ! command -v conda &>/dev/null; then
    error "conda not found. Please install Miniconda or Anaconda first:"
    echo ""
    echo "  macOS / Linux: https://docs.conda.io/en/latest/miniconda.html"
    echo ""
    echo "  Quick install (macOS, Apple Silicon):"
    echo "    curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh"
    echo "    bash Miniconda3-latest-MacOSX-arm64.sh -b -p \$HOME/miniconda3"
    echo "    \$HOME/miniconda3/bin/conda init zsh   # or: bash"
    echo ""
    echo "  Quick install (macOS, Intel):"
    echo "    curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh"
    echo "    bash Miniconda3-latest-MacOSX-x86_64.sh -b -p \$HOME/miniconda3"
    echo "    \$HOME/miniconda3/bin/conda init zsh"
    exit 1
fi

CONDA_VER=$(conda --version)
success "$CONDA_VER is ready"

# Detect platform
ARCH=$(uname -m)
OS=$(uname -s)
info "Platform: $OS / $ARCH"

# Apple Silicon: plink & admixture are x86_64-only on bioconda.
# Force osx-64 subdir so conda fetches the x86_64 packages (runs via Rosetta 2).
if [ "$OS" == "Darwin" ] && [ "$ARCH" == "arm64" ]; then
    warn "Apple Silicon (arm64) detected."
    warn "plink and admixture are x86_64 packages on bioconda — setting CONDA_SUBDIR=osx-64."
    warn "They will run transparently via Rosetta 2."
    export CONDA_SUBDIR=osx-64
    echo ""
fi

# ── Step 2: Create / update conda environment ─────────────────────────────────
header "Step 2 / 3: Create conda environment (dna-ancestry)"

ENV_NAME="dna-ancestry"

if conda env list | grep -q "^${ENV_NAME} "; then
    warn "Environment '${ENV_NAME}' already exists."
    read -p "  Update it now (conda env update)? [y/N] " -n 1 -r REPLY
    echo ""
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        info "Updating conda environment..."
        conda env update -f environment.yml --prune
        success "Environment updated"
    else
        info "Skipping update — using existing environment"
    fi
else
    info "Creating conda environment (first run may take a few minutes)..."
    conda env create -f environment.yml
    success "conda environment '${ENV_NAME}' created"
    # Pin the subdir inside the env so subsequent 'conda install' calls stay on osx-64
    if [ -n "${CONDA_SUBDIR:-}" ]; then
        conda config --env --set subdir osx-64
        info "Pinned subdir=osx-64 inside the environment."
    fi
fi

# ── Step 3: Verify tools ──────────────────────────────────────────────────────
header "Step 3 / 3: Verify tools"

PYTHON_OK=false
PLINK_OK=false
ADMIXTURE_OK=false

if conda run -n "$ENV_NAME" python --version &>/dev/null 2>&1; then
    PY_VER=$(conda run -n "$ENV_NAME" python --version 2>&1)
    success "Python ($PY_VER)"
    PYTHON_OK=true
fi

if conda run -n "$ENV_NAME" plink --version &>/dev/null 2>&1; then
    PL_VER=$(conda run -n "$ENV_NAME" plink --version 2>/dev/null | head -1 || echo "unknown")
    success "PLINK ($PL_VER)"
    PLINK_OK=true
else
    warn "PLINK not found in environment"
fi

if conda run -n "$ENV_NAME" admixture --version &>/dev/null 2>&1; then
    AD_VER=$(conda run -n "$ENV_NAME" admixture --version 2>/dev/null | head -1 || echo "unknown")
    success "ADMIXTURE ($AD_VER)"
    ADMIXTURE_OK=true

    # ── Auto-upgrade 1.3.0 → 1.3.1 on Linux x86_64 ──────────────────────────
    # ADMIXTURE 1.3.0 has a known SIGSEGV crash on modern Linux kernels
    # when allocating the genotype matrix.  1.3.1 fixes this.
    if [ "$OS" = "Linux" ] && [ "$ARCH" = "x86_64" ]; then
        AD_VER_NUM=$(conda run -n "$ENV_NAME" admixture --version 2>/dev/null \
                     | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
        if [ "$AD_VER_NUM" = "1.3.0" ]; then
            warn "ADMIXTURE 1.3.0 detected — upgrading to 1.3.1 (fixes SIGSEGV on modern kernels)..."
            CONDA_BIN=$(conda run -n "$ENV_NAME" python -c \
                        "import sys; print(sys.prefix)" 2>/dev/null)/bin
            TMP_DIR=$(mktemp -d)
            URL="https://dalexander.github.io/admixture/binaries/admixture_linux-1.3.1.tar.gz"
            if curl -fsSL "$URL" -o "$TMP_DIR/admixture131.tar.gz" 2>/dev/null; then
                tar -xzf "$TMP_DIR/admixture131.tar.gz" -C "$TMP_DIR" 2>/dev/null
                BIN=$(find "$TMP_DIR" -name admixture -type f | head -1)
                if [ -n "$BIN" ]; then
                    cp "$BIN" "$CONDA_BIN/admixture"
                    chmod +x "$CONDA_BIN/admixture"
                    success "ADMIXTURE upgraded to 1.3.1"
                else
                    warn "Could not extract admixture binary from archive; staying on 1.3.0"
                fi
            else
                warn "Failed to download ADMIXTURE 1.3.1 (no internet?); staying on 1.3.0"
            fi
            rm -rf "$TMP_DIR"
        fi
    fi
else
    warn "ADMIXTURE not found in environment"
fi

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo -e "${BOLD}${CYAN}══════════════════════════════════════════${NC}"
echo -e "${BOLD}  Summary${NC}"
echo -e "${BOLD}${CYAN}══════════════════════════════════════════${NC}"

icon() { $1 && echo -e "${GREEN}✓${NC}" || echo -e "${RED}✗${NC}"; }

echo -e "  Python     $(icon $PYTHON_OK)"
echo -e "  PLINK      $(icon $PLINK_OK)"
echo -e "  ADMIXTURE  $(icon $ADMIXTURE_OK)"
echo ""

if $PYTHON_OK && $PLINK_OK && $ADMIXTURE_OK; then
    echo -e "${GREEN}${BOLD}✓ All tools ready!${NC}"
else
    echo -e "${YELLOW}⚠  Some tools are missing. Review the warnings above.${NC}"
    echo "   After fixing, run 'dna init' for a detailed diagnosis."
fi

echo ""
echo -e "  Next steps:"
echo -e "    ${CYAN}conda activate ${ENV_NAME}${NC}"
echo -e "    ${CYAN}dna init${NC}                              # detailed environment check"
echo -e "    ${CYAN}dna download${NC}                          # download HGDP reference panel"
echo -e "    ${CYAN}dna run --vcf your.vcf --k 3,5${NC}"
echo ""
