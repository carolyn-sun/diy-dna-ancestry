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

# Detect WSL (Windows Subsystem for Linux)
IS_WSL=false
if grep -qi microsoft /proc/version 2>/dev/null; then
    IS_WSL=true
    warn "WSL (Windows Subsystem for Linux) detected."
    warn "The ADMIXTURE 1.3 64-bit binary may crash (SIGSEGV) under WSL."
    warn "setup.sh will automatically try to install a 32-bit fallback binary after env creation."
    echo ""
fi

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
else
    warn "ADMIXTURE not found in environment"
fi

# ── WSL: test admixture & install 32-bit fallback if needed ──────────────────
ADMIXTURE32_PATH=""
if $IS_WSL && $ADMIXTURE_OK; then
    header "WSL: Testing ADMIXTURE binary compatibility"

    # Quick crash test: run admixture with no args (expect exit 1, not -11/SIGSEGV)
    conda run -n "$ENV_NAME" admixture --version >/dev/null 2>&1
    ADMIX_EXIT=$?

    # Negative exit code = killed by signal (SIGSEGV = -11)
    if [ "$ADMIX_EXIT" -lt 0 ] 2>/dev/null || [ "$ADMIX_EXIT" -eq 139 ]; then
        warn "ADMIXTURE 64-bit binary crashed under WSL (exit $ADMIX_EXIT)."
        warn "Attempting to install a 32-bit static binary as a fallback..."
        echo ""

        _install_admixture32() {
            local DEST="${HOME}/bin/admixture32"
            mkdir -p "${HOME}/bin"

            # Primary: bioconda linux-32 conda package (extract binary directly)
            local PKG_URL="https://conda.anaconda.org/bioconda/linux-32/admixture-1.3.0-0.tar.bz2"
            local TMP_DIR; TMP_DIR=$(mktemp -d)

            info "Downloading 32-bit ADMIXTURE package..."
            if curl -fsSL "$PKG_URL" -o "${TMP_DIR}/admixture32.tar.bz2" 2>/dev/null; then
                # conda packages are bz2-compressed tarballs; extract the bin/admixture file
                if tar -xjf "${TMP_DIR}/admixture32.tar.bz2" -C "$TMP_DIR" \
                        --wildcards '*/admixture' 2>/dev/null; then
                    local BIN; BIN=$(find "$TMP_DIR" -name admixture -type f | head -1)
                    if [ -n "$BIN" ]; then
                        cp "$BIN" "$DEST"
                        chmod +x "$DEST"
                        rm -rf "$TMP_DIR"
                        echo "$DEST"
                        return 0
                    fi
                fi
            fi

            # Fallback: official ADMIXTURE download page (64-bit, but try anyway)
            warn "  32-bit package download failed; trying official 64-bit binary..."
            local OFFICIAL_URL="https://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz"
            if curl -fsSL "$OFFICIAL_URL" -o "${TMP_DIR}/admixture.tar.gz" 2>/dev/null; then
                tar -xzf "${TMP_DIR}/admixture.tar.gz" -C "$TMP_DIR" 2>/dev/null || true
                local BIN; BIN=$(find "$TMP_DIR" -name admixture -type f | head -1)
                if [ -n "$BIN" ]; then
                    cp "$BIN" "$DEST"
                    chmod +x "$DEST"
                    rm -rf "$TMP_DIR"
                    echo "$DEST"
                    return 0
                fi
            fi

            rm -rf "$TMP_DIR"
            return 1
        }

        # Install 32-bit support libs (Ubuntu/Debian-based WSL)
        if command -v apt-get &>/dev/null; then
            info "Installing 32-bit runtime libraries (may need sudo)..."
            sudo apt-get install -y lib32gcc-s1 libc6-i386 2>/dev/null || \
                warn "  Could not install lib32 (skipping; binary may still work)"
        fi

        ADMIXTURE32_PATH=$(_install_admixture32)
        if [ -n "$ADMIXTURE32_PATH" ]; then
            success "32-bit admixture installed: $ADMIXTURE32_PATH"
        else
            warn "Could not automatically install 32-bit admixture."
            warn "Use --nmf-fallback as an alternative (see README)."
        fi
    else
        success "ADMIXTURE 64-bit binary runs correctly under this WSL environment."
    fi
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

if [ -n "${ADMIXTURE32_PATH:-}" ]; then
    echo -e "    ${YELLOW}# WSL: use 32-bit admixture binary (64-bit crashed):${NC}"
    echo -e "    ${CYAN}dna run --vcf your.vcf --k 3,5 --admixture-bin ${ADMIXTURE32_PATH}${NC}"
elif $IS_WSL && $ADMIXTURE_OK; then
    echo -e "    ${CYAN}dna run --vcf your.vcf --k 3,5${NC}"
    echo -e "    ${YELLOW}# If ADMIXTURE crashes (SIGSEGV), try:${NC}"
    echo -e "    ${CYAN}dna run --vcf your.vcf --k 3,5 --admixture-bin ~/bin/admixture32${NC}"
    echo -e "    ${YELLOW}# or use the NMF approximation:${NC}"
    echo -e "    ${CYAN}dna run --vcf your.vcf --k 3,5 --nmf-fallback${NC}"
else
    echo -e "    ${CYAN}dna run --vcf your.vcf --k 3,5${NC}"
fi
echo ""
