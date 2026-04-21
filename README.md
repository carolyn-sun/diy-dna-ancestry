# diy-dna-ancestry

~~家用~~ DNA 祖源分析工具

A DIY tool for personal DNA ancestry analysis

## Requirements

- [Anaconda](https://www.anaconda.com/download) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
- The `.vcf` file from a DNA test

> **Apple Silicon:** `setup.sh` auto-sets `CONDA_SUBDIR=osx-64` so PLINK and ADMIXTURE
> resolve from bioconda and run transparently via Rosetta 2.

## Setup

**macOS / Linux**

```bash
bash setup.sh
conda activate dna-ancestry
```

## Pipelines

```
your.vcf
    │
    ▼
[QC & LD pruning] ──────────────── HGDP reference panel (~50k SNPs)
  geno / maf / hwe                          │
    │                                       │
    └──────────── [bmerge] ─────────────────┘
                      │
               ┌──────┴──────┐
               ▼             ▼
           [ADMIXTURE]     [PCA]
            K=3, K=5     PC1~PC10
               │             │
               └──────┬───────┘
                      ▼
               [matplotlib]
           dark-theme PCA + bar charts
```

## Full commands

```
dna init                              # check environment (conda, PLINK, ADMIXTURE)
dna download                          # download HGDP reference panel
dna run --vcf FILE [options]          # run the full pipeline
dna plot --results DIR                # re-plot from existing results
```

**`dna run` options**

| Flag          | Default      | Description                         |
| ------------- | ------------ | ----------------------------------- |
| `--vcf`       | _(required)_ | Input VCF file                      |
| `--k`         | `3,5`        | ADMIXTURE K values, comma-separated |
| `--threads`   | `4`          | Parallel threads                    |
| `--out`       | `results/`   | Output directory                    |
| `--geno`      | `0.05`       | Genotype missingness threshold      |
| `--maf`       | `0.01`       | Minimum allele frequency            |
| `--hwe`       | `1e-6`       | Hardy–Weinberg p-value cutoff       |
| `--skip-plot` | —            | Skip the plotting step              |

## Output structures

```
results/
├── pca_PC1_PC2.png       # PCA scatter plot (PC1 vs PC2)
├── pca_PC3_PC4.png       # PCA scatter plot (PC3 vs PC4)
├── admixture_K3.png      # ADMIXTURE bar chart, K=3
├── admixture_K5.png      # ADMIXTURE bar chart, K=5
├── cv_error.png          # CV error curve (best-K selection)
├── pca/
│   ├── pca.eigenvec
│   └── pca.eigenval
└── admixture/
    ├── merged.3.Q
    ├── merged.5.Q
    └── admixture_K*.log
```

## Development

```bash
conda activate dna-ancestry
pytest tests/ -v
```

## Data sources and tools

- [HGDP](https://www.internationalgenome.org/data-portal/data-collection/HGDP)
- [gnomAD](https://gnomad.broadinstitute.org/)
- [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/)
- [ADMIXTURE 1.3](https://dalexander.github.io/admixture/)
