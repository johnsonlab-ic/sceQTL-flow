# Unified sceQTL-flow Docker Image

This document describes the unified Docker image approach for sceQTL-flow, consolidating all pipeline dependencies into a single container.

## Overview

**Status:** ✅ Ready for Unified Deployment

The sceQTL-flow pipeline originally uses three separate Docker images:
- `eqtl-genotype` - Genotype processing, eQTL testing, PC optimization
- `eqtl-expression` - Single-cell data handling, pseudobulking  
- `eqtl-reports` - Report generation and visualization

The new unified approach consolidates all dependencies into a single image stored in `.devcontainer/Dockerfile`. This simplifies:
- Build pipelines
- Registry management
- CI/CD workflows
- Local development

## Image Contents

The unified image includes all necessary components:

### System Dependencies
- R 4.4.0 base
- Development tools (gcc, g++, gfortran, make, bzip2)
- HDF5, BLAS/LAPACK libraries
- Nextflow, PLINK 1.9 & 2.0, bcftools
- Python 3 with pip
- Build essentials for R package compilation

### R Packages

**Data Processing & Single-cell Analysis:**
- `Seurat` (v5) - Single-cell analysis
- `Signac` - ATAC-seq analysis
- `BPCells` - Efficient sparse matrix operations
- `dplyr`, `data.table`, `tidyverse` - Data manipulation

**Genomics & QTL Analysis:**
- `SeqArray` - GDS file handling
- `SNPRelate` - SNP analysis
- `GenomicRanges`, `BSgenome` - Genomic coordinates
- `MatrixEQTL` - eQTL testing
- `SNPlocs.Hsapiens.dbSNP155.GRCh38` - SNP annotation
- `SNPlocs.Hsapiens.dbSNP155.GRCh37` - Liftover support

**Statistical & Reporting:**
- `ggplot2`, `cowplot`, `patchwork` - Visualization
- `ggrepel`, `ggbeeswarm`, `gridExtra` - Plot enhancements
- `rmarkdown`, `tinytex` - Report generation
- `edgeR`, `presto` - Expression analysis
- `coloc`, `MendelianRandomization`, `ieugwasr` - Advanced QTL tools

### Python Packages
- `tensorqtl` - GPU-accelerated eQTL testing
- `numpy`, `pandas`, `scipy` - Scientific computing
- `Pgenlib`, `pandas-plink` - PLINK format support
- `scikit-learn`, `rpy2` - ML and R integration

### External Data
- liftOver chain file (hg19ToHg38) for coordinate conversion
- `/usr/local/src/hg19ToHg38.over.chain.gz` unpacked and ready

## Building the Image

### Quick Build

```bash
cd Images
chmod +x build.sh
./build.sh
```

### Build with Push to Registry

```bash
./build.sh --push
```

Make sure you're authenticated:
```bash
echo $CR_PAT | docker login ghcr.io -u <username> --password-stdin
```

### Build with Custom Tag

```bash
./build.sh --tag v1.0.0 --push
./build.sh --tag latest --push
./build.sh --no-cache --tag v1.0.1 --push
```

### Build with Force Rebuild

```bash
./build.sh --no-cache --push
```

## Using the Unified Image

### Update nextflow.config

Replace container specifications for all processes with the unified image:

```nextflow
profiles {
    offline {
        docker.enabled = true
        process {
            executor = "local"
            withName: /.*/ {
                container = "ghcr.io/johnsonlab-ic/sceqtl-flow:latest"
            }
        }
    }

    imperial {
        docker.enabled = false
        singularity.enabled = true
        process {
            executor = 'pbspro'
            withName: /.*/ {
                container = "docker://ghcr.io/johnsonlab-ic/sceqtl-flow:latest"
            }
        }
    }
}
```

### Test the Pipeline

```bash
cd personal
./run_sceQTL-flow.sh
```

## Verification Checklist

✅ **All genotype processes covered:**
- Genotype creation and QC
- PC optimization (coarse and fine stage)
- PC selection
- MatrixEQTL eQTL testing
- Residual calculation

✅ **All expression processes covered:**
- Pseudobulking (RNA & ATAC)
- QC and merge operations
- Cell type counting

✅ **All report generation covered:**
- RMarkdown report generation
- Visualization packages all available
- LaTeX for PDF output

✅ **Optional advanced analyses:**
- TensorQTL (GPU-optimized eQTL) 
- Coloc and MR tools for cross-trait analysis

## Image Size

The unified image is larger than individual images but provides:
- **Faster startup** - Single pull instead of multiple downloads
- **Simpler management** - One image to maintain
- **Consistent environment** - All processes use identical dependencies

Typical size: ~3.5 GB (comparable to building all 3 images separately)

## Development & Debugging

### Local Development in devcontainer

VS Code users can develop directly in the container:

```bash
# Open in devcontainer
code /path/to/sceQTL-flow
# Select "Reopen in Container" from command palette
```

The `.devcontainer/devcontainer.json` configures the development environment automatically.

### Manual Container Testing

```bash
# Build locally
docker build -t sceqtl-flow:test -f .devcontainer/Dockerfile .

# Run interactively
docker run -it sceqtl-flow:test /bin/bash

# Verify key tools
R --vanilla -e "library(Seurat); library(MatrixEQTL)" 
python3 -c "import tensorqtl; print(tensorqtl.__version__)"
which plink plink2 bcftools nextflow
```

## Migration Guide

### If You Have Existing Pipelines

1. **Backup current nextflow.config**
   ```bash
   cp nextflow.config nextflow.config.backup
   ```

2. **Build unified image**
   ```bash
   cd Images && ./build.sh --push
   ```

3. **Update nextflow.config** - Replace all three container specifications with single unified image

4. **Test run on small dataset first**
   ```bash
   cd personal && ./run_sceQTL-flow.sh
   ```

5. **Compare results** with previous runs

### Rollback

If needed, restore original configuration:
```bash
cp nextflow.config.backup nextflow.config
docker pull ghcr.io/johnsonlab-ic/eqtl-genotype:latest  # etc.
```

## Troubleshooting

### Build Fails

```bash
# Check Docker daemon
docker ps

# Clean up dangling images
docker image prune -a

# Rebuild from scratch
./build.sh --no-cache
```

### Missing Package at Runtime

Verify all packages are properly installed:

```bash
docker run ghcr.io/johnsonlab-ic/sceqtl-flow:latest R -e "
  required_pkgs <- c('Seurat', 'MatrixEQTL', 'SeqArray', 'ggplot2', 'rmarkdown')
  for (pkg in required_pkgs) {
    if (!require(pkg, character.only = TRUE)) stop(paste0('Missing: ', pkg))
  }
  cat('All packages verified!\n')
"
```

### Python Tensorqtl Issues

```bash
docker run ghcr.io/johnsonlab-ic/sceqtl-flow:latest python3 -c "
import tensorqtl
import torch
print(f'TensorQTL version: {tensorqtl.__version__}')
print(f'CUDA available: {torch.cuda.is_available()}')
"
```

## Registry Authentication

Push requires GitHub Container Registry authentication:

1. **Create Personal Access Token**
   - GitHub → Settings → Developer Settings → Personal Access Tokens
   - Scopes: `read:packages`, `write:packages`, `delete:packages`

2. **Login to Registry**
   ```bash
   export CR_PAT=your_token
   echo $CR_PAT | docker login ghcr.io -u <username> --password-stdin
   ```

3. **Build and Push**
   ```bash
   ./build.sh --push
   ```

## Version Management

Tag images by version for reproducibility:

```bash
# Release version
./build.sh --tag v1.0.0 --push

# Development version
./build.sh --tag dev-2024-02-18 --push

# Update latest tag
./build.sh --tag latest --push
```

Update `nextflow.config` to use specific versions:
```nextflow
container = "ghcr.io/johnsonlab-ic/sceqtl-flow:v1.0.0"  # Reproducible
```

## Maintenance

### Regular Updates

When updating pipeline dependencies:

1. Update `.devcontainer/Dockerfile`
2. Rebuild and test
3. Tag with new version
4. Update `nextflow.config` if deploying immediately

### Monitoring Image Usage

```bash
# Check image locally
docker image ls | grep sceqtl-flow

# Check registry
git hash-object .devcontainer/Dockerfile  # For versioning
```

## Next Steps

1. ✅ Build unified image: `cd Images && ./build.sh --push`
2. ✅ Update `nextflow.config` to use single container specification
3. ✅ Test pipeline with small dataset
4. ✅ Document any custom modifications
5. ✅ Tag with version number for releases

## References

- `.devcontainer/Dockerfile` - Unified image definition
- `Images/build.sh` - Build and deploy script
- `nextflow.config` - Pipeline configuration
- `R/genotype_functions`, `R/expression_functions` - R script locations
