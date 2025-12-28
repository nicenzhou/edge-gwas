![Logo](docs/images/EDGE_RG.jpg)

# edge-gwas: A Python Package for Flexible Genetic Encoding in GWAS

[![Version](https://img.shields.io/badge/version-0.1.1-green.svg)](https://github.com/nicenzhou/edge-gwas/releases)
[![Documentation](https://readthedocs.org/projects/edge-gwas/badge/?version=latest)](https://edge-gwas.readthedocs.io)
[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**Identify nonadditive genetic effects (recessive, dominant, over-dominant) using data-driven encoding.**

**📖 Full Documentation:** https://edge-gwas.readthedocs.io

---

## Overview

EDGE (Elastic Data-Driven Genetic Encoding) discovers nonadditive SNP effects by learning optimal genetic encodings from data, rather than assuming additive inheritance.

**Key Advantages:**
- **Flexible inheritance models**: Detects recessive, dominant, and over-dominant effects
- **Two-stage design**: Prevents overfitting via independent train/test sets
- **Population structure control**: GRM-based mixed models, PC-AiR for relatives
- **Scalable**: Handles biobank-scale data with parallel processing

---

## Quick Start

### Installation

```bash
# Install latest stable version
pip install git+https://github.com/nicenzhou/edge-gwas.git@v0.1.1

# Install external tools (PLINK2, GCTA, R packages)
edge-gwas-install-tools

# Verify installation
edge-gwas-check-tools
```

**Platforms:** Linux, macOS (Intel & Apple Silicon)

---

### Basic Example

```python
from edge_gwas import EDGEAnalysis
from edge_gwas.utils import (
    load_plink_data, prepare_phenotype_data,
    filter_variants_by_maf, calculate_pca_plink,
    attach_pcs_to_phenotype, get_pc_covariate_list,
    stratified_train_test_split
)

# 1. Load and QC data
geno, info = load_plink_data('data.bed', 'data.bim', 'data.fam')
geno = filter_variants_by_maf(geno, min_maf=0.01)
pheno = prepare_phenotype_data('pheno.txt', 'disease', ['age', 'sex'])

# 2. Control population structure
pca_df = calculate_pca_plink('data', n_pcs=10)
pheno = attach_pcs_to_phenotype(pheno, pca_df, n_pcs=10)
covariates = ['age', 'sex'] + get_pc_covariate_list(10)

# 3. Split data
train_g, test_g, train_p, test_p = stratified_train_test_split(
    geno, pheno, 'disease', test_size=0.5
)

# 4. Run EDGE
edge = EDGEAnalysis(outcome_type='binary', n_jobs=-1)
alpha_df, gwas_df = edge.run_full_analysis(
    train_g, train_p, test_g, test_p,
    outcome='disease', covariates=covariates
)

# 5. Results
print(f"Genome-wide significant hits: {(gwas_df['pval'] < 5e-8).sum()}")
print(f"Mean alpha: {alpha_df['alpha_value'].mean():.3f}")
```

---

### Full Detailed Example

For a comprehensive example covering all features including:
- Data validation and QC filtering (MAF, HWE, missingness)
- Covariate imputation and missing data handling
- GRM calculation and mixed models for related samples
- PCA with LD pruning for population structure control
- Advanced visualization and multiple output formats

**See the full detailed workflow in the documentation:** [Complete Workflows Guide](https://edge-gwas.readthedocs.io/en/latest/quickstart.html#complete-workflow-example)

---

## The EDGE Model

### Two-Stage Approach

**Stage 1: Estimate α (Training Set)**

Fit codominant model to estimate separate effects:

$$E(Y | G, COV) = \beta_0 + \beta_{Het} \cdot I_{Het} + \beta_{HA} \cdot I_{HA} + \sum_{i} \beta_{cov_i} \cdot COV_i$$

Calculate encoding parameter:

$$\alpha = \frac{\beta_{Het}}{\beta_{HA}}$$

**Stage 2: GWAS (Test Set)**

Apply learned encoding and test association:

$$E(Y | G, COV) = \beta_0 + \beta_{EDGE} \cdot G_{EDGE} + \sum_{i} \beta_{cov_i} \cdot COV_i$$

where:

$$G_{EDGE} = \begin{cases} 
0 & \text{if } G = 0 \text{ (REF/REF)} \\
\alpha & \text{if } G = α \text{ (REF/ALT)} \\
1 & \text{if } G = 1 \text{ (ALT/ALT)}
\end{cases}$$

where α is estimated from data:
- **α ≈ 0**: Recessive effect (only homozygotes affected)
- **α ≈ 0.5**: Additive effect (heterozygotes intermediate)
- **α ≈ 1**: Dominant effect (heterozygotes ≈ homozygotes)
- **α < 0 or α > 1**: Over-dominance/under-dominance

### Population Structure Control

For related samples, EDGE uses mixed models:

$$Y = X\beta + Zu + \epsilon$$

where:
- $u \sim N(0, \sigma_g^2 \cdot GRM)$ (genetic random effects)
- $\epsilon \sim N(0, \sigma_e^2 \cdot I)$ (residual errors)
- $GRM$ is the genetic relationship matrix

---

## What's New in v0.1.1

**New Features:**
- GRM support for mixed models
- Multiple PCA methods (PLINK2, PC-AiR, sklearn)
- Outcome transformations for continuous traits
- PGEN, BGEN, VCF format support
- Automated tool installation (`edge-gwas-install-tools`)
- Cross-validation and additive GWAS comparison

---

## System Requirements

| Component | Minimum | Recommended | Biobank-Scale |
|-----------|---------|-------------|---------------|
| **Python** | 3.7+ | 3.9+ | 3.9+ |
| **RAM** | 4 GB | 16 GB | 32+ GB |
| **Cores** | 1 | 4+ | 16+ |
| **Disk** | 1 GB | 10 GB | 100+ GB |

**Tips for large datasets:**
- Use approximate PCA: `calculate_pca_plink(approx=True)`
- Process chromosomes separately
- Use HPC cluster with parallel processing

---

## Dependencies

**Core packages (auto-installed):**
```
numpy, pandas, scipy, statsmodels, scikit-learn, 
matplotlib, pandas-plink, joblib
```

**Optional (for file formats):**
```bash
pip install pgenlib bgen-reader cyvcf2
```

**External tools (via `edge-gwas-install-tools`):**
- PLINK2 (PCA, LD pruning)
- GCTA (GRM calculation)
- R + GENESIS (PC-AiR)

---

## Citation

If you use edge-gwas, please cite:

```bibtex
@article{zhou2023edge,
  title={Flexibly encoded genome-wide association study identifies novel 
         nonadditive genetic risk variants for cardiometabolic traits},
  author={Zhou, Jiayan and Rico, Andre Luis Garao and Guare, Lindsay and others},
  journal={medRxiv},
  year={2023},
  doi={10.1101/2023.06.01.23290857}
}
```

**Paper:** [medRxiv](https://doi.org/10.1101/2023.06.01.23290857)

---

## Support & Contributing

**Documentation:** [FAQ](https://edge-gwas.readthedocs.io/en/latest/faq.html) | [Troubleshooting](https://edge-gwas.readthedocs.io/en/latest/troubleshooting.html)

**Questions?**
- Code: Jiayan Zhou (jyzhou@stanford.edu)
- Research: Molly Ann Hall (molly.hall@pennmedicine.upenn.edu)

**Found a bug?** [Open an issue](https://github.com/nicenzhou/edge-gwas/issues)

## License

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

[![GPL Logo](https://www.gnu.org/graphics/gplv3-88x31.png)](https://www.gnu.org/licenses/gpl-3.0)

---

*Last updated: 2025-12-28 for edge-gwas v0.1.1*

*For questions or issues, visit:* https://github.com/nicenzhou/edge-gwas/issues

**© 2022-2025 Jiayan Zhou, Molly Ann Hall, and Contributors**
