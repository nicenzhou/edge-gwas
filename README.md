# edge-gwas: A Python package for identifying nonadditive SNP effects using flexible genetic encoding

[![Version](https://img.shields.io/badge/version-0.1.1-green.svg)](https://github.com/nicenzhou/edge-gwas/releases)
[![Documentation Status](https://readthedocs.org/projects/edge-gwas/badge/?version=v0.1.1)](https://edge-gwas.readthedocs.io/en/v0.1.1/?badge=v0.1.1)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)

⚠️ **Current Version 0.1.1 - Under Public Testing**

**Recommended to use v0.1.1 - more stable and more functions.**

## Overview

EDGE GWAS (Elastic Data-Driven Encoding GWAS) identifies nonadditive genetic effects (recessive, dominant, over-dominant) using a **two-stage data-driven approach** rather than assuming additive inheritance.

### Key Features

**Core Functionality:**
- **Flexible genetic encoding**: Detects recessive, additive, dominant, and over-dominant effects
- **Two-stage design**: Reduces overfitting through train/test split
- **Binary and continuous outcomes**: With outcome transformations (log, inverse normal)

**Population Structure Control (NEW in v0.1.1):**
- **PCA calculation**: PLINK2 (exact/approximate), PC-AiR (relatedness-aware), sklearn
- **GRM support**: Calculate and use genetic relationship matrices
- **Mixed models**: Linear and logistic mixed models for related samples

**Data Handling:**
- **Multiple formats**: PLINK (.bed), PLINK2 (.pgen), BGEN, VCF
- **Comprehensive QC**: MAF, missingness, HWE, sample call rate filters
- **Cross-validation**: K-fold CV for alpha stability assessment

**Visualization & Tools:**
- Manhattan plots, QQ plots, alpha distribution plots
- Automated tool installation (`edge-gwas-install-tools`)
- Installation verification (`edge-gwas-check-tools`)

**Performance:**
- Parallel processing with multi-core support
- Memory-efficient handling of biobank-scale data
- Approximate PCA for large cohorts (>5K samples)

## Quick Start

### Installation

```bash
# Install edge-gwas
pip install git+https://github.com/nicenzhou/edge-gwas.git

# Install external tools (PLINK2, GCTA, R packages)
edge-gwas-install-tools

# Verify installation
edge-gwas-check-tools
```

**Supported Platforms:** Linux, macOS (Intel & Apple Silicon)

See [Installation Guide](https://edge-gwas.readthedocs.io/en/v0.1.1/installation.html) for details.

### Fast minimal EDGE GWAS analysis

```python
from edge_gwas import EDGEAnalysis
from edge_gwas.utils import (
    load_plink_data,
    prepare_phenotype_data,
    filter_variants_by_maf,
    calculate_pca_plink,
    attach_pcs_to_phenotype,
    get_pc_covariate_list,
    stratified_train_test_split
)
from edge_gwas.visualize import manhattan_plot, qq_plot

# 1. Load data
geno, info = load_plink_data('data.bed', 'data.bim', 'data.fam')
pheno = prepare_phenotype_data('pheno.txt', 'disease', ['age', 'sex'])

# 2. QC filtering
geno = filter_variants_by_maf(geno, min_maf=0.01)

# 3. Calculate PCA for population structure
pca_df = calculate_pca_plink('data', n_pcs=10)
pheno = attach_pcs_to_phenotype(pheno, pca_df, n_pcs=10)

# 4. Split data
train_g, test_g, train_p, test_p = stratified_train_test_split(
    geno, pheno, 'disease', test_size=0.5
)

# 5. Run EDGE analysis
edge = EDGEAnalysis(outcome_type='binary', n_jobs=-1)
covariates = ['age', 'sex'] + get_pc_covariate_list(10)

alpha_df, gwas_df = edge.run_full_analysis(
    train_g, train_p, test_g, test_p,
    outcome='disease',
    covariates=covariates
)

# 6. Visualize
manhattan_plot(gwas_df, 'manhattan.png')
lambda_gc = qq_plot(gwas_df, 'qq.png')

print(f"Significant hits: {(gwas_df['pval'] < 5e-8).sum()}")
print(f"Lambda GC: {lambda_gc:.3f}")
```

### Advanced Example: Related Samples with GRM

```python
from edge_gwas.utils import calculate_grm_gcta, load_grm_gcta, calculate_pca_pcair

# Calculate GRM
grm_prefix = calculate_grm_gcta('data', maf_threshold=0.01)
grm_matrix, grm_ids = load_grm_gcta(grm_prefix)

# Calculate PC-AiR (accounts for relatedness)
pca_df = calculate_pca_pcair('data', n_pcs=10, kinship_matrix=grm_prefix)
pheno = attach_pcs_to_phenotype(pheno, pca_df, n_pcs=10)

# Run EDGE with GRM
edge = EDGEAnalysis(outcome_type='binary')
alpha_df, gwas_df = edge.run_full_analysis(
    train_g, train_p, test_g, test_p,
    outcome='disease',
    covariates=covariates,
    grm_matrix=grm_matrix,
    grm_sample_ids=grm_ids
)
```

**More examples:** [Complete Workflows](https://edge-gwas.readthedocs.io/en/latest/examples.html)

## Statistical Model

EDGE uses flexible encoding based on separate heterozygous and homozygous effects:

**Regression Model:**

$$E(Y | G, COV) = \beta_0 + \beta_{Het} \cdot G_{Het} + \beta_{Hom} \cdot G_{Hom} + \sum_{i} \beta_{cov_i} \cdot COV_i$$

**Encoding Parameter:**

$$\alpha = \frac{\beta_{Het}}{\beta_{Hom}}$$

**Interpretation:**
- **α ≈ 0**: Recessive (only homozygotes affected)
- **α ≈ 0.5**: Additive (heterozygotes intermediate)
- **α ≈ 1**: Dominant (heterozygotes = homozygotes)
- **α < 0 or α > 1**: Over-dominance/under-recessiveness

**Learn more:** [Statistical Model Documentation](https://edge-gwas.readthedocs.io/en/latest/statistical_model.html)

## What's New in v0.1.1

### Major Changes

**Breaking Changes:**
- Removed Koalas dependency → Pure pandas (better compatibility)
- PCA functions now return DataFrames with IID as index
- All sample IDs converted to strings for consistency

**New Features:**

**Population Structure Control:**
- `calculate_pca_plink()` - PCA with PLINK2 (exact/approximate)
- `calculate_pca_pcair()` - PC-AiR for related samples  
- `calculate_pca_sklearn()` - Basic PCA
- `calculate_grm_gcta()` - GRM calculation
- `load_grm_gcta()` - Load GRM from files
- `attach_pcs_to_phenotype()` - Merge PCs with phenotypes
- `identify_related_samples()` - Find related pairs
- `filter_related_samples()` - Remove related samples

**Outcome Transformations:**
- `outcome_transform` parameter: 'log', 'log10', 'inverse_normal', 'rank_inverse_normal'

**File Format Support:**
- `load_pgen_data()` - PLINK2 format
- `load_bgen_data()` - BGEN format (UK Biobank)
- `load_vcf_data()` - VCF/VCF.GZ format

**Enhanced QC:**
- `filter_variants_by_hwe()` - Hardy-Weinberg equilibrium filtering
- `filter_samples_by_call_rate()` - Sample quality filtering
- `check_case_control_balance()` - Check case/control ratio
- `calculate_hwe_pvalues()` - Calculate HWE p-values

**Analysis Tools:**
- `additive_gwas()` - Standard additive GWAS for comparison
- `cross_validated_edge_analysis()` - K-fold cross-validation

**Command-Line Tools:**
- `edge-gwas-install-tools` - Interactive tool installer
- `edge-gwas-check-tools` - Verify tool installation

## Documentation

**Complete documentation:** https://edge-gwas.readthedocs.io/

- [Installation Guide](https://edge-gwas.readthedocs.io/en/latest/installation.html) - Installation instructions and requirements
- [Quick Start Tutorial](https://edge-gwas.readthedocs.io/en/latest/quickstart.html) - Getting started with EDGE in 5 minutes
- [User Guide](https://edge-gwas.readthedocs.io/en/latest/user_guide.html) - Comprehensive user guide and tutorials
- [API Reference](https://edge-gwas.readthedocs.io/en/latest/api_reference.html) - Complete API documentation
- [Example Workflows](https://edge-gwas.readthedocs.io/en/latest/examples.html) - Real-world analysis examples
- [Statistical Model](https://edge-gwas.readthedocs.io/en/latest/statistical_model.html) - Mathematical background and methods
- [Visualization Guide](https://edge-gwas.readthedocs.io/en/latest/visualization.html) - Creating publication-ready plots
- [Troubleshooting](https://edge-gwas.readthedocs.io/en/latest/troubleshooting.html) - Common issues and solutions
- [FAQ](https://edge-gwas.readthedocs.io/en/latest/faq.html) - Frequently asked questions
- [Citation Guide](https://edge-gwas.readthedocs.io/en/latest/citation.html) - How to cite EDGE
- [Changelog](https://edge-gwas.readthedocs.io/en/latest/changelog.html) - Version history and updates
- [Future Updates](https://edge-gwas.readthedocs.io/en/latest/futureupdates.html) - Planned features and roadmap

## System Requirements

**Minimum:**
- Python 3.7+
- 4 GB RAM
- 1 GB disk space

**Recommended:**
- Python 3.9+
- 16 GB RAM (for large datasets)
- 10 GB disk space
- Multi-core CPU

**For Biobank-Scale (>100K samples):**
- 32+ GB RAM
- Use approximate PCA: `calculate_pca_plink(approx=True)`
- Process chromosomes separately
- HPC cluster recommended

## Python Requirements

**Core Dependencies:**
- numpy >= 1.19.0
- pandas >= 1.2.0 (Koalas removed in v0.1.1)
- scipy >= 1.6.0
- statsmodels >= 0.12.0
- scikit-learn >= 0.24.0
- matplotlib >= 3.3.0
- pandas-plink >= 2.0.0
- joblib >= 1.0.0

**Optional File Format Support:**
- pgenlib >= 0.81.0 (PLINK2 format)
- bgen-reader >= 4.0.8 (BGEN format)
- cyvcf2 >= 0.30.0 (VCF format)

Install optional dependencies:

```bash
pip install pgenlib bgen-reader cyvcf2
```

## External Tools

EDGE-GWAS integrates with external tools (installed automatically via `edge-gwas-install-tools`):

**PLINK2:**
- Fast PCA calculation and LD pruning
- Download: https://www.cog-genomics.org/plink/2.0/

**GCTA:**
- Genetic relationship matrix calculation
- Download: https://yanglab.westlake.edu.cn/software/gcta/

**R + GENESIS:**
- PC-AiR for relatedness-aware PCA
- R packages: GENESIS, SNPRelate, gdsfmt

**Verify installation:**

```bash
edge-gwas-check-tools
```

## Citation

If you use edge-gwas in your research, please cite:

```bibtex
@article{zhou2023edgegwas,
  title={Flexibly encoded genome-wide association study identifies novel nonadditive 
         genetic risk variants for cardiometabolic traits},
  author={Zhou, Jiayan and Rico, Andre Luis Garao and Guare, Lindsay and others},
  journal={medRxiv},
  year={2023},
  doi={10.1101/2023.06.01.23290857}
}
```

## Version History

### v0.1.1 (2025-12-25 - Current)
**Breaking Changes:**
- Removed Koalas dependency → Pure pandas
- PCA functions now return DataFrames with IID index
- Sample ID handling standardized (all converted to strings)

**New Features:**
- GRM support for mixed models (linear and logistic)
- Multiple PCA methods (PLINK2, PC-AiR, sklearn)
- Outcome transformations for continuous traits
- Support for PGEN, BGEN, VCF formats
- Enhanced QC functions (HWE, call rate filtering)
- Cross-validation and additive GWAS comparison
- Automated tool installation system

### v0.1.0 (2025-12-24 - Deprecated)
- Initial packaged release
- Two-stage EDGE analysis framework
- Support for binary and continuous outcomes
- PLINK format support
- Basic visualization
- Complete documentation

### v0.0.0 (2024-04-02 - Deprecated)
- Original standalone scripts
- No longer maintained
- Available at: https://github.com/nicenzhou/EDGE (reference only)

**Migrate to v0.1.1 for continued support**

## Support

**Contact:**
- **Code questions:** Jiayan Zhou - jyzhou@stanford.edu
- **Research questions:** Molly Ann Hall - molly.hall@pennmedicine.upenn.edu

## License

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

[![GPL Logo](https://www.gnu.org/graphics/gplv3-88x31.png)](https://www.gnu.org/licenses/gpl-3.0)

---

*Last updated: 2025-12-25 for edge-gwas v0.1.1*

*For questions or issues, visit:* https://github.com/nicenzhou/edge-gwas/issues
