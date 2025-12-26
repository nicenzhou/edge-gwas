# edge-gwas: A Python package for identifying nonadditive SNP effects using flexible genetic encoding

[![Documentation Status](https://readthedocs.org/projects/edge-gwas/badge/?version=latest)](https://edge-gwas.readthedocs.io/en/latest/?badge=latest)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)

⚠️ **Current Version 0.1.1 - Under Public Testing**

## Overview

EDGE GWAS (Elastic Data-Driven Encoding GWAS) identifies nonadditive genetic effects (recessive, dominant, over-dominant) using a **two-stage data-driven approach** rather than assuming additive inheritance.

### Key Features

- **Flexible genetic encoding**: Detects recessive, additive, dominant, and over-dominant effects
- **Two-stage design**: Reduces overfitting through train/test split
- **Population structure support**: 
  - PCA calculation with PLINK2 (exact and approximate methods)
  - PC-AiR for relatedness-aware PCA
  - GRM calculation with GCTA
- **Multiple file formats**: PLINK (.bed), PGEN, BGEN, VCF
- **Comprehensive QC**: Built-in filters for MAF, missingness, HWE
- **Rich visualization**: Manhattan plots, QQ plots, alpha distribution plots
- **Parallel processing**: Multi-core support for large-scale analyses

## Quick Start

### Installation

**Quick install (recommended):**

```bash
# Clone and install with all dependencies
git clone https://github.com/nicenzhou/edge-gwas.git
cd edge-gwas
chmod +x install.sh
./install.sh
```

**Manual install:**

```bash
# Install Python package
pip install git+https://github.com/nicenzhou/edge-gwas.git

# Install external tools (PLINK2, GCTA, R packages)
edge-gwas-install-tools

# Verify installation
edge-gwas-check-tools
```

See [INSTALL.md](INSTALL.md) for detailed installation instructions.

### Minimal Example

```python
from edge_gwas import EDGEAnalysis
from edge_gwas.utils import (
    load_plink_data, 
    prepare_phenotype_data, 
    stratified_train_test_split,
    calculate_pca_plink,
    attach_pcs_to_phenotype
)
from edge_gwas.visualization import manhattan_plot, qq_plot

# Load genotype data
geno, info = load_plink_data('data.bed', 'data.bim', 'data.fam')

# Calculate PCs for population structure
pca_df = calculate_pca_plink('data', n_pcs=10, approx=True)

# Load and prepare phenotype data
pheno = prepare_phenotype_data('pheno.txt', 'disease', ['age', 'sex'])

# Attach PCs to phenotype
pheno = attach_pcs_to_phenotype(pheno, pca_df, n_pcs=10)

# Split data
train_g, test_g, train_p, test_p = stratified_train_test_split(
    geno, pheno, 'disease', test_size=0.5
)

# Run EDGE analysis
edge = EDGEAnalysis(outcome_type='binary', n_jobs=8)
alpha_df, gwas_df = edge.run_full_analysis(
    train_g, train_p, test_g, test_p,
    outcome='disease', 
    covariates=['age', 'sex', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']
)

# Visualize results
manhattan_plot(gwas_df, output_file='manhattan.png')
qq_plot(gwas_df, output_file='qq.png')
```

**For complete examples and detailed workflows, see:** https://edge-gwas.readthedocs.io/en/latest/examples.html

## Statistical Model

EDGE employs a flexible encoding approach based on a regression model that separately estimates effects for heterozygous and homozygous alternate genotypes:

**Equation 1: Regression Model**

$$E(Y | SNP_{Het}, SNP_{HA}, COV_i) = \beta_0 + \beta_{Het} \cdot SNP_{Het} + \beta_{HA} \cdot SNP_{HA} + \sum_{i} \beta_{cov_i} \cdot COV_i$$

Where:
- $Y$ = phenotype/outcome
- $SNP_{Het}$ = indicator for heterozygous genotype
- $SNP_{HA}$ = indicator for homozygous alternate genotype
- $COV_i$ = covariates (including population structure PCs)
- $\beta_{Het}$ = effect size for heterozygous genotype
- $\beta_{HA}$ = effect size for homozygous alternate genotype
- $\beta_{cov_i}$ = effect sizes for covariates
  
**Equation 2: Encoding Parameter**

$$\alpha = \frac{\beta_{Het}}{\beta_{HA}}$$

Where:
- $\alpha$ = encoding parameter representing the ratio of heterozygous to homozygous alternate effects
- Interpretation:
  - $\alpha \approx 0$: Recessive model
  - $\alpha \approx 0.5$: Additive model
  - $\alpha \approx 1$: Dominant model
  - Other values: Over-dominant or under-dominant effects

**Learn more:** https://edge-gwas.readthedocs.io/en/latest/statistical_model.html

## Population Structure Correction

EDGE-GWAS provides multiple methods for population structure correction:

### 1. Basic PCA (sklearn)
Fast, simple PCA without relatedness correction:
```python
from edge_gwas.utils import calculate_pca_sklearn
pca_df = calculate_pca_sklearn(genotype_df, n_pcs=10)
```

### 2. PLINK2 PCA
Fast, LD-pruned PCA with optional approximate method for large cohorts:
```python
from edge_gwas.utils import calculate_pca_plink

# Exact PCA (smaller cohorts)
pca_df = calculate_pca_plink('data', n_pcs=10, approx=False)

# Approximate PCA (large cohorts >5000 samples)
pca_df = calculate_pca_plink('data', n_pcs=10, approx=True, approx_samples=5000)
```

### 3. PC-AiR
Relatedness-aware PCA using GENESIS:
```python
from edge_gwas.utils import calculate_pca_pcair, calculate_grm_gcta

# Calculate GRM
grm_prefix = calculate_grm_gcta('data')

# Run PC-AiR
pca_df = calculate_pca_pcair('data', n_pcs=10, kinship_matrix=grm_prefix)
```

### 4. GRM for Mixed Models
Calculate genetic relationship matrix:
```python
from edge_gwas.utils import calculate_grm_gcta, load_grm_gcta

# Calculate GRM
grm_prefix = calculate_grm_gcta('data', maf_threshold=0.01)

# Load GRM
grm_matrix, sample_ids = load_grm_gcta(grm_prefix)
```

## Documentation

**Complete documentation available at:** https://edge-gwas.readthedocs.io/

- [Installation Guide](https://edge-gwas.readthedocs.io/en/latest/installation.html)
- [Quick Start Tutorial](https://edge-gwas.readthedocs.io/en/latest/quickstart.html)
- [User Guide](https://edge-gwas.readthedocs.io/en/latest/user_guide.html)
- [API Reference](https://edge-gwas.readthedocs.io/en/latest/api_reference.html)
- [Example Workflows](https://edge-gwas.readthedocs.io/en/latest/examples.html)
- [Statistical Model](https://edge-gwas.readthedocs.io/en/latest/statistical_model.html)
- [Population Structure](https://edge-gwas.readthedocs.io/en/latest/population_structure.html)

## External Tools

EDGE-GWAS integrates with several external tools (installed automatically):

- **PLINK2**: Fast PCA calculation and data QC
  - Download: https://www.cog-genomics.org/plink/2.0/
- **GCTA**: Genetic relationship matrix calculation
  - Download: https://yanglab.westlake.edu.cn/software/gcta/
- **R + GENESIS**: PC-AiR for relatedness-aware PCA
  - R packages: GENESIS, SNPRelate, gdsfmt

These are installed automatically during package installation. To verify:
```bash
edge-gwas-check-tools
```

## Version History

### v0.1.1 (2025-12-25 - Current)
- **Breaking change**: Replaced Koalas with pandas for better compatibility
- Added population structure functions:
  - `calculate_pca_plink()` with exact and approximate methods
  - `calculate_pca_pcair()` for relatedness-aware PCA
  - `calculate_grm_gcta()` for genetic relationship matrices
  - `attach_pcs_to_phenotype()` helper function
  - `identify_related_samples()` and `filter_related_samples()`
- Automatic installation of external tools (PLINK2, GCTA, R packages)
- Support for multiple genetic file formats (PGEN, BGEN, VCF)
- Enhanced QC functions (HWE, call rate filtering)
- Command-line tools: `edge-gwas-check-tools`, `edge-gwas-install-tools`

### v0.1.0 (2025-12-24)
- Initial packaged release
- Complete two-stage EDGE analysis
- Support for binary and quantitative outcomes
- Comprehensive documentation on Read the Docs

### v0.0.0 (2024-04-02 - Deprecated)
- Original standalone scripts
- **No longer maintained**
- Available at: https://github.com/nicenzhou/EDGE (reference only)

**Migrate to v0.1.1 for continued support**

## Requirements

### Python Requirements
- Python 3.7+
- numpy >= 1.19.0
- pandas >= 1.2.0
- scipy >= 1.6.0
- statsmodels >= 0.12.0
- scikit-learn >= 0.24.0
- matplotlib >= 3.3.0
- pandas-plink >= 2.0.0

### External Tools (Installed Automatically)
- PLINK2 (for PCA calculation)
- GCTA (for GRM calculation)
- R with GENESIS, SNPRelate, gdsfmt (for PC-AiR)

### Optional Dependencies
- pgenlib (for PGEN format)
- bgen-reader (for BGEN format)
- cyvcf2 (for VCF format)

See [requirements.txt](requirements.txt) for complete list.

## System Requirements

**Minimum:**
- 4 GB RAM
- 1 GB disk space

**Recommended:**
- 16 GB RAM (for large datasets)
- 10 GB disk space
- Multi-core CPU

**For Large-Scale Analyses (>100K samples):**
- 32+ GB RAM
- Use approximate PCA: `calculate_pca_plink(approx=True)`
- HPC cluster recommended

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

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Areas for Contribution
- Additional file format support
- Performance optimizations
- Enhanced visualization options
- Documentation improvements
- Bug reports and fixes

## Support

- **Documentation**: https://edge-gwas.readthedocs.io/
- **GitHub Issues**: https://github.com/nicenzhou/edge-gwas/issues
- **Discussions**: https://github.com/nicenzhou/edge-gwas/discussions

## Contact

**Corresponding Author:**  
Molly Ann Hall - molly.hall@pennmedicine.upenn.edu

**For questions about the code:**  
Jiayan Zhou - jyzhou@stanford.edu

## License

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

[![GPL Logo](https://www.gnu.org/graphics/gplv3-88x31.png)](https://www.gnu.org/licenses/gpl-3.0)

## Acknowledgments

- Original EDGE method development: Jiayan Zhou, Molly Hall, and collaborators
- Built with support from pandas-plink, scikit-learn, and statsmodels
- External tools: PLINK2 (Chang et al.), GCTA (Yang et al.), GENESIS (Gogarten et al.)
