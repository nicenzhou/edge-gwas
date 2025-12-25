# edge-gwas: A Python package for identifying nonadditive SNP effects using flexible genetic encoding

[![Documentation Status](https://readthedocs.org/projects/edge-gwas/badge/?version=latest)](https://edge-gwas.readthedocs.io/en/latest/?badge=latest)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)

⚠️ **Current Version 0.1.0 - Under Public Testing**

## Overview

EDGE GWAS (Elastic Data-Driven Encoding GWAS) identifies nonadditive genetic effects (recessive, dominant, over-dominant) using a **two-stage data-driven approach** rather than assuming additive inheritance.  

## Quick Start

### Installation

```bash
pip install git+https://github.com/nicenzhou/edge-gwas.git
```

### Minimal Example

```python
from edge_gwas import EDGEAnalysis
from edge_gwas.utils import load_plink_data, prepare_phenotype_data, stratified_train_test_split
from edge_gwas.visualize import manhattan_plot, qq_plot

# Load data
geno, info = load_plink_data('data.bed', 'data.bim', 'data.fam')
pheno = prepare_phenotype_data('pheno.txt', 'disease', ['age', 'sex', 'PC1', 'PC2'])

# Split data
train_g, test_g, train_p, test_p = stratified_train_test_split(geno, pheno, 'disease')

# Run EDGE analysis
edge = EDGEAnalysis(outcome_type='binary')
alpha_df, gwas_df = edge.run_full_analysis(
    train_g, train_p, test_g, test_p,
    outcome='disease', 
    covariates=['age', 'sex', 'PC1', 'PC2']
)

# Visualize
manhattan_plot(gwas_df, 'manhattan.png')
qq_plot(gwas_df, 'qq.png')
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
- $COV_i$ = covariates
- $\beta_{Het}$ = effect size for heterozygous genotype
- $\beta_{HA}$ = effect size for homozygous alternate genotype
- $\beta_{cov_i}$ = effect sizes for covariates
  
**Equation 2: Encoding Parameter**

$$\alpha = \frac{\beta_{Het}}{\beta_{HA}}$$

Where:
- $$\alpha$$ = encoding parameter representing the ratio of heterozygous to homozygous alternate effects
- $$\beta_{Het}$$ = effect size for heterozygous genotype (from Equation 1)
- $$\beta_{HA}$$ = effect size for homozygous alternate genotype (from Equation 1)

**Learn more:** https://edge-gwas.readthedocs.io/en/latest/statistical_model.html

## Documentation

**Complete documentation available at:** https://edge-gwas.readthedocs.io/

- [Installation Guide](https://edge-gwas.readthedocs.io/en/latest/installation.html)
- [Quick Start Tutorial](https://edge-gwas.readthedocs.io/en/latest/quickstart.html)
- [User Guide](https://edge-gwas.readthedocs.io/en/latest/user_guide.html)
- [API Reference](https://edge-gwas.readthedocs.io/en/latest/api_reference.html)
- [Example Workflows](https://edge-gwas.readthedocs.io/en/latest/examples.html)
- [Statistical Model](https://edge-gwas.readthedocs.io/en/latest/statistical_model.html)
  
## Version History

### v0.1.0 (Current - Public Testing)
- Initial packaged release
- Complete two-stage EDGE analysis
- Support for binary and quantitative outcomes
- Comprehensive documentation on Read the Docs

### v0.0.0 (Deprecated)
- Original standalone scripts
- **No longer maintained**
- Available at: https://github.com/nicenzhou/EDGE (reference only)

**Migrate to v0.1.0 for continued support**

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

## Contact

**Corresponding Author:**  
Molly Ann Hall - molly.hall@pennmedicine.upenn.edu

**For questions about the code:**  
Jiayan Zhou - jyzhou@stanford.edu

## License
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

[![GPL Logo](https://www.gnu.org/graphics/gplv3-88x31.png)](https://www.gnu.org/licenses/gpl-3.0)
