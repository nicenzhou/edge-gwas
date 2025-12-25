# EDGE GWAS Implementation Guide

## Version 0.0.0 (Under Public Testing)
⚠️ **Note:** This package is currently under active development and public testing.

## Overview

EDGE-GWAS (Encoding Deviation Genotypic Effects GWAS) identifies nonadditive SNP effects using flexible genetic encoding, rather than assuming additive inheritance.

**Key Features:**
- Two-stage analysis: calculate alpha on training data, apply to test data
- Detects recessive, dominant, and over-dominant effects
- Handles binary and quantitative outcomes
- Support for PLINK format data

## Installation

```bash
# Upgrade pip first
pip install --upgrade pip

# Install edge-gwas
pip install git+https://github.com/nicenzhou/edge-gwas.git

# Verify (use python3 on Mac/Linux)
python3 -c "from edge_gwas import EDGEAnalysis; print('✓ Installed successfully')"
```

## Quick Start

### 1. Import and Initialize

```python
import pandas as pd
from edge_gwas import EDGEAnalysis, load_plink_data, prepare_phenotype_data

# Initialize
edge = EDGEAnalysis(
    outcome_type='binary',    # or 'continuous'
    maf_threshold=0.01,
    verbose=True
)
```

### 2. Load Data

```python
# Load PLINK format genotype data
genotype_data, variant_info = load_plink_data(
    plink_prefix='path/to/plink_files',
    chromosome=None  # None for all, or specify chromosome number
)

# Load phenotype data
# Index: sample IDs matching genotype data
# Columns: outcome + covariates
phenotype_df = pd.read_csv('phenotype_file.csv', index_col=0)

# Define outcome and covariates
outcome = 'disease'  # your outcome variable name
covariates = ['age', 'sex', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 
              'PC6', 'PC7', 'PC8', 'PC9', 'PC10']

# Prepare phenotype (QC and formatting)
phenotype_df = prepare_phenotype_data(
    phenotype_df,
    outcome=outcome,
    covariates=covariates,
    remove_outliers=True  # For continuous outcomes
)
```

### 3. Run Two-Stage Analysis

```python
# Split into training and test sets
train_samples = phenotype_df.sample(frac=0.5, random_state=42).index
test_samples = phenotype_df.index.difference(train_samples)

# Run complete analysis
alpha_df, gwas_df = edge.run_full_analysis(
    train_genotype=genotype_data.loc[train_samples],
    train_phenotype=phenotype_df.loc[train_samples],
    test_genotype=genotype_data.loc[test_samples],
    test_phenotype=phenotype_df.loc[test_samples],
    outcome=outcome,
    covariates=covariates,
    variant_info=variant_info,
    output_prefix='edge_results'
)

# Check results
print(f"Tested {len(gwas_df)} variants")
print(f"Significant (p < 5e-8): {(gwas_df['pval'] < 5e-8).sum()}")
```

## Key Output Files

**Alpha values** (`edge_results_alpha.txt`):
- `variant_id`: SNP ID
- `alpha_value`: Encoding parameter (β_het/β_hom)
- `eaf`: Effect allele frequency
- `coef_het`, `coef_hom`: Coefficients

**GWAS results** (`edge_results_gwas.txt`):
- `variant_id`: SNP ID
- `coef`: Effect coefficient
- `pval`: P-value
- `alpha_value`: Applied alpha
- `std_err`, `stat`: Statistics

## Visualization

### Manhattan Plot

```python
import matplotlib.pyplot as plt
import numpy as np

def manhattan_plot(gwas_df, output='manhattan.png'):
    gwas_df['-log10p'] = -np.log10(gwas_df['pval'])
    
    fig, ax = plt.subplots(figsize=(14, 6))
    colors = ['#1f77b4', '#ff7f0e']
    x_pos = 0
    
    for chrom in sorted(gwas_df['chr'].unique()):
        data = gwas_df[gwas_df['chr'] == chrom]
        ax.scatter(x_pos + np.arange(len(data)), data['-log10p'], 
                  c=colors[int(chrom) % 2], s=2, alpha=0.7)
        x_pos += len(data)
    
    ax.axhline(-np.log10(5e-8), color='red', linestyle='--', label='p=5e-8')
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('-log₁₀(p-value)')
    ax.set_title('EDGE GWAS Manhattan Plot')
    plt.savefig(output, dpi=300, bbox_inches='tight')

manhattan_plot(gwas_df)
```

### QQ Plot

```python
from scipy import stats

def qq_plot(gwas_df, output='qq_plot.png'):
    pvals = gwas_df['pval'].dropna()
    pvals = pvals[pvals > 0]
    
    n = len(pvals)
    observed = -np.log10(np.sort(pvals))
    expected = -np.log10(np.arange(1, n + 1) / (n + 1))
    
    # Genomic inflation factor
    chisq = stats.chi2.ppf(1 - pvals, df=1)
    lambda_gc = np.median(chisq) / stats.chi2.ppf(0.5, df=1)
    
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(expected, observed, s=10, alpha=0.5)
    max_val = max(expected.max(), observed.max())
    ax.plot([0, max_val], [0, max_val], 'r--', linewidth=2)
    ax.text(0.05, 0.95, f'λ = {lambda_gc:.3f}', transform=ax.transAxes,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax.set_xlabel('Expected -log₁₀(p)')
    ax.set_ylabel('Observed -log₁₀(p)')
    plt.savefig(output, dpi=300, bbox_inches='tight')

qq_plot(gwas_df)
```

## Citation

Zhou, J., et al. (2023). Flexibly encoded genome-wide association study identifies novel nonadditive genetic risk variants for cardiometabolic traits. *medRxiv*, 2023.06.01.23290857. https://doi.org/10.1101/2023.06.01.23290857

```bibtex
@article{zhou2023edgegwas,
  title={Flexibly encoded genome-wide association study identifies novel nonadditive genetic risk variants for cardiometabolic traits},
  author={Zhou, Jiayan and Rico, Andre Luis Garao and Guare, Lindsay and others},
  journal={medRxiv},
  year={2023},
  doi={10.1101/2023.06.01.23290857}
}
```

## Contact

**Corresponding Author:** Molly Ann Hall - molly.hall@pennmedicine.upenn.edu  
**Code Questions:** Jiayan Zhou - jyzhou@stanford.edu

## License

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
