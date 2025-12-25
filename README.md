# EDGE GWAS Python Package Implementation Guide

## Version 0.1.0 (Under Public Testing)

⚠️ **Note:** This package is currently under active development and public testing.

## Overview

EDGE-GWAS (Elastic Data-Driven Encoding GWAS) identifies nonadditive SNP effects using flexible genetic encoding, rather than assuming additive inheritance.

**Key Features:**
- Two-stage analysis: calculate alpha on training data, apply to test data
- Detects under-recessive, recessive, additive, dominant, and over-dominant effects on a continuous scale
- Handles binary and quantitative outcomes
- Support for PLINK format data
- Built-in visualization functions

## Installation

```bash
# Install edge-gwas with all dependencies
pip install git+https://github.com/nicenzhou/edge-gwas.git
```

```bash
# Install edge-gwas on Mac/Linux
# Clone the repository
git clone https://github.com/nicenzhou/edge-gwas.git
cd edge-gwas

# Install dependencies
pip install -r requirements.txt

# Install the package in development mode
pip install -e .
```

```python
# Verify installation in Python/Jupyter Notebook
from edge_gwas import EDGEAnalysis, manhattan_plot, qq_plot
```

```python
# Verify installation use python3 on Mac/Linux
python3 -c "from edge_gwas import EDGEAnalysis, manhattan_plot, qq_plot; print('✓ Installed successfully')"
```

## Quick Start

### 1. Import and Initialize

```python
import pandas as pd
from edge_gwas import (
    EDGEAnalysis, 
    load_plink_data, 
    prepare_phenotype_data,
    manhattan_plot,
    qq_plot,
    plot_alpha_distribution
)

# Initialize
edge = EDGEAnalysis(
    outcome_type='binary', # or 'continuous'
    n_jobs = 8, # number of core
    max_iter=1000, # maximum iterations for model convergence
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

#### Option A: One-Step (Recommended for Quick Analysis)

```python
# Split into training and test sets
train_samples = phenotype_df.sample(frac=0.5, random_state=42).index
test_samples = phenotype_df.index.difference(train_samples)

# Run complete analysis in one step
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

print(f"Tested {len(gwas_df)} variants")
print(f"Significant (p < 5e-8): {(gwas_df['pval'] < 5e-8).sum()}")
```

#### Option B: Two-Step (More Control and Flexibility)

```python
# Split data
train_samples = phenotype_df.sample(frac=0.5, random_state=42).index
test_samples = phenotype_df.index.difference(train_samples)

# Step 1: Calculate alpha values on training data
alpha_df = edge.calculate_alpha(
    genotype_data=genotype_data.loc[train_samples],
    phenotype_df=phenotype_df.loc[train_samples],
    outcome=outcome,
    covariates=covariates,
    variant_info=variant_info
)
alpha_df.to_csv('alpha_values.txt', sep='\t', index=False)
print(f"Alpha calculated for {len(alpha_df)} variants")

# Step 2: Apply alpha values on test data
gwas_df = edge.apply_alpha(
    genotype_data=genotype_data.loc[test_samples],
    phenotype_df=phenotype_df.loc[test_samples],
    outcome=outcome,
    covariates=covariates,
    alpha_values=alpha_df
)
gwas_df.to_csv('gwas_results.txt', sep='\t', index=False)
print(f"Significant (p < 5e-8): {(gwas_df['pval'] < 5e-8).sum()}")
```

#### Option C: Use Pre-calculated Alpha Values

```python
# Load previously calculated alpha values
alpha_df = pd.read_csv('alpha_values.txt', sep='\t')

# Apply to new test data
gwas_df = edge.apply_alpha(
    genotype_data=new_genotype_data,
    phenotype_df=new_phenotype_df,
    outcome=outcome,
    covariates=covariates,
    alpha_values=alpha_df
)
```

### 4. Create Visualizations

```python
# Manhattan plot
manhattan_plot(
    gwas_df, 
    output='manhattan.png',
    title='My EDGE GWAS Study',
    sig_threshold=5e-8
)

# QQ plot with genomic inflation factor
lambda_gc = qq_plot(
    gwas_df, 
    output='qq_plot.png',
    title='EDGE GWAS QQ Plot'
)

# Alpha distribution plot
plot_alpha_distribution(
    alpha_df,
    output='alpha_distribution.png',
    bins=50
)

print(f"Genomic inflation factor (λ): {lambda_gc:.3f}")
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

## Visualization Functions

### Using Built-in Functions (Recommended)

```python
from edge_gwas import manhattan_plot, qq_plot, plot_alpha_distribution

# Manhattan plot
manhattan_plot(
    gwas_df,                    # GWAS results DataFrame
    output='manhattan.png',     # Output filename
    title='EDGE GWAS',          # Plot title
    sig_threshold=5e-8,         # Significance threshold
    figsize=(14, 6),            # Figure size
    colors=['#1f77b4', '#ff7f0e']  # Chromosome colors
)

# QQ plot
lambda_gc = qq_plot(
    gwas_df,                    # GWAS results DataFrame
    output='qq_plot.png',       # Output filename
    title='QQ Plot',            # Plot title
    figsize=(8, 8)              # Figure size
)

# Alpha distribution
plot_alpha_distribution(
    alpha_df,                   # Alpha values DataFrame
    output='alpha_dist.png',    # Output filename
    bins=50,                    # Number of bins
    figsize=(10, 6)             # Figure size
)
```

### Multiple Chromosomes

```python
# Combine results from multiple chromosomes
gwas_chr1 = pd.read_csv('edge_chr1_gwas.txt', sep='\t')
gwas_chr2 = pd.read_csv('edge_chr2_gwas.txt', sep='\t')
gwas_chr3 = pd.read_csv('edge_chr3_gwas.txt', sep='\t')

# Pass as list to manhattan_plot
manhattan_plot(
    [gwas_chr1, gwas_chr2, gwas_chr3],
    output='all_chromosomes_manhattan.png'
)

# Or concatenate first
all_results = pd.concat([gwas_chr1, gwas_chr2, gwas_chr3])
qq_plot(all_results, output='all_chromosomes_qq.png')
```

### Custom Visualization (Alternative)

If you prefer to customize further:

```python
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

def custom_manhattan_plot(gwas_df, output='manhattan.png'):
    """Custom Manhattan plot with additional features"""
    gwas_df['-log10p'] = -np.log10(gwas_df['pval'])
    
    fig, ax = plt.subplots(figsize=(14, 6))
    colors = ['#1f77b4', '#ff7f0e']
    x_pos = 0
    
    for chrom in sorted(gwas_df['chr'].unique()):
        data = gwas_df[gwas_df['chr'] == chrom]
        ax.scatter(x_pos + np.arange(len(data)), data['-log10p'], 
                  c=colors[int(chrom) % 2], s=2, alpha=0.7)
        x_pos += len(data)
    
    # Significance lines
    ax.axhline(-np.log10(5e-8), color='red', linestyle='--', label='p=5e-8')
    ax.axhline(-np.log10(1e-5), color='blue', linestyle=':', label='p=1e-5')
    
    ax.set_xlabel('Chromosome', fontsize=12)
    ax.set_ylabel('-log₁₀(p-value)', fontsize=12)
    ax.set_title('EDGE GWAS Manhattan Plot', fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches='tight')
    plt.close()

def custom_qq_plot(gwas_df, output='qq_plot.png'):
    """Custom QQ plot with genomic inflation factor"""
    pvals = gwas_df['pval'].dropna()
    pvals = pvals[pvals > 0]
    
    n = len(pvals)
    observed = -np.log10(np.sort(pvals))
    expected = -np.log10(np.arange(1, n + 1) / (n + 1))
    
    # Genomic inflation factor
    chisq = stats.chi2.ppf(1 - pvals, df=1)
    lambda_gc = np.median(chisq) / stats.chi2.ppf(0.5, df=1)
    
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(expected, observed, s=10, alpha=0.5, color='blue')
    
    # Diagonal line
    max_val = max(expected.max(), observed.max())
    ax.plot([0, max_val], [0, max_val], 'r--', linewidth=2, label='Expected')
    
    # Lambda annotation
    ax.text(0.05, 0.95, f'λ = {lambda_gc:.3f}', 
            transform=ax.transAxes, fontsize=12, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    ax.set_xlabel('Expected -log₁₀(p)', fontsize=12)
    ax.set_ylabel('Observed -log₁₀(p)', fontsize=12)
    ax.set_title('QQ Plot', fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches='tight')
    plt.close()
    
    return lambda_gc

# Usage
custom_manhattan_plot(gwas_df)
lambda_gc = custom_qq_plot(gwas_df)
print(f"λ = {lambda_gc:.3f}")
```

## Complete Analysis Example

```python
"""Complete EDGE GWAS workflow with visualizations"""

from edge_gwas import (
    EDGEAnalysis,
    load_plink_data,
    prepare_phenotype_data,
    manhattan_plot,
    qq_plot,
    plot_alpha_distribution
)
import pandas as pd

# 1. Load data
genotype_data, variant_info = load_plink_data('data/genotypes')
phenotype_df = pd.read_csv('data/phenotypes.csv', index_col=0)

# 2. Prepare phenotype
phenotype_df = prepare_phenotype_data(
    phenotype_df,
    outcome='disease',
    covariates=['age', 'sex'] + [f'PC{i}' for i in range(1, 11)],
    remove_outliers=True
)

# 3. Initialize EDGE
edge = EDGEAnalysis(outcome_type='binary', verbose=True)

# 4. Split data
train_idx = phenotype_df.sample(frac=0.5, random_state=42).index
test_idx = phenotype_df.index.difference(train_idx)

# 5. Run analysis
alpha_df, gwas_df = edge.run_full_analysis(
    train_genotype=genotype_data.loc[train_idx],
    train_phenotype=phenotype_df.loc[train_idx],
    test_genotype=genotype_data.loc[test_idx],
    test_phenotype=phenotype_df.loc[test_idx],
    outcome='disease',
    covariates=['age', 'sex'] + [f'PC{i}' for i in range(1, 11)],
    output_prefix='results/edge_analysis'
)

# 6. Create visualizations
manhattan_plot(gwas_df, output='results/manhattan.png')
lambda_gc = qq_plot(gwas_df, output='results/qq_plot.png')
plot_alpha_distribution(alpha_df, output='results/alpha_dist.png')

# 7. Summary
print(f"\n=== Analysis Summary ===")
print(f"Total variants tested: {len(gwas_df)}")
print(f"Significant hits (p<5e-8): {(gwas_df['pval'] < 5e-8).sum()}")
print(f"Genomic inflation (λ): {lambda_gc:.3f}")
print(f"Mean alpha: {alpha_df['alpha_value'].mean():.3f}")
print(f"Median alpha: {alpha_df['alpha_value'].median():.3f}")
```

## Version History 

### Version 0.1.0 (Current - Under Public Testing)
**Release Date:** December 2025  
**Status:** ⚠️ Public Testing Phase

**Initial Release Features:**
- Two-stage EDGE analysis (calculate alpha + apply alpha)
- Support for binary and quantitative outcomes
- PLINK format data loading
- Built-in visualization functions (Manhattan plot, QQ plot, alpha distribution)
- Genomic inflation factor calculation

**Available Functions:**
- `EDGEAnalysis` class with methods:
  - `calculate_alpha()` - Calculate encoding parameters from training data
  - `apply_alpha()` - Apply encoding to test data
  - `run_full_analysis()` - Complete two-stage workflow
  - `get_skipped_snps()` - Retrieve SNPs with convergence issues
- `load_plink_data()` - Load genotype data from PLINK files
- `prepare_phenotype_data()` - QC and format phenotype data
- `manhattan_plot()` - Create Manhattan plots
- `qq_plot()` - Create QQ plots with lambda calculation
- `plot_alpha_distribution()` - Visualize alpha value distribution

**Known Issues:**
- Package under active development
- Documentation is being expanded
- Additional test coverage needed

### Version 0.0.0 (Plain Python Codes)
**Release Date:** April 2024  
**Status:** Closed, no further maintenance

**Legacy Features:**
- Basic EDGE analysis implementation
- Standalone Python scripts
- No package structure
- Manual dependency management

**Notes:**
- This version consisted of individual Python scripts
- Required manual setup and configuration
- No pip installation available
- Users should migrate to version 0.1.0 or later

**Repository:**
- Available at: https://github.com/nicenzhou/EDGE (original repository)
- For reference only, recommended for new projects in a Python notebook

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

[![GPL Logo](https://www.gnu.org/graphics/gplv3-88x31.png)](https://www.gnu.org/licenses/gpl-3.0)
