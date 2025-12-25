# EDGE GWAS Python Package Implementation Guide

## Current Version 0.1.0 (Under Testing)

⚠️ **Note:** This package is currently under active development and public testing.

> **Previous Version:**
> The original EDGE implementation (v0.0.0, plain Python scripts) is available at [https://github.com/nicenzhou/EDGE](https://github.com/nicenzhou/EDGE).
> The version 0.0.0 is **no longer maintained**; users are encouraged to migrate to v0.1.0.

## Overview

EDGE-GWAS (Elastic Data-Driven Encoding GWAS) identifies nonadditive SNP effects using flexible genetic encoding, rather than assuming additive inheritance.

**Key Features:**
- Two-stage analysis: calculate alpha on training data, apply to test data
- Detects under-recessive, recessive, additive, dominant, and over-dominant effects on a continuous scale
- Handles binary and quantitative outcomes
- Support for PLINK format data
- Built-in visualization functions

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

```python
""" Fast minimal EDGE GWAS analysis """
from edge_gwas import EDGEAnalysis
from edge_gwas.utils import load_plink_data, prepare_phenotype_data, stratified_train_test_split
from edge_gwas.visualize import manhattan_plot, qq_plot

# Load data
geno, info = load_plink_data('data.bed', 'data.bim', 'data.fam')
pheno = prepare_phenotype_data('pheno.txt', 'disease', ['age', 'sex', 'PC1', 'PC2'])

# Split and analyze
train_g, test_g, train_p, test_p = stratified_train_test_split(geno, pheno, 'disease')
edge = EDGEAnalysis(outcome_type='binary')
alpha_df, gwas_df = edge.run_full_analysis(train_g, train_p, test_g, test_p, 'disease', ['age', 'sex', 'PC1', 'PC2'])

# Visualize
manhattan_plot(gwas_df, 'manhattan.png')
qq_plot(gwas_df, 'qq.png')
```

## EDGE-GWAS Package Function Reference

### Complete Workflow with All Options

This section demonstrates a complete EDGE GWAS analysis workflow, showing different options at each step.

#### 1. Import and Initialize

```python
import pandas as pd
from edge_gwas import EDGEAnalysis
from edge_gwas.utils import (
    load_plink_data,
    prepare_phenotype_data,
    filter_variants_by_maf,
    filter_variants_by_missing,
    stratified_train_test_split,
    calculate_genomic_inflation,
    merge_alpha_with_gwas
)
from edge_gwas.visualize import (
    manhattan_plot,
    qq_plot,
    plot_alpha_distribution
)
from edge_gwas.io_handlers import (
    save_results,
    create_summary_report
)

# Initialize EDGE Analysis
edge = EDGEAnalysis(
    outcome_type='binary',    # 'binary' or 'continuous'
    n_jobs=8,                 # Number of CPU cores (-1 for all cores)
    max_iter=1000,            # Maximum iterations for model convergence
    verbose=True              # Print progress information
)
```

#### 2. Load and Prepare Data

```python
# Load genotype data from PLINK files
genotype_data, variant_info = load_plink_data(
    bed_file='path/to/data.bed',
    bim_file='path/to/data.bim',
    fam_file='path/to/data.fam',
    verbose=True
)

# Load phenotype data from file
phenotype_df = prepare_phenotype_data(
    phenotype_file='path/to/phenotypes.txt',
    outcome_col='disease',
    covariate_cols=['age', 'sex', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 
                    'PC6', 'PC7', 'PC8', 'PC9', 'PC10'],
    sample_id_col='IID',
    sep='\t',
    log_transform_outcome=False  # Set True for continuous outcomes if needed
)

# Define outcome and covariates
outcome = 'disease'
covariates = ['age', 'sex'] + [f'PC{i}' for i in range(1, 11)]

print(f"Loaded {len(genotype_data)} samples, {genotype_data.shape[1]} variants")
print(f"Phenotype data: {len(phenotype_df)} samples")
```

#### 3. Quality Control Filtering

```python
# Filter variants by MAF (remove rare variants)
genotype_data = filter_variants_by_maf(
    genotype_data,
    min_maf=0.01,
    verbose=True
)

# Filter variants by missing rate
genotype_data = filter_variants_by_missing(
    genotype_data,
    max_missing=0.1,
    verbose=True
)

print(f"After QC: {genotype_data.shape[1]} variants retained")
```

#### 4. Split Data into Training and Test Sets

**Option A: Using utility function with stratification**

```python
train_geno, test_geno, train_pheno, test_pheno = stratified_train_test_split(
    genotype_df=genotype_data,
    phenotype_df=phenotype_df,
    outcome_col=outcome,
    test_size=0.5,
    random_state=42,
    is_binary=True  # Stratify by outcome for balanced split
)
```

**Option B: Manual split**

```python
# Simple random split
train_samples = phenotype_df.sample(frac=0.5, random_state=42).index
test_samples = phenotype_df.index.difference(train_samples)

train_geno = genotype_data.loc[train_samples]
test_geno = genotype_data.loc[test_samples]
train_pheno = phenotype_df.loc[train_samples]
test_pheno = phenotype_df.loc[test_samples]
```

#### 5. Run EDGE Analysis

**Option A: One-Step Complete Analysis (Recommended for Quick Analysis)**

```python
# Run complete two-stage analysis in one call
alpha_df, gwas_df = edge.run_full_analysis(
    train_genotype=train_geno,
    train_phenotype=train_pheno,
    test_genotype=test_geno,
    test_phenotype=test_pheno,
    outcome=outcome,
    covariates=covariates,
    variant_info=variant_info,
    output_prefix='edge_results'  # Automatically saves results
)

print(f"Tested {len(gwas_df)} variants")
print(f"Significant (p < 5e-8): {(gwas_df['pval'] < 5e-8).sum()}")
```

**Option B: Two-Step Analysis (More Control and Flexibility)**

```python
# Step 1: Calculate alpha values on training data
print("Step 1: Calculating alpha values on training data...")
alpha_df = edge.calculate_alpha(
    genotype_data=train_geno,
    phenotype_df=train_pheno,
    outcome=outcome,
    covariates=covariates,
    variant_info=variant_info
)

# Save alpha values
alpha_df.to_csv('alpha_values.txt', sep='\t', index=False)
print(f"Alpha calculated for {len(alpha_df)} variants")
print(f"Mean alpha: {alpha_df['alpha_value'].mean():.3f}")

# Step 2: Apply alpha values on test data
print("\nStep 2: Applying alpha values on test data...")
gwas_df = edge.apply_alpha(
    genotype_data=test_geno,
    phenotype_df=test_pheno,
    outcome=outcome,
    covariates=covariates,
    alpha_values=alpha_df
)

# Save GWAS results
gwas_df.to_csv('gwas_results.txt', sep='\t', index=False)
print(f"Tested {len(gwas_df)} variants")
print(f"Significant (p < 5e-8): {(gwas_df['pval'] < 5e-8).sum()}")
```

**Option C: Use Pre-calculated Alpha Values**

```python
# Load previously calculated alpha values
alpha_df = pd.read_csv('alpha_values.txt', sep='\t')

# Or use load_alpha_values utility
from edge_gwas.io_handlers import load_alpha_values
alpha_df = load_alpha_values('alpha_values.txt')

# Apply to new test data
gwas_df = edge.apply_alpha(
    genotype_data=new_test_geno,
    phenotype_df=new_test_pheno,
    outcome=outcome,
    covariates=covariates,
    alpha_values=alpha_df
)
```

#### 6. Quality Control and Diagnostics

```python
# Calculate genomic inflation factor
lambda_gc = calculate_genomic_inflation(gwas_df['pval'])
print(f"Genomic inflation factor (λ): {lambda_gc:.3f}")

# Check for skipped SNPs
skipped_snps = edge.get_skipped_snps()
if len(skipped_snps) > 0:
    print(f"Warning: {len(skipped_snps)} SNPs skipped due to convergence issues")
    
# Merge alpha values with GWAS results for comprehensive output
merged_df = merge_alpha_with_gwas(gwas_df, alpha_df)
print(f"\nMerged results columns: {list(merged_df.columns)}")
```

#### 7. Create Visualizations

```python
# Manhattan plot
manhattan_plot(
    gwas_df,
    output='manhattan.png',
    title='EDGE GWAS Manhattan Plot',
    sig_threshold=5e-8,
    figsize=(14, 6)
)

# QQ plot with genomic inflation factor
lambda_gc = qq_plot(
    gwas_df,
    output='qq_plot.png',
    title='EDGE GWAS QQ Plot',
    figsize=(8, 8)
)

# Alpha distribution plot
plot_alpha_distribution(
    alpha_df,
    output='alpha_distribution.png',
    bins=50,
    figsize=(10, 6)
)

print(f"Genomic inflation factor (λ): {lambda_gc:.3f}")
```

#### 8. Save Results and Generate Report

```python
# Save results with standardized format
output_files = save_results(
    gwas_df=gwas_df,
    alpha_df=alpha_df,
    output_prefix='my_edge_gwas',
    save_alpha=True
)
print(f"Results saved to: {output_files}")

# Generate summary report
report = create_summary_report(
    gwas_df=gwas_df,
    alpha_df=alpha_df,
    significance_threshold=5e-8,
    output_file='summary_report.txt'
)
print(report)

# Save merged results
merged_df.to_csv('edge_gwas_complete_results.txt', sep='\t', index=False)
print(f"\nAnalysis complete! Files saved:")
print(f"  - Alpha values: {output_files['alpha']}")
print(f"  - GWAS results: {output_files['gwas']}")
print(f"  - Merged results: edge_gwas_complete_results.txt")
print(f"  - Manhattan plot: manhattan.png")
print(f"  - QQ plot: qq_plot.png")
print(f"  - Alpha distribution: alpha_distribution.png")
print(f"  - Summary report: summary_report.txt")
```

---

### Function Reference Summary

#### Core Module (`edge_gwas.core`)

**Class: `EDGEAnalysis`**

Main class for performing EDGE GWAS analysis.

```python
edge = EDGEAnalysis(
    outcome_type='binary',  # 'binary' or 'continuous'
    n_jobs=-1,              # Number of parallel jobs (-1 = all cores)
    max_iter=1000,          # Max iterations for model fitting
    verbose=True            # Print progress messages
)
```

**Methods:**

| Method | Purpose | Key Input | Output |
|--------|---------|-----------|--------|
| `calculate_alpha()` | Calculate encoding parameters from training data | Genotype, phenotype, covariates | Alpha DataFrame |
| `apply_alpha()` | Apply alpha values to test data for GWAS | Genotype, phenotype, alpha values | GWAS results DataFrame |
| `run_full_analysis()` | Complete two-stage workflow | Train/test genotype & phenotype | (alpha_df, gwas_df) |
| `get_skipped_snps()` | Get list of SNPs that failed convergence | None | List of variant IDs |

---

#### Utilities Module (`edge_gwas.utils`)

**Data Loading Functions:**

| Function | Purpose | Key Input | Output |
|----------|---------|-----------|--------|
| `load_plink_data()` | Load PLINK binary files | .bed, .bim, .fam paths | (genotype_df, variant_info_df) |
| `prepare_phenotype_data()` | Load and prepare phenotype data | File path, column names | Phenotype DataFrame |

**Data Processing Functions:**

| Function | Purpose | Key Input | Output |
|----------|---------|-----------|--------|
| `stratified_train_test_split()` | Split data into train/test sets | Genotype/phenotype DataFrames | (train_g, test_g, train_p, test_p) |
| `filter_variants_by_maf()` | Filter by minor allele frequency | Genotype DataFrame, MAF threshold | Filtered DataFrame |
| `filter_variants_by_missing()` | Filter by missingness rate | Genotype DataFrame, missing threshold | Filtered DataFrame |
| `merge_alpha_with_gwas()` | Merge GWAS and alpha results | GWAS and alpha DataFrames | Merged DataFrame |

**Statistical Functions:**

| Function | Purpose | Key Input | Output |
|----------|---------|-----------|--------|
| `calculate_genomic_inflation()` | Calculate lambda (λ) | P-values Series | Lambda (float) |
| `qq_plot_data()` | Prepare QQ plot data | P-values Series | (expected, observed) arrays |

---

#### Visualization Module (`edge_gwas.visualize`)

| Function | Purpose | Key Input | Output |
|----------|---------|-----------|--------|
| `manhattan_plot()` | Create Manhattan plot | GWAS DataFrame | Saves PNG file |
| `qq_plot()` | Create QQ plot | GWAS DataFrame | Saves PNG, returns lambda |
| `plot_alpha_distribution()` | Plot alpha histogram | Alpha DataFrame | Saves PNG file |

---

#### I/O Handlers Module (`edge_gwas.io_handlers`)

| Function | Purpose | Key Input | Output |
|----------|---------|-----------|--------|
| `save_results()` | Save GWAS and alpha results | DataFrames, output prefix | Dictionary of file paths |
| `load_alpha_values()` | Load pre-calculated alpha values | File path | Alpha DataFrame |
| `format_gwas_output()` | Format results for publication | GWAS DataFrame | Formatted DataFrame |
| `create_summary_report()` | Generate text summary | GWAS and alpha DataFrames | Summary report string |

---

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
- Automatic quality control filters (MAF, HWE, missingness)

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
- Available at: [https://github.com/nicenzhou/EDGE](https://github.com/nicenzhou/EDGE) (original repository)
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
