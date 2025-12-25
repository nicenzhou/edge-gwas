# EDGE GWAS Implementation Guide

## Overview

EDGE (Encoding-based Differential Genetic Effects) is a GWAS method that identifies optimal genetic encoding models for each variant, testing all possible inheritance patterns instead of assuming additive effects.

**Key Features:**
- Flexible genetic encoding (recessive, additive, dominant, etc.)
- Handles binary and quantitative outcomes
- Publication-ready visualizations

## Package Structure

<pre>
edge-gwas/
├── edge_gwas/
│   ├── __init__.py
│   ├── core.py
│   ├── utils.py
│   └── io_handlers.py
├── tests/
│   ├── __init__.py
│   └── test_core.py
├── examples/
│   ├── example_binary_outcome.py
│   └── example_continuous_outcome.py
├── setup.py
├── requirements.txt
├── README.md
├── LICENSE
└── .gitignore
</pre>

## Statistical Model

**Regression Model:**

$$E(Y | SNP_{Het}, SNP_{HA}, COV_i) = \beta_0 + \beta_{Het} \cdot SNP_{Het} + \beta_{HA} \cdot SNP_{HA} + \sum_{i} \beta_{cov_i} \cdot COV_i$$

**Encoding Parameter:**

$$\alpha = \frac{\beta_{Het}}{\beta_{HA}}$$

## Installation

```bash
pip install git+https://github.com/G2lab/UKBB_GWAS_pipeline.git
```

## Quick Start

### 1. Load Data

```python
import pandas as pd
import dxdata
from UKBB_GWAS_pipeline.edge_pipeline import EdgeGWAS

# Connect to UK Biobank
dataset_id = "your_dataset_id"
participant = dxdata.load_dataset(id=dataset_id)

# Define fields
fields = ["eid", "p31", "p21022"] + [f"p22009_a{i}" for i in range(1, 11)] + ["p41270"]

# Retrieve as pandas DataFrame
participant_df = participant.retrieve_fields(
    fields=fields, 
    engine=dxdata.connect()
).to_pandas()
```

### 2. Prepare Phenotype

```python
# Binary trait (e.g., Type 2 Diabetes)
participant_df['diabetes_t2'] = participant_df['p41270'].apply(
    lambda x: 1 if 'E11' in str(x) else 0
)

# Rename covariates
participant_df.rename(columns={
    'p31': 'sex',
    'p21022': 'age',
    **{f'p22009_a{i}': f'PC{i}' for i in range(1, 11)}
}, inplace=True)

# Define covariates
covariates = ['sex', 'age'] + [f'PC{i}' for i in range(1, 11)]

# QC: Remove missing data
participant_df = participant_df.dropna(subset=['diabetes_t2'] + covariates)
```

### 3. Run EDGE

```python
# Initialize
edge = EdgeGWAS(
    phenotype_name='diabetes_t2',
    phenotype_type='binary',
    covariates=covariates,
    output_dir='/path/to/output',
    maf_threshold=0.01
)

# Run genome-wide
for chrom in range(1, 23):
    print(f"Processing chromosome {chrom}...")
    edge.run_chromosome(
        pheno_df=participant_df,
        chromosome=chrom,
        dataset_id=dataset_id
    )
```

## Output Files

### Main Results: `edge_results_chr{N}.txt`

| Column | Description |
|--------|-------------|
| CHR | Chromosome |
| SNP | rsID |
| BP | Position |
| A1 | Effect allele |
| A2 | Reference allele |
| MAF | Minor allele frequency |
| BETA_HET | Heterozygous effect |
| BETA_HA | Homozygous alternate effect |
| ALPHA | Encoding parameter |
| ENCODING | Inheritance model |
| P_EDGE | EDGE p-value |
| P_ADD | Additive p-value |

## Visualization

### Manhattan Plot

```python
import matplotlib.pyplot as plt
import numpy as np

def manhattan_plot(results_files, output='manhattan.png'):
    # Combine results
    results = pd.concat([pd.read_csv(f, sep='\t') for f in results_files])
    results['-log10p'] = -np.log10(results['P_EDGE'])
    
    # Plot
    fig, ax = plt.subplots(figsize=(14, 6))
    colors = ['#1f77b4', '#ff7f0e']
    x_pos = 0
    
    for chrom in range(1, 23):
        data = results[results['CHR'] == chrom]
        ax.scatter(x_pos + np.arange(len(data)), data['-log10p'], 
                  c=colors[chrom % 2], s=2, alpha=0.7)
        x_pos += len(data)
    
    ax.axhline(-np.log10(5e-8), color='red', linestyle='--', label='p=5e-8')
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('-log₁₀(p-value)')
    ax.set_title('EDGE GWAS Manhattan Plot')
    plt.savefig(output, dpi=300, bbox_inches='tight')

# Usage
files = [f'/path/to/output/edge_results_chr{i}.txt' for i in range(1, 23)]
manhattan_plot(files)
```

### QQ Plot

```python
from scipy import stats

def qq_plot(results_files, output='qq_plot.png'):
    # Combine and filter
    results = pd.concat([pd.read_csv(f, sep='\t') for f in results_files])
    pvals = results['P_EDGE'].dropna()
    pvals = pvals[pvals > 0]
    
    # Calculate expected vs observed
    n = len(pvals)
    observed = -np.log10(np.sort(pvals))
    expected = -np.log10(np.arange(1, n + 1) / (n + 1))
    
    # Genomic inflation factor
    chisq = stats.chi2.ppf(1 - pvals, df=1)
    lambda_gc = np.median(chisq) / stats.chi2.ppf(0.5, df=1)
    
    # Plot
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(expected, observed, s=10, alpha=0.5)
    max_val = max(expected.max(), observed.max())
    ax.plot([0, max_val], [0, max_val], 'r--', linewidth=2)
    ax.text(0.05, 0.95, f'λ = {lambda_gc:.3f}', transform=ax.transAxes,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax.set_xlabel('Expected -log₁₀(p)')
    ax.set_ylabel('Observed -log₁₀(p)')
    ax.set_title('QQ Plot')
    plt.savefig(output, dpi=300, bbox_inches='tight')
    print(f"Genomic inflation factor: {lambda_gc:.3f}")

# Usage
files = [f'/path/to/output/edge_results_chr{i}.txt' for i in range(1, 23)]
qq_plot(files)
```

## Best Practices

1. **QC Filters:** MAF ≥ 0.01, genotype missingness < 5%, HWE p > 1e-6
2. **Covariates:** Always include sex, age, and 10 PCs
3. **Sample Size:** Cases ≥ 1000, Controls ≥ 5000 for binary traits
4. **Significance:** Genome-wide threshold p < 5×10⁻⁸
5. **λ Interpretation:** λ = 1.0-1.05 acceptable; >1.1 suggests stratification

## Applications

EDGE has been successfully applied in large-scale genomic studies:

- **UK Biobank (UKB)**: Genome-wide analyses across diverse phenotypes
- **Million Veteran Program (MVP)**: Large-scale veteran health genomics research

### Implementation Code

For the complete EDGE methodology and detailed implementation, see:

**Original EDGE Method:**
- Repository: [https://github.com/nicenzhou/EDGE](https://github.com/nicenzhou/EDGE)
- Publication: [EDGE GWAS_Preprint](https://doi.org/10.1101/2023.06.01.23290857)

## Citation

If you use EDGE in your research, please cite: 
Zhou, J., Rico, A. L. G., Guare, L., Million Veteran Program, Chang, K. M., Tsao, P. S., Assimes, T. L., Verma, S. S., & Hall, M. A. (2023). Flexibly encoded genome-wide association study identifies novel nonadditive genetic risk variants for cardiometabolic traits. *medRxiv*, 2023.06.01.23290857. https://doi.org/10.1101/2023.06.01.23290857

**BibTeX:**
```bibtex
@article{zhou2023edgegwas,
  title={Flexibly encoded genome-wide association study identifies novel nonadditive genetic risk variants for cardiometabolic traits},
  author={Zhou, Jiayan and Rico, Andre Luis Garao and Guare, Lindsay and Million Veteran Program and Chang, Kyong-Mi and Tsao, Philip S and Assimes, Themistocles L and Verma, Shefali Setia and Hall, Molly Ann},
  journal={medRxiv},
  pages={2023--06},
  year={2023},
  publisher={Cold Spring Harbor Laboratory Press},
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

