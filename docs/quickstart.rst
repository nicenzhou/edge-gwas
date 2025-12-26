.. _quickstart:

Quick Start Guide
=================

This guide demonstrates a complete EDGE GWAS analysis in minutes using v0.1.1 features.

Minimal Example
----------------------------

Fast minimal EDGE GWAS analysis:

.. code-block:: python

   from edge_gwas import EDGEAnalysis
   from edge_gwas.utils import (
       load_plink_data,
       prepare_phenotype_data,
       stratified_train_test_split
   )
   from edge_gwas.visualize import manhattan_plot, qq_plot
   
   # Load data
   geno, info = load_plink_data('data.bed', 'data.bim', 'data.fam')
   pheno = prepare_phenotype_data(
       'pheno.txt', 'disease', ['age', 'sex', 'PC1', 'PC2']
   )
   
   # Split data
   train_g, test_g, train_p, test_p = stratified_train_test_split(
       geno, pheno, 'disease', test_size=0.5
   )
   
   # Run analysis
   edge = EDGEAnalysis(outcome_type='binary', n_jobs=-1)
   alpha_df, gwas_df = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       'disease', ['age', 'sex', 'PC1', 'PC2']
   )
   
   # Visualize
   manhattan_plot(gwas_df, 'manhattan.png')
   lambda_gc = qq_plot(gwas_df, 'qq.png')
   print(f"Lambda GC: {lambda_gc:.3f}")

Complete Example with v0.1.1 Features
---------------------------------------------------

Full workflow with population structure control:

.. code-block:: python

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
   print("Loading data...")
   geno, info = load_plink_data('data.bed', 'data.bim', 'data.fam')
   pheno = prepare_phenotype_data('pheno.txt', 'disease', ['age', 'sex'])
   
   # 2. QC filtering
   print("Applying QC...")
   geno = filter_variants_by_maf(geno, min_maf=0.01)
   
   # 3. Calculate PCA (NEW in v0.1.1)
   print("Calculating PCA...")
   pca_df = calculate_pca_plink('data', n_pcs=10)
   pheno = attach_pcs_to_phenotype(pheno, pca_df, n_pcs=10)
   
   # 4. Split data
   print("Splitting data...")
   train_g, test_g, train_p, test_p = stratified_train_test_split(
       geno, pheno, 'disease', test_size=0.5
   )
   
   # 5. Run EDGE
   print("Running EDGE analysis...")
   edge = EDGEAnalysis(outcome_type='binary', n_jobs=-1, verbose=True)
   
   covariates = ['age', 'sex'] + get_pc_covariate_list(10)
   
   alpha_df, gwas_df = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       outcome='disease',
       covariates=covariates,
       output_prefix='results/edge'
   )
   
   # 6. Visualize
   print("Creating plots...")
   manhattan_plot(gwas_df, 'results/manhattan.png')
   lambda_gc = qq_plot(gwas_df, 'results/qq.png')
   
   print(f"\nResults:")
   print(f"  Variants tested: {len(gwas_df):,}")
   print(f"  Significant (p < 5e-8): {(gwas_df['pval'] < 5e-8).sum()}")
   print(f"  Lambda GC: {lambda_gc:.3f}")

Step-by-Step Workflow
----------------------

1. Install and Setup
~~~~~~~~~~~~~~~~~~~~

First, ensure edge-gwas and tools are installed:

.. code-block:: bash

   # Install edge-gwas
   pip install git+https://github.com/nicenzhou/edge-gwas.git
   
   # Install external tools (optional, NEW in v0.1.1)
   edge-gwas-install-tools
   
   # Verify installation
   edge-gwas-check-tools

2. Import Packages
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import pandas as pd
   from edge_gwas import EDGEAnalysis
   from edge_gwas.utils import (
       load_plink_data,
       prepare_phenotype_data,
       filter_variants_by_maf,
       filter_variants_by_hwe,
       calculate_pca_plink,
       attach_pcs_to_phenotype,
       get_pc_covariate_list,
       stratified_train_test_split
   )
   from edge_gwas.visualize import manhattan_plot, qq_plot, plot_alpha_distribution

3. Load Data
~~~~~~~~~~~~

**Load genotype data (multiple formats supported):**

.. code-block:: python

   # Option 1: PLINK binary format (most common)
   genotype_df, variant_info = load_plink_data(
       bed_file='path/to/data.bed',
       bim_file='path/to/data.bim',
       fam_file='path/to/data.fam'
   )
   
   # Option 2: PLINK2 format (NEW in v0.1.1)
   from edge_gwas.utils import load_pgen_data
   genotype_df, variant_info = load_pgen_data(
       pgen_file='data.pgen',
       pvar_file='data.pvar',
       psam_file='data.psam'
   )
   
   # Option 3: BGEN format (NEW in v0.1.1)
   from edge_gwas.utils import load_bgen_data
   genotype_df, variant_info = load_bgen_data(
       bgen_file='data.bgen',
       sample_file='data.sample'
   )
   
   # Option 4: VCF format (NEW in v0.1.1)
   from edge_gwas.utils import load_vcf_data
   genotype_df, variant_info = load_vcf_data(
       vcf_file='data.vcf.gz',
       dosage=True
   )

**Load phenotype data:**

.. code-block:: python

   phenotype_df = prepare_phenotype_data(
       phenotype_file='path/to/phenotypes.txt',
       outcome_col='disease',
       covariate_cols=['age', 'sex'],
       sample_id_col='IID',
       sep='\t'
   )
   
   print(f"Loaded {len(genotype_df)} samples, {genotype_df.shape[1]} variants")
   print(f"Loaded {len(phenotype_df)} samples with phenotypes")

4. Quality Control
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Filter by MAF
   genotype_df = filter_variants_by_maf(
       genotype_df, min_maf=0.01, verbose=True
   )
   
   # Filter by HWE (NEW in v0.1.1)
   from edge_gwas.utils import filter_variants_by_hwe
   genotype_df = filter_variants_by_hwe(
       genotype_df, hwe_threshold=1e-6, verbose=True
   )
   
   # Check case/control balance (NEW in v0.1.1)
   from edge_gwas.utils import check_case_control_balance
   balance = check_case_control_balance(phenotype_df, 'disease', verbose=True)

5. Population Structure Control (NEW in v0.1.1)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Option A: Calculate PCA (for unrelated samples):**

.. code-block:: python

   # Calculate PCA using PLINK2
   pca_df = calculate_pca_plink(
       plink_prefix='path/to/data',
       n_pcs=10,
       maf_threshold=0.01,
       ld_r2=0.2,
       verbose=True
   )
   
   # Attach PCs to phenotype
   phenotype_df = attach_pcs_to_phenotype(
       phenotype_df, pca_df, n_pcs=10, verbose=True
   )

**Option B: Calculate PC-AiR (for related samples):**

.. code-block:: python

   from edge_gwas.utils import calculate_pca_pcair
   
   # PC-AiR accounts for relatedness
   pca_df = calculate_pca_pcair(
       plink_prefix='path/to/data',
       n_pcs=10,
       kin_threshold=0.0884,
       verbose=True
   )
   
   phenotype_df = attach_pcs_to_phenotype(
       phenotype_df, pca_df, n_pcs=10
   )

**Option C: Use GRM for mixed model (for related samples):**

.. code-block:: python

   from edge_gwas.utils import calculate_grm_gcta, load_grm_gcta
   
   # Calculate GRM
   grm_prefix = calculate_grm_gcta(
       plink_prefix='path/to/data',
       output_prefix='output/grm',
       maf_threshold=0.01,
       verbose=True
   )
   
   # Load GRM
   grm_matrix, grm_ids = load_grm_gcta('output/grm', verbose=True)

6. Initialize EDGE Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # For binary outcomes (case-control)
   edge = EDGEAnalysis(
       outcome_type='binary',
       n_jobs=-1,  # Use all CPU cores
       verbose=True
   )
   
   # For continuous outcomes with transformation (NEW in v0.1.1)
   edge = EDGEAnalysis(
       outcome_type='continuous',
       outcome_transform='rank_inverse_normal',  # Transform outcome
       n_jobs=-1,
       verbose=True
   )

7. Split Data
~~~~~~~~~~~~~

.. code-block:: python

   # Split with stratification
   train_g, test_g, train_p, test_p = stratified_train_test_split(
       genotype_df=genotype_df,
       phenotype_df=phenotype_df,
       outcome_col='disease',
       test_size=0.5,
       random_state=42,
       is_binary=True  # Set False for continuous outcomes
   )
   
   print(f"Training: {len(train_g)} samples")
   print(f"Testing: {len(test_g)} samples")

8. Run EDGE Analysis
~~~~~~~~~~~~~~~~~~~~

**Basic analysis:**

.. code-block:: python

   # Build covariate list
   covariates = ['age', 'sex'] + get_pc_covariate_list(10)
   
   # Run complete analysis
   alpha_df, gwas_df = edge.run_full_analysis(
       train_genotype=train_g,
       train_phenotype=train_p,
       test_genotype=test_g,
       test_phenotype=test_p,
       outcome='disease',
       covariates=covariates,
       variant_info=variant_info,
       output_prefix='results/edge'
   )

**With GRM (NEW in v0.1.1):**

.. code-block:: python

   # Run with GRM for population structure control
   alpha_df, gwas_df = edge.run_full_analysis(
       train_genotype=train_g,
       train_phenotype=train_p,
       test_genotype=test_g,
       test_phenotype=test_p,
       outcome='disease',
       covariates=covariates,
       grm_matrix=grm_matrix,
       grm_sample_ids=grm_ids,
       output_prefix='results/edge_with_grm'
   )

9. Examine Results
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Summary statistics
   print(f"\nResults Summary:")
   print(f"  Total variants tested: {len(gwas_df):,}")
   print(f"  Significant (p < 5e-8): {(gwas_df['pval'] < 5e-8).sum()}")
   print(f"  Suggestive (p < 1e-5): {((gwas_df['pval'] >= 5e-8) & (gwas_df['pval'] < 1e-5)).sum()}")
   
   # Alpha statistics
   print(f"\nAlpha Value Statistics:")
   print(f"  Mean: {alpha_df['alpha_value'].mean():.3f}")
   print(f"  Median: {alpha_df['alpha_value'].median():.3f}")
   print(f"  Std: {alpha_df['alpha_value'].std():.3f}")
   
   # Top hits
   top_hits = gwas_df.nsmallest(10, 'pval')
   print(f"\nTop 10 hits:")
   print(top_hits[['variant_id', 'chr', 'pos', 'pval', 'coef', 'alpha_value']])

10. Visualize Results
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Manhattan plot
   manhattan_plot(
       gwas_df,
       output='results/manhattan.png',
       title='EDGE GWAS Manhattan Plot',
       sig_threshold=5e-8,
       suggestive_threshold=1e-5
   )
   
   # QQ plot with genomic inflation
   from edge_gwas.utils import calculate_genomic_inflation
   lambda_gc = qq_plot(
       gwas_df,
       output='results/qq_plot.png',
       title='QQ Plot'
   )
   
   print(f"Genomic inflation factor (λ): {lambda_gc:.3f}")
   
   # Alpha distribution
   plot_alpha_distribution(
       alpha_df,
       output='results/alpha_distribution.png',
       title='Alpha Value Distribution'
   )

Understanding the Output
------------------------

Alpha Values DataFrame
~~~~~~~~~~~~~~~~~~~~~~

The alpha DataFrame contains encoding parameters calculated from training data:

.. code-block:: python

   print(alpha_df.head())
   
   #    variant_id  alpha_value    eaf  coef_het  coef_hom  ...
   # 0       rs123         0.45  0.235      0.23      0.51  ...
   # 1       rs456         0.52  0.412      0.31      0.60  ...

**Key columns:**

* ``variant_id``: SNP identifier
* ``alpha_value``: Encoding parameter (β_het/β_hom)
  
  * α ≈ 0: Recessive model
  * α ≈ 0.5: Additive model
  * α ≈ 1: Dominant model

* ``eaf``: Effect allele frequency
* ``coef_het``: Heterozygous coefficient (β_het)
* ``coef_hom``: Homozygous coefficient (β_hom)
* ``pval_het``, ``pval_hom``: P-values for het/hom effects

GWAS Results DataFrame
~~~~~~~~~~~~~~~~~~~~~~

The GWAS DataFrame contains association test results:

.. code-block:: python

   print(gwas_df.head())
   
   #    variant_id  chr      pos       pval    coef  std_err  alpha_value
   # 0       rs123    1  1234567  1.23e-08  0.0451   0.0089         0.45
   # 1       rs456    1  2345678  3.45e-06  0.0312   0.0067         0.52

**Key columns:**

* ``variant_id``: SNP identifier
* ``chr``: Chromosome number
* ``pos``: Genomic position
* ``pval``: Association p-value
* ``coef``: Effect size coefficient
* ``std_err``: Standard error of coefficient
* ``stat``: Test statistic
* ``alpha_value``: Applied encoding parameter
* ``n_samples``: Number of samples in analysis

Interpreting Alpha Values
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Categorize inheritance patterns
   def categorize_alpha(alpha):
       if alpha < 0.3:
           return 'Recessive'
       elif alpha > 0.7:
           return 'Dominant'
       else:
           return 'Additive'
   
   alpha_df['inheritance'] = alpha_df['alpha_value'].apply(categorize_alpha)
   
   print("\nInheritance pattern distribution:")
   print(alpha_df['inheritance'].value_counts())
   
   # Recessive     1234
   # Additive      5678
   # Dominant       890

Quick Recipes
-------------

Recipe 1: Continuous Trait with Transformation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # For skewed continuous traits (e.g., gene expression)
   edge = EDGEAnalysis(
       outcome_type='continuous',
       outcome_transform='rank_inverse_normal'
   )
   
   alpha_df, gwas_df = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       outcome='expression_level',
       covariates=['age', 'sex', 'PC1', 'PC2', 'PC3']
   )

Recipe 2: Analysis with Related Samples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from edge_gwas.utils import (
       calculate_grm_gcta,
       load_grm_gcta,
       calculate_pca_pcair,
       attach_pcs_to_phenotype
   )
   
   # Calculate GRM
   grm_prefix = calculate_grm_gcta('data')
   grm_matrix, grm_ids = load_grm_gcta(grm_prefix)
   
   # Calculate PC-AiR (accounts for relatedness)
   pca_df = calculate_pca_pcair('data', n_pcs=10, kinship_matrix=grm_prefix)
   pheno = attach_pcs_to_phenotype(pheno, pca_df, n_pcs=10)
   
   # Run EDGE with GRM
   edge = EDGEAnalysis(outcome_type='binary')
   alpha_df, gwas_df = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       outcome='disease',
       covariates=['age', 'sex'] + get_pc_covariate_list(10),
       grm_matrix=grm_matrix,
       grm_sample_ids=grm_ids
   )

Recipe 3: Compare EDGE with Additive Model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from edge_gwas.utils import additive_gwas
   import numpy as np
   
   # Run EDGE
   alpha_df, edge_results = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       'disease', ['age', 'sex']
   )
   
   # Run standard additive GWAS
   additive_results = additive_gwas(
       test_g, test_p, 'disease', ['age', 'sex'], 'binary'
   )
   
   # Compare
   comparison = edge_results[['variant_id', 'pval']].merge(
       additive_results[['variant_id', 'pval']],
       on='variant_id', suffixes=('_edge', '_additive')
   )
   
   # Find where EDGE is better
   comparison['edge_advantage'] = (
       -np.log10(comparison['pval_edge']) - 
       -np.log10(comparison['pval_additive'])
   )
   
   better = comparison[comparison['edge_advantage'] > 2].sort_values(
       'edge_advantage', ascending=False
   )
   print(f"EDGE shows >100-fold improvement for {len(better)} variants")

Recipe 4: Cross-Validation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from edge_gwas.utils import cross_validated_edge_analysis
   
   # 5-fold cross-validation
   avg_alpha, meta_gwas, all_alpha, all_gwas = cross_validated_edge_analysis(
       genotype_df=geno,
       phenotype_df=pheno,
       outcome='disease',
       covariates=['age', 'sex'] + get_pc_covariate_list(10),
       outcome_type='binary',
       n_folds=5,
       n_jobs=8
   )
   
   # Check alpha stability
   unstable = avg_alpha[avg_alpha['alpha_std'] > 0.3]
   print(f"Unstable alpha variants: {len(unstable)}")

Recipe 5: Multi-Format Data Loading
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from edge_gwas.utils import load_pgen_data, load_bgen_data, load_vcf_data
   
   # Load different formats - all return same structure
   
   # PLINK2 format
   geno1, info1 = load_pgen_data('data.pgen', 'data.pvar', 'data.psam')
   
   # BGEN format (UK Biobank)
   geno2, info2 = load_bgen_data('data.bgen', 'data.sample')
   
   # VCF format
   geno3, info3 = load_vcf_data('data.vcf.gz', dosage=True)
   
   # All can be used identically in EDGE analysis
   edge = EDGEAnalysis(outcome_type='binary')
   alpha_df, gwas_df = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       'disease', ['age', 'sex']
   )

Common Workflow Patterns
------------------------

Pattern 1: Basic Case-Control Study
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # 1. Load and QC
   geno, info = load_plink_data('data.bed', 'data.bim', 'data.fam')
   geno = filter_variants_by_maf(geno, min_maf=0.01)
   
   # 2. Load phenotype
   pheno = prepare_phenotype_data('pheno.txt', 'disease', ['age', 'sex'])
   
   # 3. Add PCs
   pca_df = calculate_pca_plink('data', n_pcs=10)
   pheno = attach_pcs_to_phenotype(pheno, pca_df, n_pcs=10)
   
   # 4. Split and analyze
   train_g, test_g, train_p, test_p = stratified_train_test_split(
       geno, pheno, 'disease'
   )
   
   edge = EDGEAnalysis(outcome_type='binary', n_jobs=-1)
   alpha_df, gwas_df = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       'disease', ['age', 'sex'] + get_pc_covariate_list(10)
   )
   
   # 5. Visualize
   manhattan_plot(gwas_df, 'manhattan.png')
   qq_plot(gwas_df, 'qq.png')

Pattern 2: Quantitative Trait Study
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Use outcome transformation for non-normal traits
   edge = EDGEAnalysis(
       outcome_type='continuous',
       outcome_transform='rank_inverse_normal',
       n_jobs=-1
   )
   
   alpha_df, gwas_df = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       outcome='BMI',
       covariates=['age', 'sex'] + get_pc_covariate_list(10)
   )

Pattern 3: Biobank-Scale Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # For large cohorts with relatedness
   
   # 1. Calculate GRM (once)
   grm_prefix = calculate_grm_gcta('data', method='grm-sparse')
   grm_matrix, grm_ids = load_grm_gcta(grm_prefix)
   
   # 2. Calculate PC-AiR (accounts for relatedness)
   pca_df = calculate_pca_pcair(
       'data', n_pcs=20, kinship_matrix=grm_prefix, approx=True
   )
   pheno = attach_pcs_to_phenotype(pheno, pca_df, n_pcs=20)
   
   # 3. Run EDGE with GRM
   edge = EDGEAnalysis(outcome_type='binary', n_jobs=-1)
   alpha_df, gwas_df = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       'disease',
       ['age', 'sex', 'batch'] + get_pc_covariate_list(20),
       grm_matrix=grm_matrix,
       grm_sample_ids=grm_ids
   )

Troubleshooting Quick Fixes
----------------------------

**Problem: "No common samples found"**

.. code-block:: python

   # Fix: Ensure IDs are strings
   geno.index = geno.index.astype(str)
   pheno.index = pheno.index.astype(str)
   pca_df.index = pca_df.index.astype(str)

**Problem: High genomic inflation (λ > 1.05)**

.. code-block:: python

   # Fix: Add more PCs or use GRM
   
   # Option 1: Add more PCs
   pca_df = calculate_pca_plink('data', n_pcs=20)  # Increase from 10
   
   # Option 2: Use GRM
   grm_matrix, grm_ids = load_grm_gcta('grm_prefix')
   alpha_df, gwas_df = edge.run_full_analysis(
       ..., grm_matrix=grm_matrix, grm_sample_ids=grm_ids
   )

**Problem: Memory errors**

.. code-block:: python

   # Fix: Reduce parallel jobs or process by chromosome
   edge = EDGEAnalysis(outcome_type='binary', n_jobs=4)  # Reduce from -1

**Problem: PLINK2/GCTA not found**

.. code-block:: bash

   # Fix: Run installer
   edge-gwas-install-tools
   
   # Or add to PATH
   export PATH="$HOME/.local/bin:$PATH"

**Problem: Slow performance**

.. code-block:: python

   # Fix: Use all cores and approximate PCA for large data
   edge = EDGEAnalysis(outcome_type='binary', n_jobs=-1)
   
   pca_df = calculate_pca_plink(
       'data', n_pcs=10, approx=True, approx_samples=5000
   )

What's New in v0.1.1?
---------------------

If you're upgrading from v0.1.0, here's what's new:

**1. No more Koalas! Pure pandas:**

.. code-block:: python

   # OLD (v0.1.0)
   import databricks.koalas as ks
   df = data.to_koalas()
   
   # NEW (v0.1.1)
   import pandas as pd
   df = data  # Already pandas!

**2. Population structure control:**

.. code-block:: python

   # Calculate PCA
   pca_df = calculate_pca_plink('data', n_pcs=10)
   
   # Calculate GRM
   grm_matrix, grm_ids = load_grm_gcta('grm_prefix')
   
   # Use in analysis
   alpha_df, gwas_df = edge.run_full_analysis(
       ..., grm_matrix=grm_matrix, grm_sample_ids=grm_ids
   )

**3. Outcome transformations:**

.. code-block:: python

   edge = EDGEAnalysis(
       outcome_type='continuous',
       outcome_transform='rank_inverse_normal'  # NEW!
   )

**4. Multiple file formats:**

.. code-block:: python

   from edge_gwas.utils import load_pgen_data, load_bgen_data, load_vcf_data
   
   geno, info = load_pgen_data('data.pgen', 'data.pvar', 'data.psam')

**5. Enhanced QC:**

.. code-block:: python

   from edge_gwas.utils import (
       filter_variants_by_hwe,
       filter_samples_by_call_rate,
       check_case_control_balance
   )

**6. Comparison tools:**

.. code-block:: python

   from edge_gwas.utils import additive_gwas, cross_validated_edge_analysis

---

*Last updated: 2025-12-25 for edge-gwas v0.1.1*

*For questions or issues, visit:* https://github.com/nicenzhou/edge-gwas/issues
