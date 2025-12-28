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

   from edge_gwas import EDGEAnalysis, download_test_files
   from edge_gwas.utils import (
       load_plink_data, prepare_phenotype_data,
       filter_genotype_data, validate_and_align_data,
       calculate_grm_gcta, load_grm_gcta,
       calculate_pca_plink, attach_pcs_to_phenotype,
       stratified_train_test_split, get_pc_covariate_list
   )
   from edge_gwas.visualize import manhattan_plot, qq_plot, plot_alpha_distribution
   
   # 1. Download test data
   download_test_files()
   
   # 2. Load data
   geno, info = load_plink_data('tests/test.bed', 'tests/test.bim', 'tests/test.fam')
   pheno = prepare_phenotype_data('tests/test.phen', outcome_col='disease', 
                                   covariate_cols=['age'], sample_id_col='IID')
   
   # 3. Quality control
   geno_qc, pheno_qc = filter_genotype_data(
       geno, pheno,
       min_maf=0.01,
       max_missing_per_variant=0.05,
       min_call_rate_per_sample=0.95
   )
   
   # 4. Population structure control
   grm_prefix = calculate_grm_gcta('tests/test', output_prefix='grm')
   grm_matrix, grm_ids = load_grm_gcta('grm')
   pca_df = calculate_pca_plink('tests/test', n_pcs=10)
   pheno_qc = attach_pcs_to_phenotype(pheno_qc, pca_df, n_pcs=10, drop_na=True)
   
   # 5. Train/test split
   train_g, test_g, train_p, test_p = stratified_train_test_split(
       geno_qc, pheno_qc, outcome_col='disease', test_size=0.5, random_state=42
   )
   
   # 6. Run EDGE analysis
   edge = EDGEAnalysis(outcome_type='binary')
   covariates = ['age'] + get_pc_covariate_list(10)
   
   alpha_df, gwas_df = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       outcome='disease',
       covariates=covariates,
       variant_info=info,
       grm_matrix=grm_matrix,
       grm_sample_ids=grm_ids,
       output_prefix='results/edge'
   )
   
   # 7. Visualize
   manhattan_plot(gwas_df, output='results/manhattan.png')
   lambda_gc = qq_plot(gwas_df, output='results/qq.png')
   plot_alpha_distribution(alpha_df, output='results/alpha_dist.png')
   
   # 8. Summary
   print(f"Variants tested: {len(gwas_df)}")
   print(f"Significant (p<5e-8): {(gwas_df['pval'] < 5e-8).sum()}")
   print(f"Genomic inflation: λ = {lambda_gc:.3f}")

Step-by-Step Workflow
----------------------

1. Install and Setup
~~~~~~~~~~~~~~~~~~~~

First, ensure edge-gwas and tools are installed:

.. code-block:: bash

   pip install git+https://github.com/nicenzhou/edge-gwas.git
   
   edge-gwas-install-tools
   
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

Load genotype data (multiple formats supported):

.. code-block:: python

   genotype_df, variant_info = load_plink_data(
       bed_file='path/to/data.bed',
       bim_file='path/to/data.bim',
       fam_file='path/to/data.fam'
   )
   
   from edge_gwas.utils import load_pgen_data
   genotype_df, variant_info = load_pgen_data(
       pgen_file='data.pgen',
       pvar_file='data.pvar',
       psam_file='data.psam'
   )
   
   from edge_gwas.utils import load_bgen_data
   genotype_df, variant_info = load_bgen_data(
       bgen_file='data.bgen',
       sample_file='data.sample'
   )
   
   from edge_gwas.utils import load_vcf_data
   genotype_df, variant_info = load_vcf_data(
       vcf_file='data.vcf.gz',
       dosage=True
   )

Load phenotype data:

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

   genotype_df = filter_variants_by_maf(
       genotype_df, min_maf=0.01, verbose=True
   )
   
   from edge_gwas.utils import filter_variants_by_hwe
   genotype_df = filter_variants_by_hwe(
       genotype_df, hwe_threshold=1e-6, verbose=True
   )
   
   from edge_gwas.utils import check_case_control_balance
   balance = check_case_control_balance(phenotype_df, 'disease', verbose=True)

5. Population Structure Control (NEW in v0.1.1)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Option A: Calculate PCA (for unrelated samples):

.. code-block:: python

   pca_df = calculate_pca_plink(
       plink_prefix='path/to/data',
       n_pcs=10,
       maf_threshold=0.01,
       ld_r2=0.2,
       verbose=True
   )
   
   phenotype_df = attach_pcs_to_phenotype(
       phenotype_df, pca_df, n_pcs=10, verbose=True
   )

Option B: Calculate PC-AiR (for related samples):

.. code-block:: python

   from edge_gwas.utils import calculate_pca_pcair
   
   pca_df = calculate_pca_pcair(
       plink_prefix='path/to/data',
       n_pcs=10,
       kin_threshold=0.0884,
       verbose=True
   )
   
   phenotype_df = attach_pcs_to_phenotype(
       phenotype_df, pca_df, n_pcs=10
   )

Option C: Use GRM for mixed model (for related samples):

.. code-block:: python

   from edge_gwas.utils import calculate_grm_gcta, load_grm_gcta
   
   grm_prefix = calculate_grm_gcta(
       plink_prefix='path/to/data',
       output_prefix='output/grm',
       maf_threshold=0.01,
       verbose=True
   )
   
   grm_matrix, grm_ids = load_grm_gcta('output/grm', verbose=True)


6. Initialize EDGE Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: python

   edge = EDGEAnalysis(
       outcome_type='binary',
       n_jobs=-1,
       verbose=True
   )
   
   edge = EDGEAnalysis(
       outcome_type='continuous',
       outcome_transform='rank_inverse_normal',
       n_jobs=-1,
       verbose=True
   )

7. Split Data
~~~~~~~~~~~~~

.. code-block:: python

   train_g, test_g, train_p, test_p = stratified_train_test_split(
       genotype_df=genotype_df,
       phenotype_df=phenotype_df,
       outcome_col='disease',
       test_size=0.5,
       random_state=42,
       is_binary=True
   )
   
   print(f"Training: {len(train_g)} samples")
   print(f"Testing: {len(test_g)} samples")

8. Run EDGE Analysis
~~~~~~~~~~~~~~~~~~~~

Basic analysis:

.. code-block:: python

   covariates = ['age', 'sex'] + get_pc_covariate_list(10)
   
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

With GRM (NEW in v0.1.1):

.. code-block:: python

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

   print(f"\nResults Summary:")
   print(f"  Total variants tested: {len(gwas_df):,}")
   print(f"  Significant (p < 5e-8): {(gwas_df['pval'] < 5e-8).sum()}")
   print(f"  Suggestive (p < 1e-5): {((gwas_df['pval'] >= 5e-8) & (gwas_df['pval'] < 1e-5)).sum()}")
   
   print(f"\nAlpha Value Statistics:")
   print(f"  Mean: {alpha_df['alpha_value'].mean():.3f}")
   print(f"  Median: {alpha_df['alpha_value'].median():.3f}")
   print(f"  Std: {alpha_df['alpha_value'].std():.3f}")
   
   top_hits = gwas_df.nsmallest(10, 'pval')
   print(f"\nTop 10 hits:")
   print(top_hits[['variant_id', 'chr', 'pos', 'pval', 'coef', 'alpha_value']])

10. Visualize Results
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   manhattan_plot(
       gwas_df,
       output='results/manhattan.png',
       title='EDGE GWAS Manhattan Plot',
       sig_threshold=5e-8,
       suggestive_threshold=1e-5
   )
   
   from edge_gwas.utils import calculate_genomic_inflation
   lambda_gc = qq_plot(
       gwas_df,
       output='results/qq_plot.png',
       title='QQ Plot'
   )
   
   print(f"Genomic inflation factor (λ): {lambda_gc:.3f}")
   
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

Key columns:

* variant_id: SNP identifier
* alpha_value: Encoding parameter (β_het/β_hom)
  
  * α ≈ 0: Recessive model
  * α ≈ 0.5: Additive model
  * α ≈ 1: Dominant model

* eaf: Effect allele frequency
* coef_het: Heterozygous coefficient (β_het)
* coef_hom: Homozygous coefficient (β_hom)
* pval_het, pval_hom: P-values for het/hom effects

GWAS Results DataFrame
~~~~~~~~~~~~~~~~~~~~~~

The GWAS DataFrame contains association test results:

.. code-block:: python

   print(gwas_df.head())

Key columns:

* variant_id: SNP identifier
* chr: Chromosome number
* pos: Genomic position
* pval: Association p-value
* coef: Effect size coefficient
* std_err: Standard error of coefficient
* stat: Test statistic
* alpha_value: Applied encoding parameter
* n_samples: Number of samples in analysis

Interpreting Alpha Values
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

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

See Also
--------

**Documentation:**

* `Documentation Home <index.html>`_ - Home
* :ref:`installation` - Installation instructions and requirements
* :ref:`quickstart` - Getting started guide with simple examples
* :ref:`user_guide` - Comprehensive user guide and tutorials
* :ref:`api_reference` - Complete API documentation
* :ref:`examples` - Example analyses and case studies
* :ref:`visualization` - Plotting and visualization guide
* :ref:`statistical_model` - Statistical methods and mathematical background
* :ref:`troubleshooting` - Troubleshooting guide and common issues
* :ref:`faq` - Frequently asked questions
* :ref:`citation` - How to cite EDGE in publications
* :ref:`changelog` - Version history and release notes
* :ref:`futureupdates` - Planned features and roadmap

---

*Last updated: 2025-12-25 for edge-gwas v0.1.1*

*For questions or issues, visit:* https://github.com/nicenzhou/edge-gwas/issues
