.. _quickstart:

Quick Start Guide
=================

This guide demonstrates a complete EDGE GWAS analysis in minutes.

Minimal Example
---------------

Fast minimal EDGE GWAS analysis:

.. code-block:: python

   from edge_gwas import EDGEAnalysis
   from edge_gwas.utils import load_plink_data, prepare_phenotype_data, stratified_train_test_split
   from edge_gwas.visualize import manhattan_plot, qq_plot
   
   # Load data
   geno, info = load_plink_data('data.bed', 'data.bim', 'data.fam')
   pheno = prepare_phenotype_data('pheno.txt', 'disease', ['age', 'sex', 'PC1', 'PC2'])
   
   # Split data
   train_g, test_g, train_p, test_p = stratified_train_test_split(
       geno, pheno, 'disease'
   )
   
   # Run analysis
   edge = EDGEAnalysis(outcome_type='binary')
   alpha_df, gwas_df = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       'disease', ['age', 'sex', 'PC1', 'PC2']
   )
   
   # Visualize
   manhattan_plot(gwas_df, 'manhattan.png')
   qq_plot(gwas_df, 'qq.png')

Step-by-Step Workflow
----------------------

1. Import and Initialize
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import pandas as pd
   from edge_gwas import EDGEAnalysis
   from edge_gwas.utils import (
       load_plink_data,
       prepare_phenotype_data,
       stratified_train_test_split
   )
   from edge_gwas.visualize import manhattan_plot, qq_plot
   
   # Initialize EDGE Analysis
   edge = EDGEAnalysis(
       outcome_type='binary',    # 'binary' or 'continuous'
       n_jobs=8,                 # Number of CPU cores
       verbose=True              # Print progress
   )

2. Load Data
~~~~~~~~~~~~

.. code-block:: python

   # Load genotype data
   genotype_data, variant_info = load_plink_data(
       bed_file='path/to/data.bed',
       bim_file='path/to/data.bim',
       fam_file='path/to/data.fam'
   )
   
   # Load phenotype data
   phenotype_df = prepare_phenotype_data(
       phenotype_file='path/to/phenotypes.txt',
       outcome_col='disease',
       covariate_cols=['age', 'sex', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5'],
       sample_id_col='IID'
   )

3. Split Data
~~~~~~~~~~~~~

.. code-block:: python

   train_geno, test_geno, train_pheno, test_pheno = stratified_train_test_split(
       genotype_df=genotype_data,
       phenotype_df=phenotype_df,
       outcome_col='disease',
       test_size=0.5,
       random_state=42
   )

4. Run EDGE Analysis
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # One-step complete analysis
   alpha_df, gwas_df = edge.run_full_analysis(
       train_genotype=train_geno,
       train_phenotype=train_pheno,
       test_genotype=test_geno,
       test_phenotype=test_pheno,
       outcome='disease',
       covariates=['age', 'sex', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5'],
       output_prefix='edge_results'
   )
   
   print(f"Tested {len(gwas_df)} variants")
   print(f"Significant (p < 5e-8): {(gwas_df['pval'] < 5e-8).sum()}")

5. Visualize Results
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Manhattan plot
   manhattan_plot(gwas_df, 'manhattan.png', title='EDGE GWAS')
   
   # QQ plot
   lambda_gc = qq_plot(gwas_df, 'qq_plot.png')
   print(f"Genomic inflation factor (λ): {lambda_gc:.3f}")

Understanding the Output
------------------------

Alpha Values DataFrame
~~~~~~~~~~~~~~~~~~~~~~

* ``variant_id``: SNP identifier
* ``alpha_value``: Encoding parameter (β_het/β_hom)
* ``eaf``: Effect allele frequency
* ``coef_het``: Heterozygous coefficient
* ``coef_hom``: Homozygous coefficient

GWAS Results DataFrame
~~~~~~~~~~~~~~~~~~~~~~

* ``variant_id``: SNP identifier
* ``chr``: Chromosome
* ``pos``: Position
* ``pval``: P-value
* ``coef``: Effect coefficient
* ``std_err``: Standard error
* ``alpha_value``: Applied encoding parameter

Next Steps
----------

* Read the :ref:`user_guide` for detailed information
* Check :ref:`examples` for complete workflows
* Explore :ref:`api_reference` for all functions
* Learn about the :ref:`statistical_model`
