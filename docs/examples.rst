.. _examples:

Example Workflows
=================

This page provides complete, real-world examples of EDGE GWAS analysis including new v0.1.1 features.

Example 1: Basic Binary Outcome Analysis with PCA
--------------------------------------------------

Complete workflow for case-control study with population structure control.

.. code-block:: python

   import pandas as pd
   from edge_gwas import EDGEAnalysis
   from edge_gwas.utils import (
       load_plink_data,
       prepare_phenotype_data,
       stratified_train_test_split,
       filter_variants_by_maf,
       calculate_pca_plink,
       attach_pcs_to_phenotype,
       get_pc_covariate_list,
       calculate_genomic_inflation
   )
   from edge_gwas.visualize import manhattan_plot, qq_plot
   
   # 1. Load data
   print("Loading genotype data...")
   genotype_df, variant_info = load_plink_data(
       bed_file='data/chr1.bed',
       bim_file='data/chr1.bim',
       fam_file='data/chr1.fam'
   )
   
   print("Loading phenotype data...")
   phenotype_df = prepare_phenotype_data(
       phenotype_file='data/phenotypes.txt',
       outcome_col='disease_status',
       covariate_cols=['age', 'sex'],
       sample_id_col='IID'
   )
   
   # 2. Calculate PCA for population structure (NEW in v0.1.1)
   print("Calculating principal components...")
   pca_df = calculate_pca_plink(
       plink_prefix='data/chr1',
       n_pcs=10,
       maf_threshold=0.01,
       ld_r2=0.2
   )
   
   # Attach PCs to phenotype
   phenotype_df = attach_pcs_to_phenotype(
       phenotype_df, pca_df, n_pcs=10
   )
   
   # 3. Quality control
   print("Applying QC filters...")
   genotype_df = filter_variants_by_maf(genotype_df, min_maf=0.01)
   
   # 4. Split data
   print("Splitting into train/test sets...")
   train_g, test_g, train_p, test_p = stratified_train_test_split(
       genotype_df,
       phenotype_df,
       outcome_col='disease_status',
       test_size=0.5,
       random_state=42
   )
   
   # 5. Run EDGE analysis
   print("Running EDGE analysis...")
   edge = EDGEAnalysis(outcome_type='binary', n_jobs=8, verbose=True)
   
   # Build covariates list including PCs
   covariates = ['age', 'sex'] + get_pc_covariate_list(10)
   
   alpha_df, gwas_df = edge.run_full_analysis(
       train_genotype=train_g,
       train_phenotype=train_p,
       test_genotype=test_g,
       test_phenotype=test_p,
       outcome='disease_status',
       covariates=covariates,
       variant_info=variant_info,
       output_prefix='chr1_edge_results'
   )
   
   # 6. Check results
   print(f"\nTotal variants tested: {len(gwas_df)}")
   print(f"Significant hits (p < 5e-8): {(gwas_df['pval'] < 5e-8).sum()}")
   
   lambda_gc = calculate_genomic_inflation(gwas_df['pval'])
   print(f"Genomic inflation factor: {lambda_gc:.3f}")
   
   # 7. Visualize
   manhattan_plot(gwas_df, 'chr1_manhattan.png', title='Chromosome 1 EDGE GWAS')
   qq_plot(gwas_df, 'chr1_qq.png')
   
   # 8. Extract top hits
   top_hits = gwas_df[gwas_df['pval'] < 5e-8].sort_values('pval')
   print(f"\nTop hits:\n{top_hits[['variant_id', 'chr', 'pos', 'pval', 'coef']]}")

Example 2: Quantitative Trait with Outcome Transformation
----------------------------------------------------------

Analysis of continuous phenotype with rank-based inverse normal transformation.

.. code-block:: python

   from edge_gwas import EDGEAnalysis
   from edge_gwas.utils import (
       load_plink_data,
       prepare_phenotype_data,
       calculate_pca_plink,
       attach_pcs_to_phenotype,
       stratified_train_test_split
   )
   
   # Load data
   genotype_df, variant_info = load_plink_data('data.bed', 'data.bim', 'data.fam')
   
   # Prepare quantitative phenotype (e.g., BMI)
   phenotype_df = prepare_phenotype_data(
       phenotype_file='phenotypes.txt',
       outcome_col='BMI',
       covariate_cols=['age', 'sex'],
       sample_id_col='IID'
   )
   
   # Calculate PCA
   pca_df = calculate_pca_plink('data', n_pcs=10)
   phenotype_df = attach_pcs_to_phenotype(phenotype_df, pca_df, n_pcs=10)
   
   # Split data
   train_g, test_g, train_p, test_p = stratified_train_test_split(
       genotype_df, phenotype_df, 'BMI', test_size=0.5, is_binary=False
   )
   
   # Run analysis with continuous outcome and transformation (NEW in v0.1.1)
   edge = EDGEAnalysis(
       outcome_type='continuous',
       outcome_transform='rank_inverse_normal',  # NEW: Apply RINT
       n_jobs=-1
   )
   
   alpha_df, gwas_df = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       outcome='BMI',
       covariates=['age', 'sex'] + [f'PC{i}' for i in range(1, 11)]
   )
   
   # Visualize
   from edge_gwas.visualize import manhattan_plot, qq_plot, plot_alpha_distribution
   manhattan_plot(gwas_df, 'bmi_manhattan.png')
   qq_plot(gwas_df, 'bmi_qq.png')
   plot_alpha_distribution(alpha_df, 'bmi_alpha_dist.png')

Example 3: Analysis with GRM for Related Samples
-------------------------------------------------

Control for population structure and relatedness using GRM.

.. code-block:: python

   from edge_gwas import EDGEAnalysis
   from edge_gwas.utils import (
       load_plink_data,
       prepare_phenotype_data,
       calculate_grm_gcta,
       load_grm_gcta,
       calculate_pca_pcair,
       attach_pcs_to_phenotype,
       stratified_train_test_split
   )
   
   # Load data
   genotype_df, variant_info = load_plink_data('data.bed', 'data.bim', 'data.fam')
   phenotype_df = prepare_phenotype_data(
       'phenotypes.txt', 'disease', ['age', 'sex']
   )
   
   # Calculate GRM for population structure control (NEW in v0.1.1)
   print("Calculating GRM...")
   grm_prefix = calculate_grm_gcta(
       plink_prefix='data',
       output_prefix='output/grm',
       maf_threshold=0.01,
       max_threads=8
   )
   
   # Load GRM
   grm_matrix, grm_ids = load_grm_gcta('output/grm')
   
   # Calculate PC-AiR (accounts for relatedness)  (NEW in v0.1.1)
   print("Calculating PC-AiR...")
   pca_df = calculate_pca_pcair(
       plink_prefix='data',
       n_pcs=10,
       kinship_matrix='output/grm',
       kin_threshold=0.0884
   )
   
   # Attach PCs to phenotype
   phenotype_df = attach_pcs_to_phenotype(phenotype_df, pca_df, n_pcs=10)
   
   # Split data
   train_g, test_g, train_p, test_p = stratified_train_test_split(
       genotype_df, phenotype_df, 'disease', test_size=0.5
   )
   
   # Run EDGE with GRM (NEW in v0.1.1)
   edge = EDGEAnalysis(outcome_type='binary', n_jobs=8)
   
   alpha_df, gwas_df = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       outcome='disease',
       covariates=['age', 'sex'] + [f'PC{i}' for i in range(1, 11)],
       grm_matrix=grm_matrix,
       grm_sample_ids=grm_ids,
       output_prefix='results/edge_with_grm'
   )
   
   print(f"Analysis complete with GRM control")

Example 4: Multi-Format Data Loading
-------------------------------------

Load data from different genetic file formats.

.. code-block:: python

   from edge_gwas.utils import (
       load_plink_data,
       load_pgen_data,
       load_bgen_data,
       load_vcf_data
   )
   
   # Option 1: Load PLINK format (.bed/.bim/.fam)
   print("Loading PLINK format...")
   geno1, info1 = load_plink_data(
       bed_file='data.bed',
       bim_file='data.bim',
       fam_file='data.fam'
   )
   
   # Option 2: Load PLINK 2 format (.pgen/.pvar/.psam) - NEW in v0.1.1
   print("Loading PGEN format...")
   geno2, info2 = load_pgen_data(
       pgen_file='data.pgen',
       pvar_file='data.pvar',
       psam_file='data.psam'
   )
   
   # Option 3: Load BGEN format - NEW in v0.1.1
   print("Loading BGEN format...")
   geno3, info3 = load_bgen_data(
       bgen_file='data.bgen',
       sample_file='data.sample'  # Optional
   )
   
   # Option 4: Load VCF format - NEW in v0.1.1
   print("Loading VCF format...")
   geno4, info4 = load_vcf_data(
       vcf_file='data.vcf.gz',
       dosage=True  # Use dosages (DS field) if available
   )
   
   # All formats return the same structure:
   # - genotype_df: samples × variants DataFrame
   # - variant_info: variant information DataFrame
   
   print(f"PLINK: {geno1.shape}")
   print(f"PGEN: {geno2.shape}")
   print(f"BGEN: {geno3.shape}")
   print(f"VCF: {geno4.shape}")
   
   # Continue with standard EDGE workflow using any format
   from edge_gwas import EDGEAnalysis
   edge = EDGEAnalysis(outcome_type='binary')
   # ... rest of analysis

Example 5: Comprehensive Quality Control Pipeline
--------------------------------------------------

Detailed QC with all available filters.

.. code-block:: python

   from edge_gwas import EDGEAnalysis
   from edge_gwas.utils import (
       load_plink_data,
       prepare_phenotype_data,
       filter_variants_by_maf,
       filter_variants_by_missing,
       filter_variants_by_hwe,
       filter_samples_by_call_rate,
       check_case_control_balance,
       identify_related_samples,
       filter_related_samples,
       calculate_grm_gcta,
       load_grm_gcta,
       stratified_train_test_split
   )
   
   # Load data
   genotype_df, variant_info = load_plink_data('data.bed', 'data.bim', 'data.fam')
   phenotype_df = prepare_phenotype_data(
       'phenotypes.txt', 'disease', ['age', 'sex']
   )
   
   print(f"Starting: {genotype_df.shape[1]} variants, {genotype_df.shape[0]} samples")
   
   # QC Step 1: Filter by MAF
   genotype_df = filter_variants_by_maf(genotype_df, min_maf=0.01, verbose=True)
   
   # QC Step 2: Filter by missingness
   genotype_df = filter_variants_by_missing(genotype_df, max_missing=0.05, verbose=True)
   
   # QC Step 3: Filter by HWE (NEW in v0.1.1)
   genotype_df = filter_variants_by_hwe(genotype_df, hwe_threshold=1e-6, verbose=True)
   
   # QC Step 4: Filter samples by call rate (NEW in v0.1.1)
   genotype_df, phenotype_df = filter_samples_by_call_rate(
       genotype_df, phenotype_df, min_call_rate=0.95, verbose=True
   )
   
   # QC Step 5: Check case/control balance (NEW in v0.1.1)
   balance = check_case_control_balance(phenotype_df, 'disease', verbose=True)
   
   if balance['ratio'] < 0.1 or balance['ratio'] > 10:
       print(f"Warning: Imbalanced case/control ratio: {balance['ratio']:.2f}")
   
   # QC Step 6: Check for related samples (NEW in v0.1.1)
   print("\nChecking for related samples...")
   grm_prefix = calculate_grm_gcta('data', output_prefix='qc/grm')
   grm_matrix, grm_ids = load_grm_gcta('qc/grm')
   
   related_pairs = identify_related_samples(
       grm_matrix, grm_ids, threshold=0.0884, verbose=True
   )
   
   if len(related_pairs) > 0:
       print(f"Found {len(related_pairs)} related pairs")
       # Option: Remove related samples
       phenotype_df = filter_related_samples(
           phenotype_df, grm_matrix, grm_ids,
           threshold=0.0884, method='greedy', verbose=True
       )
       genotype_df = genotype_df.loc[phenotype_df.index]
   
   print(f"\nAfter QC: {genotype_df.shape[1]} variants, {genotype_df.shape[0]} samples")
   
   # Continue with analysis
   train_g, test_g, train_p, test_p = stratified_train_test_split(
       genotype_df, phenotype_df, 'disease', test_size=0.5
   )
   
   edge = EDGEAnalysis(outcome_type='binary', n_jobs=8)
   alpha_df, gwas_df = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       outcome='disease',
       covariates=['age', 'sex']
   )

Example 6: Comparing EDGE with Additive GWAS
---------------------------------------------

Compare EDGE results with traditional additive model.

.. code-block:: python

   import pandas as pd
   import numpy as np
   from edge_gwas import EDGEAnalysis
   from edge_gwas.utils import additive_gwas
   import matplotlib.pyplot as plt
   
   # Run EDGE analysis
   edge = EDGEAnalysis(outcome_type='binary', n_jobs=8)
   alpha_df, edge_gwas_results = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       outcome='disease',
       covariates=['age', 'sex'] + [f'PC{i}' for i in range(1, 11)]
   )
   
   # Run standard additive GWAS (NEW in v0.1.1)
   print("Running additive GWAS for comparison...")
   additive_results = additive_gwas(
       genotype_df=test_g,
       phenotype_df=test_p,
       outcome='disease',
       covariates=['age', 'sex'] + [f'PC{i}' for i in range(1, 11)],
       outcome_type='binary'
   )
   
   # Merge results for comparison
   comparison = edge_gwas_results[['variant_id', 'pval', 'coef']].merge(
       additive_results[['variant_id', 'pval', 'coef']],
       on='variant_id',
       suffixes=('_edge', '_additive')
   )
   
   # Add alpha values
   comparison = comparison.merge(
       alpha_df[['variant_id', 'alpha_value']],
       on='variant_id',
       how='left'
   )
   
   # Calculate -log10 p-values
   comparison['-log10p_edge'] = -np.log10(comparison['pval_edge'])
   comparison['-log10p_additive'] = -np.log10(comparison['pval_additive'])
   comparison['edge_advantage'] = comparison['-log10p_edge'] - comparison['-log10p_additive']
   
   # Find variants where EDGE shows advantage
   edge_better = comparison[comparison['edge_advantage'] > 2].sort_values(
       'edge_advantage', ascending=False
   )
   
   print(f"\nVariants where EDGE shows >100-fold p-value improvement:")
   print(edge_better[['variant_id', 'pval_edge', 'pval_additive', 
                       'alpha_value', 'edge_advantage']].head(10))
   
   # Analyze by inheritance pattern
   def categorize_alpha(alpha):
       if pd.isna(alpha):
           return 'Unknown'
       elif alpha < 0.3:
           return 'Recessive'
       elif alpha > 0.7:
           return 'Dominant'
       else:
           return 'Additive'
   
   comparison['inheritance'] = comparison['alpha_value'].apply(categorize_alpha)
   
   print("\nEDGE advantage by inheritance pattern:")
   print(comparison.groupby('inheritance')['edge_advantage'].agg(['mean', 'median', 'std']))
   
   # Visualize comparison
   fig, axes = plt.subplots(2, 2, figsize=(14, 12))
   
   # Scatter plot
   axes[0, 0].scatter(comparison['-log10p_additive'], comparison['-log10p_edge'], 
                      alpha=0.5, s=10, c='blue')
   axes[0, 0].plot([0, comparison['-log10p_additive'].max()], 
                   [0, comparison['-log10p_additive'].max()], 
                   'r--', label='y=x')
   axes[0, 0].set_xlabel('Additive -log10(p)')
   axes[0, 0].set_ylabel('EDGE -log10(p)')
   axes[0, 0].set_title('EDGE vs Additive GWAS')
   axes[0, 0].legend()
   axes[0, 0].grid(True, alpha=0.3)
   
   # Advantage distribution
   axes[0, 1].hist(comparison['edge_advantage'], bins=50, edgecolor='black', alpha=0.7)
   axes[0, 1].axvline(0, color='red', linestyle='--', label='No difference')
   axes[0, 1].set_xlabel('EDGE Advantage (-log10p difference)')
   axes[0, 1].set_ylabel('Number of variants')
   axes[0, 1].set_title('Distribution of EDGE Advantage')
   axes[0, 1].legend()
   axes[0, 1].grid(True, alpha=0.3, axis='y')
   
   # Alpha vs advantage
   axes[1, 0].scatter(comparison['alpha_value'], comparison['edge_advantage'], 
                      alpha=0.5, s=10)
   axes[1, 0].axhline(0, color='red', linestyle='--')
   axes[1, 0].axvline(0.5, color='green', linestyle='--', label='Additive (α=0.5)')
   axes[1, 0].set_xlabel('Alpha value')
   axes[1, 0].set_ylabel('EDGE Advantage')
   axes[1, 0].set_title('Alpha Value vs EDGE Advantage')
   axes[1, 0].legend()
   axes[1, 0].grid(True, alpha=0.3)
   
   # Boxplot by inheritance
   inheritance_order = ['Recessive', 'Additive', 'Dominant']
   comparison_clean = comparison[comparison['inheritance'].isin(inheritance_order)]
   data_to_plot = [comparison_clean[comparison_clean['inheritance'] == i]['edge_advantage'].values 
                    for i in inheritance_order]
   axes[1, 1].boxplot(data_to_plot, labels=inheritance_order)
   axes[1, 1].axhline(0, color='red', linestyle='--')
   axes[1, 1].set_ylabel('EDGE Advantage')
   axes[1, 1].set_title('EDGE Advantage by Inheritance Pattern')
   axes[1, 1].grid(True, alpha=0.3, axis='y')
   
   plt.tight_layout()
   plt.savefig('edge_vs_additive_comprehensive.png', dpi=300)
   plt.close()
   
   print("\nComparison plot saved to 'edge_vs_additive_comprehensive.png'")

Example 7: Cross-Validation for Model Stability
------------------------------------------------

Use k-fold cross-validation to assess alpha stability.

.. code-block:: python

   import pandas as pd
   from edge_gwas.utils import cross_validated_edge_analysis
   
   # Run cross-validated analysis (NEW in v0.1.1)
   print("Running 5-fold cross-validation...")
   avg_alpha, meta_gwas, all_alpha, all_gwas = cross_validated_edge_analysis(
       genotype_df=genotype_df,
       phenotype_df=phenotype_df,
       outcome='disease',
       covariates=['age', 'sex'] + [f'PC{i}' for i in range(1, 11)],
       outcome_type='binary',
       n_folds=5,
       n_jobs=8,
       random_state=42
   )
   
   # Assess alpha stability
   print("\nAlpha value stability across folds:")
   print(f"Mean alpha std: {avg_alpha['alpha_std'].mean():.3f}")
   print(f"Median alpha std: {avg_alpha['alpha_std'].median():.3f}")
   
   # Find unstable variants
   unstable_variants = avg_alpha[avg_alpha['alpha_std'] > 0.3]
   print(f"\nVariants with unstable alpha (std > 0.3): {len(unstable_variants)}")
   
   if len(unstable_variants) > 0:
       print("\nMost unstable variants:")
       print(unstable_variants.nlargest(10, 'alpha_std')[
           ['variant_id', 'alpha_mean', 'alpha_std']
       ])
   
   # Find stable and significant variants
   stable_significant = meta_gwas[
       (meta_gwas['pval'] < 5e-8) &
       (meta_gwas['variant_id'].isin(avg_alpha[avg_alpha['alpha_std'] < 0.2]['variant_id']))
   ]
   
   print(f"\nStable and significant variants: {len(stable_significant)}")
   
   # Save results
   avg_alpha.to_csv('cv_averaged_alpha.txt', sep='\t', index=False)
   meta_gwas.to_csv('cv_meta_gwas.txt', sep='\t', index=False)
   all_alpha.to_csv('cv_all_alpha_folds.txt', sep='\t', index=False)
   all_gwas.to_csv('cv_all_gwas_folds.txt', sep='\t', index=False)
   
   # Visualize
   from edge_gwas.visualize import manhattan_plot, qq_plot
   manhattan_plot(meta_gwas, 'cv_manhattan.png', title='Cross-Validated EDGE GWAS')
   qq_plot(meta_gwas, 'cv_qq.png')
   
   # Plot alpha stability
   import matplotlib.pyplot as plt
   
   fig, axes = plt.subplots(1, 2, figsize=(14, 6))
   
   # Alpha stability distribution
   axes[0].hist(avg_alpha['alpha_std'], bins=50, edgecolor='black', alpha=0.7)
   axes[0].axvline(0.3, color='red', linestyle='--', label='Instability threshold')
   axes[0].set_xlabel('Alpha Standard Deviation')
   axes[0].set_ylabel('Number of Variants')
   axes[0].set_title('Alpha Stability Across Folds')
   axes[0].legend()
   axes[0].grid(True, alpha=0.3, axis='y')
   
   # Alpha mean vs std
   axes[1].scatter(avg_alpha['alpha_mean'], avg_alpha['alpha_std'], alpha=0.5, s=10)
   axes[1].axhline(0.3, color='red', linestyle='--', label='Instability threshold')
   axes[1].set_xlabel('Mean Alpha Value')
   axes[1].set_ylabel('Alpha Standard Deviation')
   axes[1].set_title('Alpha Mean vs Variability')
   axes[1].legend()
   axes[1].grid(True, alpha=0.3)
   
   plt.tight_layout()
   plt.savefig('alpha_stability.png', dpi=300)
   plt.close()

Example 8: Complete Analysis with All v0.1.1 Features
------------------------------------------------------

Comprehensive workflow using GRM, PCA, outcome transformation, and quality control.

.. code-block:: python

   from edge_gwas import EDGEAnalysis
   from edge_gwas.utils import (
       load_plink_data,
       prepare_phenotype_data,
       filter_variants_by_maf,
       filter_variants_by_missing,
       filter_variants_by_hwe,
       filter_samples_by_call_rate,
       calculate_grm_gcta,
       load_grm_gcta,
       calculate_pca_pcair,
       attach_pcs_to_phenotype,
       get_pc_covariate_list,
       stratified_train_test_split
   )
   from edge_gwas.visualize import manhattan_plot, qq_plot, plot_alpha_distribution
   
   print("="*70)
   print("EDGE-GWAS v0.1.1 - Complete Analysis Pipeline")
   print("="*70)
   
   # Step 1: Load data
   print("\n[1/9] Loading genotype data...")
   genotype_df, variant_info = load_plink_data('data.bed', 'data.bim', 'data.fam')
   print(f"  Loaded: {genotype_df.shape[0]} samples, {genotype_df.shape[1]} variants")
   
   print("\n[2/9] Loading phenotype data...")
   phenotype_df = prepare_phenotype_data(
       'phenotypes.txt',
       outcome_col='trait',
       covariate_cols=['age', 'sex', 'batch'],
       sample_id_col='IID'
   )
   print(f"  Loaded: {len(phenotype_df)} samples")
   
   # Step 2: Quality control
   print("\n[3/9] Applying quality control filters...")
   genotype_df = filter_variants_by_maf(genotype_df, min_maf=0.01, verbose=True)
   genotype_df = filter_variants_by_missing(genotype_df, max_missing=0.05, verbose=True)
   genotype_df = filter_variants_by_hwe(genotype_df, hwe_threshold=1e-6, verbose=True)
   genotype_df, phenotype_df = filter_samples_by_call_rate(
       genotype_df, phenotype_df, min_call_rate=0.95, verbose=True
   )
   
   # Step 3: Calculate GRM
   print("\n[4/9] Calculating genetic relationship matrix...")
   grm_prefix = calculate_grm_gcta(
       plink_prefix='data',
       output_prefix='output/grm',
       maf_threshold=0.01,
       max_threads=8,
       verbose=True
   )
   grm_matrix, grm_ids = load_grm_gcta('output/grm', verbose=True)
   
   # Step 4: Calculate PCA with PC-AiR (accounts for relatedness)
   print("\n[5/9] Calculating PC-AiR (accounts for relatedness)...")
   pca_df = calculate_pca_pcair(
       plink_prefix='data',
       n_pcs=10,
       kinship_matrix='output/grm',
       kin_threshold=0.0884,
       verbose=True
   )
   
   # Attach PCs to phenotype
   phenotype_df = attach_pcs_to_phenotype(
       phenotype_df, pca_df, n_pcs=10, verbose=True
   )
   
   # Step 5: Split data
   print("\n[6/9] Splitting data into train/test sets...")
   train_g, test_g, train_p, test_p = stratified_train_test_split(
       genotype_df, phenotype_df,
       outcome_col='trait',
       test_size=0.5,
       is_binary=False,
       random_state=42
   )
   
   # Step 6: Run EDGE analysis with all features
   print("\n[7/9] Running EDGE analysis...")
   print("  Features enabled:")
   print("    - Outcome transformation: rank_inverse_normal")
   print("    - GRM control: Yes")
   print("    - Population structure PCs: 10")
   
   edge = EDGEAnalysis(
       outcome_type='continuous',
       outcome_transform='rank_inverse_normal',
       n_jobs=-1,
       verbose=True
   )
   
   covariates = ['age', 'sex', 'batch'] + get_pc_covariate_list(10)
   
   alpha_df, gwas_df = edge.run_full_analysis(
       train_genotype=train_g,
       train_phenotype=train_p,
       test_genotype=test_g,
       test_phenotype=test_p,
       outcome='trait',
       covariates=covariates,
       variant_info=variant_info,
       grm_matrix=grm_matrix,
       grm_sample_ids=grm_ids,
       output_prefix='output/edge_complete'
   )
   
   # Step 7: Generate visualizations
   print("\n[8/9] Generating visualizations...")
   manhattan_plot(
       gwas_df,
       output='output/manhattan.png',
       title='EDGE GWAS - Complete Analysis',
       sig_threshold=5e-8,
       suggestive_threshold=1e-5
   )
   
   from edge_gwas.utils import calculate_genomic_inflation
   lambda_gc = qq_plot(
       gwas_df,
       output='output/qq_plot.png',
       title='QQ Plot - EDGE GWAS'
   )
   
   plot_alpha_distribution(
       alpha_df,
       output='output/alpha_distribution.png',
       title='Alpha Value Distribution'
   )
   
   # Step 8: Summarize results
   print("\n[9/9] Analysis Summary")
   print("="*70)
   print(f"Total variants tested: {len(gwas_df):,}")
   print(f"Genome-wide significant (p < 5e-8): {(gwas_df['pval'] < 5e-8).sum()}")
   print(f"Suggestive (p < 1e-5): {((gwas_df['pval'] >= 5e-8) & (gwas_df['pval'] < 1e-5)).sum()}")
   print(f"Genomic inflation factor (λ): {lambda_gc:.3f}")
   
   print(f"\nAlpha statistics:")
   print(f"  Mean: {alpha_df['alpha_value'].mean():.3f}")
   print(f"  Median: {alpha_df['alpha_value'].median():.3f}")
   print(f"  Std: {alpha_df['alpha_value'].std():.3f}")
   
   # Categorize by inheritance
   def categorize_alpha(alpha):
       if pd.isna(alpha):
           return 'Unknown'
       elif alpha < 0.3:
           return 'Recessive'
       elif 0.3 <= alpha < 0.7:
           return 'Additive'
       else:
           return 'Dominant'
   
   alpha_df['inheritance'] = alpha_df['alpha_value'].apply(categorize_alpha)
   print(f"\nInheritance pattern distribution:")
   print(alpha_df['inheritance'].value_counts())
   
   # Top hits
   top_hits = gwas_df[gwas_df['pval'] < 5e-8].sort_values('pval')
   if len(top_hits) > 0:
       print(f"\nTop {min(10, len(top_hits))} significant hits:")
       print(top_hits[['variant_id', 'chr', 'pos', 'pval', 'coef', 'alpha_value']].head(10))
   
   print("\nOutput files:")
   print("  - output/edge_complete_alpha_values.csv")
   print("  - output/edge_complete_gwas_results.csv")
   print("  - output/manhattan.png")
   print("  - output/qq_plot.png")
   print("  - output/alpha_distribution.png")
   print("="*70)

Example 9: Multi-Chromosome Genome-Wide Analysis
-------------------------------------------------

Efficient genome-wide analysis processing chromosomes separately.

.. code-block:: python

   import pandas as pd
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
   
   # Load phenotype once
   phenotype_df = prepare_phenotype_data(
       'phenotypes.txt',
       outcome_col='disease',
       covariate_cols=['age', 'sex'],
       sample_id_col='IID'
   )
   
   # Calculate PCA from chromosome 1 (representative)
   print("Calculating PCA from chromosome 1...")
   pca_df = calculate_pca_plink('data/chr1', n_pcs=10)
   phenotype_df = attach_pcs_to_phenotype(phenotype_df, pca_df, n_pcs=10)
   
   # Initialize EDGE
   edge = EDGEAnalysis(outcome_type='binary', n_jobs=8, verbose=True)
   covariates = ['age', 'sex'] + get_pc_covariate_list(10)
   
   all_alpha = []
   all_gwas = []
   
   # Process each chromosome
   for chrom in range(1, 23):
       print(f"\n{'='*70}")
       print(f"Processing Chromosome {chrom}")
       print(f"{'='*70}")
       
       try:
           # Load chromosome data
           genotype_df, variant_info = load_plink_data(
               bed_file=f'data/chr{chrom}.bed',
               bim_file=f'data/chr{chrom}.bim',
               fam_file=f'data/chr{chrom}.fam'
           )
           
           # QC
           genotype_df = filter_variants_by_maf(genotype_df, min_maf=0.01)
           
           # Split data
           train_g, test_g, train_p, test_p = stratified_train_test_split(
               genotype_df, phenotype_df,
               outcome_col='disease',
               test_size=0.5,
               random_state=42
           )
           
           # Run EDGE
           alpha_df, gwas_df = edge.run_full_analysis(
               train_g, train_p, test_g, test_p,
               outcome='disease',
               covariates=covariates,
               variant_info=variant_info,
               output_prefix=f'output/chr{chrom}_edge'
           )
           
           # Store results
           all_alpha.append(alpha_df)
           all_gwas.append(gwas_df)
           
           print(f"Chromosome {chrom} complete: {len(gwas_df):,} variants tested")
           sig = (gwas_df['pval'] < 5e-8).sum()
           if sig > 0:
               print(f"  Significant hits: {sig}")
           
           # Clean up memory
           del genotype_df, train_g, test_g
           import gc
           gc.collect()
           
       except Exception as e:
           print(f"Error processing chromosome {chrom}: {e}")
           continue
   
   # Combine all results
   print(f"\n{'='*70}")
   print("Combining results from all chromosomes...")
   print(f"{'='*70}")
   
   combined_alpha = pd.concat(all_alpha, ignore_index=True)
   combined_gwas = pd.concat(all_gwas, ignore_index=True)
   
   # Save combined results
   combined_alpha.to_csv('output/genome_wide_alpha.txt', sep='\t', index=False)
   combined_gwas.to_csv('output/genome_wide_gwas.txt', sep='\t', index=False)
   
   # Generate genome-wide visualizations
   manhattan_plot(
       combined_gwas,
       output='output/genome_wide_manhattan.png',
       title='Genome-Wide EDGE GWAS'
   )
   
   lambda_gc = qq_plot(
       combined_gwas,
       output='output/genome_wide_qq.png',
       title='Genome-Wide QQ Plot'
   )
   
   # Summary
   total_variants = len(combined_gwas)
   sig_hits = combined_gwas[combined_gwas['pval'] < 5e-8]
   
   print(f"\nGenome-Wide Analysis Summary:")
   print(f"  Total variants tested: {total_variants:,}")
   print(f"  Significant hits (p < 5e-8): {len(sig_hits)}")
   print(f"  Genomic inflation (λ): {lambda_gc:.3f}")
   
   if len(sig_hits) > 0:
       print(f"\n  Chromosomes with significant hits:")
       chr_counts = sig_hits['chr'].value_counts().sort_index()
       for chr_num, count in chr_counts.items():
           print(f"    Chr {chr_num}: {count} hits")
       
       print(f"\n  Top 10 genome-wide hits:")
       print(sig_hits.nlargest(10, 'pval', keep='first')[
           ['variant_id', 'chr', 'pos', 'pval', 'coef', 'alpha_value']
       ])

Example 10: Batch Processing for Very Large Datasets
-----------------------------------------------------

Memory-efficient processing for biobank-scale data.

.. code-block:: python

   import pandas as pd
   import gc
   from edge_gwas import EDGEAnalysis
   from edge_gwas.utils import (
       load_plink_data,
       prepare_phenotype_data,
       filter_variants_by_maf,
       stratified_train_test_split
   )
   
   def process_chromosome_batch(chrom_list, phenotype_df, output_dir):
       """Process a batch of chromosomes efficiently."""
       
       edge = EDGEAnalysis(outcome_type='binary', n_jobs=-1, verbose=True)
       batch_results = []
       
       for chrom in chrom_list:
           print(f"\nProcessing chromosome {chrom}...")
           
           # Load chromosome
           geno, info = load_plink_data(
               f'data/chr{chrom}.bed',
               f'data/chr{chrom}.bim',
               f'data/chr{chrom}.fam'
           )
           
           # QC and split
           geno = filter_variants_by_maf(geno, min_maf=0.01, verbose=False)
           train_g, test_g, train_p, test_p = stratified_train_test_split(
               geno, phenotype_df, 'disease', test_size=0.5, random_state=42
           )
           
           # Analyze
           alpha_df, gwas_df = edge.run_full_analysis(
               train_g, train_p, test_g, test_p,
               outcome='disease',
               covariates=['age', 'sex'] + [f'PC{i}' for i in range(1, 11)],
               variant_info=info,
               output_prefix=f'{output_dir}/chr{chrom}'
           )
           
           batch_results.append({
               'chrom': chrom,
               'n_variants': len(gwas_df),
               'n_significant': (gwas_df['pval'] < 5e-8).sum()
           })
           
           # Aggressive memory cleanup
           del geno, info, train_g, test_g, train_p, test_p, alpha_df, gwas_df
           gc.collect()
           
           print(f"  Completed: {batch_results[-1]['n_variants']:,} variants")
       
       return pd.DataFrame(batch_results)
   
   # Load phenotype
   phenotype_df = prepare_phenotype_data(
       'phenotypes.txt',
       outcome_col='disease',
       covariate_cols=['age', 'sex'] + [f'PC{i}' for i in range(1, 11)]
   )
   
   # Process in small batches to minimize memory usage
   batches = [
       [1, 2, 3, 4, 5],
       [6, 7, 8, 9, 10],
       [11, 12, 13, 14, 15],
       [16, 17, 18, 19, 20],
       [21, 22]
   ]
   
   all_batch_results = []
   
   for i, batch in enumerate(batches, 1):
       print(f"\n{'='*70}")
       print(f"Processing Batch {i}/{len(batches)}: Chromosomes {batch}")
       print(f"{'='*70}")
       
       batch_results = process_chromosome_batch(batch, phenotype_df, 'output')
       all_batch_results.append(batch_results)
       
       # Clean up between batches
       gc.collect()
   
   # Combine batch summaries
   summary = pd.concat(all_batch_results, ignore_index=True)
   summary.to_csv('output/batch_processing_summary.txt', sep='\t', index=False)
   
   print(f"\n{'='*70}")
   print("Batch Processing Complete")
   print(f"{'='*70}")
   print(summary)
   print(f"\nTotal variants: {summary['n_variants'].sum():,}")
   print(f"Total significant: {summary['n_significant'].sum()}")

Example 11: Using Pre-Calculated Alpha Values
----------------------------------------------

Apply alpha values from one cohort to another.

.. code-block:: python

   import pandas as pd
   from edge_gwas import EDGEAnalysis
   from edge_gwas.io_handlers import load_alpha_values
   from edge_gwas.utils import load_plink_data, prepare_phenotype_data
   from edge_gwas.visualize import manhattan_plot, qq_plot
   
   # Load pre-calculated alpha values from discovery cohort
   print("Loading pre-calculated alpha values...")
   alpha_df = load_alpha_values('discovery_cohort_alpha_values.txt')
   print(f"  Loaded alpha values for {len(alpha_df):,} variants")
   
   # Load replication cohort data
   print("\nLoading replication cohort genotype data...")
   test_genotype, variant_info = load_plink_data(
       'replication_cohort.bed',
       'replication_cohort.bim',
       'replication_cohort.fam'
   )
   
   print("Loading replication cohort phenotype data...")
   test_phenotype = prepare_phenotype_data(
       'replication_phenotypes.txt',
       outcome_col='disease',
       covariate_cols=['age', 'sex'] + [f'PC{i}' for i in range(1, 11)]
   )
   
   # Find common variants between discovery and replication
   common_variants = set(alpha_df['variant_id']) & set(test_genotype.columns)
   print(f"\nCommon variants: {len(common_variants):,}")
   
   # Subset to common variants
   test_genotype = test_genotype[list(common_variants)]
   alpha_df = alpha_df[alpha_df['variant_id'].isin(common_variants)]
   
   # Apply alpha values to replication cohort
   print("\nApplying EDGE encoding to replication cohort...")
   edge = EDGEAnalysis(outcome_type='binary', n_jobs=8, verbose=True)
   
   gwas_df = edge.apply_alpha(
       genotype_data=test_genotype,
       phenotype_df=test_phenotype,
       outcome='disease',
       covariates=['age', 'sex'] + [f'PC{i}' for i in range(1, 11)],
       alpha_values=alpha_df
   )
   
   # Save results
   gwas_df.to_csv('replication_gwas_results.txt', sep='\t', index=False)
   
   # Visualize
   manhattan_plot(gwas_df, 'replication_manhattan.png', 
                  title='Replication Cohort - EDGE GWAS')
   lambda_gc = qq_plot(gwas_df, 'replication_qq.png')
   
   # Summary
   print(f"\nReplication Analysis Summary:")
   print(f"  Variants tested: {len(gwas_df):,}")
   print(f"  Significant hits (p < 5e-8): {(gwas_df['pval'] < 5e-8).sum()}")
   print(f"  Genomic inflation (λ): {lambda_gc:.3f}")
   
   # Compare with discovery cohort if available
   if 'discovery_gwas_results.txt' in os.listdir('.'):
       discovery_gwas = pd.read_csv('discovery_gwas_results.txt', sep='\t')
       
       # Merge discovery and replication results
       comparison = discovery_gwas[['variant_id', 'pval']].merge(
           gwas_df[['variant_id', 'pval']],
           on='variant_id',
           suffixes=('_discovery', '_replication')
       )
       
       # Find variants significant in both
       both_sig = comparison[
           (comparison['pval_discovery'] < 5e-8) &
           (comparison['pval_replication'] < 0.05)
       ]
       
       print(f"\nVariants significant in discovery and replicated (p < 0.05):")
       print(f"  Count: {len(both_sig)}")
       
       if len(both_sig) > 0:
           print("\nTop replicated variants:")
           print(both_sig.sort_values('pval_replication')[
               ['variant_id', 'pval_discovery', 'pval_replication']
           ].head(10))

Example 12: Detailed Alpha Interpretation and Reporting
--------------------------------------------------------

Deep dive into alpha values and their biological interpretation.

.. code-block:: python

   import pandas as pd
   import numpy as np
   import matplotlib.pyplot as plt
   from edge_gwas import EDGEAnalysis
   from edge_gwas.visualize import plot_alpha_distribution
   
   # After running EDGE analysis
   edge = EDGEAnalysis(outcome_type='binary', n_jobs=8)
   alpha_df, gwas_df = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       outcome='disease',
       covariates=['age', 'sex'] + [f'PC{i}' for i in range(1, 11)]
   )
   
   # Define inheritance model categories
   def categorize_alpha(alpha):
       """Categorize alpha by inheritance model with detailed bins."""
       if pd.isna(alpha):
           return 'Unknown'
       elif alpha < 0:
           return 'Under-recessive'
       elif 0 <= alpha < 0.2:
           return 'Strong recessive'
       elif 0.2 <= alpha < 0.4:
           return 'Recessive'
       elif 0.4 <= alpha < 0.6:
           return 'Additive'
       elif 0.6 <= alpha < 0.8:
           return 'Dominant'
       elif 0.8 <= alpha <= 1.2:
           return 'Strong dominant'
       else:
           return 'Over-dominant'
   
   alpha_df['inheritance_detailed'] = alpha_df['alpha_value'].apply(categorize_alpha)
   
   # Summary statistics
   print("\n" + "="*70)
   print("ALPHA VALUE INTERPRETATION REPORT")
   print("="*70)
   
   print("\n1. OVERALL DISTRIBUTION")
   print("-" * 70)
   print(f"Total variants: {len(alpha_df):,}")
   print(f"Mean alpha: {alpha_df['alpha_value'].mean():.3f}")
   print(f"Median alpha: {alpha_df['alpha_value'].median():.3f}")
   print(f"Std alpha: {alpha_df['alpha_value'].std():.3f}")
   print(f"Min alpha: {alpha_df['alpha_value'].min():.3f}")
   print(f"Max alpha: {alpha_df['alpha_value'].max():.3f}")
   
   print("\n2. INHERITANCE MODEL DISTRIBUTION")
   print("-" * 70)
   inheritance_counts = alpha_df['inheritance_detailed'].value_counts().sort_index()
   for model, count in inheritance_counts.items():
       pct = (count / len(alpha_df)) * 100
       print(f"{model:20s}: {count:6,} ({pct:5.1f}%)")
   
   # Analyze by effect allele frequency
   print("\n3. ALPHA BY EFFECT ALLELE FREQUENCY")
   print("-" * 70)
   
   eaf_bins = [0, 0.05, 0.1, 0.25, 0.5, 1.0]
   eaf_labels = ['<5%', '5-10%', '10-25%', '25-50%', '>50%']
   alpha_df['eaf_bin'] = pd.cut(alpha_df['eaf'], bins=eaf_bins, labels=eaf_labels)
   
   eaf_summary = alpha_df.groupby('eaf_bin')['alpha_value'].agg(['count', 'mean', 'std'])
   print(eaf_summary)
   
   # Find extreme variants
   print("\n4. EXTREME ALPHA VALUES")
   print("-" * 70)
   
   # Most recessive
   most_recessive = alpha_df.nsmallest(10, 'alpha_value')
   print("\nMost recessive variants (lowest alpha):")
   print(most_recessive[['variant_id', 'alpha_value', 'eaf', 'coef_het', 'coef_hom']])
   
   # Most dominant
   most_dominant = alpha_df.nlargest(10, 'alpha_value')
   print("\nMost dominant variants (highest alpha):")
   print(most_dominant[['variant_id', 'alpha_value', 'eaf', 'coef_het', 'coef_hom']])
   
   # Merge with GWAS results to find significant nonadditive effects
   from edge_gwas.utils import merge_alpha_with_gwas
   merged = merge_alpha_with_gwas(gwas_df, alpha_df)
   
   print("\n5. SIGNIFICANT NONADDITIVE EFFECTS")
   print("-" * 70)
   
   # Nonadditive: alpha < 0.3 or alpha > 0.7
   sig_nonadditive = merged[
       (merged['pval'] < 5e-8) &
       ((merged['alpha_value'] < 0.3) | (merged['alpha_value'] > 0.7))
   ].sort_values('pval')
   
   print(f"\nGenome-wide significant nonadditive variants: {len(sig_nonadditive)}")
   
   if len(sig_nonadditive) > 0:
       print("\nTop significant nonadditive effects:")
       print(sig_nonadditive[[
           'variant_id', 'chr', 'pos', 'pval', 'coef',
           'alpha_value', 'inheritance_detailed'
       ]].head(10))
   
   # Create comprehensive visualization
   fig = plt.figure(figsize=(16, 12))
   gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
   
   # 1. Alpha distribution histogram
   ax1 = fig.add_subplot(gs[0, :2])
   ax1.hist(alpha_df['alpha_value'], bins=100, edgecolor='black', alpha=0.7)
   ax1.axvline(0.5, color='red', linestyle='--', linewidth=2, label='Additive (α=0.5)')
   ax1.axvline(0.3, color='orange', linestyle='--', label='Recessive threshold')
   ax1.axvline(0.7, color='orange', linestyle='--', label='Dominant threshold')
   ax1.set_xlabel('Alpha Value', fontsize=12)
   ax1.set_ylabel('Number of Variants', fontsize=12)
   ax1.set_title('Distribution of Alpha Values', fontsize=14, fontweight='bold')
   ax1.legend()
   ax1.grid(True, alpha=0.3, axis='y')
   
   # 2. Inheritance model pie chart
   ax2 = fig.add_subplot(gs[0, 2])
   inheritance_simple = alpha_df['inheritance_detailed'].apply(
       lambda x: 'Recessive' if 'recessive' in x.lower()
       else 'Dominant' if 'dominant' in x.lower()
       else 'Additive'
   )
   counts = inheritance_simple.value_counts()
   colors = ['#ff9999', '#66b3ff', '#99ff99']
   ax2.pie(counts.values, labels=counts.index, autopct='%1.1f%%',
           colors=colors, startangle=90)
   ax2.set_title('Inheritance Models', fontsize=14, fontweight='bold')
   
   # 3. Alpha vs EAF scatter
   ax3 = fig.add_subplot(gs[1, :2])
   scatter = ax3.scatter(alpha_df['eaf'], alpha_df['alpha_value'],
                         alpha=0.3, s=10, c=alpha_df['alpha_value'],
                         cmap='RdYlBu_r', vmin=0, vmax=1)
   ax3.axhline(0.5, color='red', linestyle='--', linewidth=2)
   ax3.set_xlabel('Effect Allele Frequency', fontsize=12)
   ax3.set_ylabel('Alpha Value', fontsize=12)
   ax3.set_title('Alpha Value vs Effect Allele Frequency', fontsize=14, fontweight='bold')
   ax3.grid(True, alpha=0.3)
   plt.colorbar(scatter, ax=ax3, label='Alpha Value')
   
   # 4. Boxplot by EAF bins
   ax4 = fig.add_subplot(gs[1, 2])
   alpha_df_clean = alpha_df.dropna(subset=['eaf_bin'])
   data_to_plot = [alpha_df_clean[alpha_df_clean['eaf_bin'] == bin]['alpha_value'].values
                    for bin in eaf_labels if bin in alpha_df_clean['eaf_bin'].values]
   labels_to_use = [bin for bin in eaf_labels if bin in alpha_df_clean['eaf_bin'].values]
   ax4.boxplot(data_to_plot, labels=labels_to_use)
   ax4.axhline(0.5, color='red', linestyle='--', linewidth=2)
   ax4.set_xlabel('EAF Bin', fontsize=12)
   ax4.set_ylabel('Alpha Value', fontsize=12)
   ax4.set_title('Alpha by EAF', fontsize=14, fontweight='bold')
   ax4.grid(True, alpha=0.3, axis='y')
   
   # 5. Coefficient comparison
   ax5 = fig.add_subplot(gs[2, 0])
   ax5.scatter(alpha_df['coef_hom'], alpha_df['coef_het'], alpha=0.3, s=10)
   ax5.plot([alpha_df['coef_hom'].min(), alpha_df['coef_hom'].max()],
            [alpha_df['coef_hom'].min() * 0.5, alpha_df['coef_hom'].max() * 0.5],
            'r--', label='α=0.5 (additive)')
   ax5.set_xlabel('Homozygous Coefficient', fontsize=12)
   ax5.set_ylabel('Heterozygous Coefficient', fontsize=12)
   ax5.set_title('Het vs Hom Coefficients', fontsize=14, fontweight='bold')
   ax5.legend()
   ax5.grid(True, alpha=0.3)
   
   # 6. P-value vs Alpha for significant hits
   ax6 = fig.add_subplot(gs[2, 1])
   sig_variants = merged[merged['pval'] < 1e-5]
   scatter = ax6.scatter(sig_variants['alpha_value'],
                         -np.log10(sig_variants['pval']),
                         alpha=0.5, s=20, c=sig_variants['alpha_value'],
                         cmap='RdYlBu_r', vmin=0, vmax=1)
   ax6.axvline(0.5, color='red', linestyle='--', linewidth=2)
   ax6.axhline(-np.log10(5e-8), color='green', linestyle='--',
               label='Genome-wide sig.')
   ax6.set_xlabel('Alpha Value', fontsize=12)
   ax6.set_ylabel('-log10(p-value)', fontsize=12)
   ax6.set_title('Significance vs Alpha', fontsize=14, fontweight='bold')
   ax6.legend()
   ax6.grid(True, alpha=0.3)
   
   # 7. QQ plot colored by inheritance
   ax7 = fig.add_subplot(gs[2, 2])
   merged_sorted = merged.sort_values('pval')
   expected = -np.log10(np.arange(1, len(merged_sorted) + 1) / (len(merged_sorted) + 1))
   observed = -np.log10(merged_sorted['pval'])
   
   # Color by inheritance
   colors_dict = {'Recessive': 'blue', 'Additive': 'green', 'Dominant': 'red'}
   inheritance_simple_merged = merged_sorted['alpha_value'].apply(
       lambda x: 'Recessive' if x < 0.4 else 'Dominant' if x > 0.6 else 'Additive'
   )
   
   for inh_type, color in colors_dict.items():
       mask = inheritance_simple_merged == inh_type
       ax7.scatter(expected[mask], observed[mask], alpha=0.5, s=10,
                  c=color, label=inh_type)
   
   ax7.plot([0, expected.max()], [0, expected.max()], 'k--', label='y=x')
   ax7.set_xlabel('Expected -log10(p)', fontsize=12)
   ax7.set_ylabel('Observed -log10(p)', fontsize=12)
   ax7.set_title('QQ Plot by Inheritance', fontsize=14, fontweight='bold')
   ax7.legend()
   ax7.grid(True, alpha=0.3)
   
   plt.suptitle('EDGE GWAS - Comprehensive Alpha Analysis',
                fontsize=16, fontweight='bold', y=0.995)
   plt.savefig('alpha_comprehensive_analysis.png', dpi=300, bbox_inches='tight')
   plt.close()
   
   print("\n" + "="*70)
   print("Comprehensive visualization saved to: alpha_comprehensive_analysis.png")
   print("="*70)
   
   # Save detailed report
   report_file = 'alpha_interpretation_report.txt'
   with open(report_file, 'w') as f:
       f.write("="*70 + "\n")
       f.write("EDGE GWAS - ALPHA VALUE INTERPRETATION REPORT\n")
       f.write("="*70 + "\n\n")
       
       f.write("1. OVERALL STATISTICS\n")
       f.write("-" * 70 + "\n")
       f.write(f"Total variants: {len(alpha_df):,}\n")
       f.write(f"Mean alpha: {alpha_df['alpha_value'].mean():.3f}\n")
       f.write(f"Median alpha: {alpha_df['alpha_value'].median():.3f}\n")
       f.write(f"Std alpha: {alpha_df['alpha_value'].std():.3f}\n")
       f.write(f"Range: [{alpha_df['alpha_value'].min():.3f}, {alpha_df['alpha_value'].max():.3f}]\n\n")
       
       f.write("2. INHERITANCE MODEL DISTRIBUTION\n")
       f.write("-" * 70 + "\n")
       for model, count in inheritance_counts.items():
           pct = (count / len(alpha_df)) * 100
           f.write(f"{model:20s}: {count:6,} ({pct:5.1f}%)\n")
       f.write("\n")
       
       f.write("3. INTERPRETATION GUIDE\n")
       f.write("-" * 70 + "\n")
       f.write("Alpha Range    | Inheritance    | Interpretation\n")
       f.write("-" * 70 + "\n")
       f.write("< 0            | Under-recessive| Heterozygote disadvantage\n")
       f.write("0.0 - 0.2      | Strong recessive| Strong recessive effect\n")
       f.write("0.2 - 0.4      | Recessive      | Recessive inheritance\n")
       f.write("0.4 - 0.6      | Additive       | Standard additive model\n")
       f.write("0.6 - 0.8      | Dominant       | Dominant inheritance\n")
       f.write("0.8 - 1.2      | Strong dominant| Strong dominant effect\n")
       f.write("> 1.2          | Over-dominant  | Heterozygote advantage\n\n")
       
       f.write("4. SIGNIFICANT NONADDITIVE EFFECTS\n")
       f.write("-" * 70 + "\n")
       f.write(f"Total genome-wide significant: {len(sig_nonadditive)}\n\n")
       
       if len(sig_nonadditive) > 0:
           f.write("Top 10 significant nonadditive variants:\n")
           f.write(sig_nonadditive[[
               'variant_id', 'chr', 'pos', 'pval', 'alpha_value', 'inheritance_detailed'
           ]].head(10).to_string())
   
   print(f"\nDetailed report saved to: {report_file}")

Tips and Best Practices
------------------------

Data Preparation
~~~~~~~~~~~~~~~~

1. **Sample ID Consistency** (NEW in v0.1.1):

   .. code-block:: python
   
      # Ensure all IDs are strings for consistent matching
      genotype_df.index = genotype_df.index.astype(str)
      phenotype_df.index = phenotype_df.index.astype(str)
      
      # Check for common samples
      common = set(genotype_df.index) & set(phenotype_df.index)
      print(f"Common samples: {len(common)}")

2. **Population Structure Control** (NEW in v0.1.1):

   .. code-block:: python
   
      # For unrelated samples: use standard PCA
      pca_df = calculate_pca_plink('genotypes', n_pcs=10)
      
      # For related samples: use PC-AiR
      pca_df = calculate_pca_pcair('genotypes', n_pcs=10)
      
      # For biobank data: use GRM + PC-AiR
      grm_prefix = calculate_grm_gcta('genotypes')
      pca_df = calculate_pca_pcair('genotypes', kinship_matrix=grm_prefix)

3. **Outcome Transformation** (NEW in v0.1.1):

   .. code-block:: python
   
      # Check distribution
      import matplotlib.pyplot as plt
      plt.hist(phenotype_df['trait'], bins=50)
      plt.title('Raw Trait Distribution')
      plt.show()
      
      # For skewed continuous traits: use rank-based inverse normal
      edge = EDGEAnalysis(
          outcome_type='continuous',
          outcome_transform='rank_inverse_normal'
      )
      
      # For log-distributed traits (e.g., gene expression)
      edge = EDGEAnalysis(
          outcome_type='continuous',
          outcome_transform='log'
      )

Quality Control Best Practices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Comprehensive QC Pipeline**:

   .. code-block:: python
   
      # Recommended order of QC steps
      # 1. Sample call rate (>95%)
      geno, pheno = filter_samples_by_call_rate(geno, pheno, min_call_rate=0.95)
      
      # 2. Variant missingness (<5%)
      geno = filter_variants_by_missing(geno, max_missing=0.05)
      
      # 3. MAF filtering (>1% for common variants)
      geno = filter_variants_by_maf(geno, min_maf=0.01)
      
      # 4. HWE test (p > 1e-6 in controls)
      geno = filter_variants_by_hwe(geno, hwe_threshold=1e-6)
      
      # 5. Check for relatedness
      grm_matrix, grm_ids = load_grm_gcta('grm_prefix')
      related = identify_related_samples(grm_matrix, grm_ids, threshold=0.0884)

2. **Case-Control Balance Check**:

   .. code-block:: python
   
      balance = check_case_control_balance(pheno_df, 'disease')
      
      # Warn if severely imbalanced
      if balance['ratio'] < 0.1 or balance['ratio'] > 10:
          print(f"WARNING: Imbalanced ratio {balance['ratio']:.2f}")
          print("Consider using SAIGE or other methods for imbalanced data")

Analysis Strategy
~~~~~~~~~~~~~~~~~

1. **Train/Test Split**:

   - Use 50/50 split for balanced power
   - Always stratify binary outcomes
   - Use same random seed for reproducibility

2. **Covariate Selection**:

   .. code-block:: python
   
      # Recommended covariates
      basic_covariates = ['age', 'sex', 'batch']
      
      # Add PCs (10 recommended, more for diverse populations)
      pc_covariates = get_pc_covariate_list(10)
      
      # Combine
      all_covariates = basic_covariates + pc_covariates

3. **Multiple Testing**:

   - Genome-wide significance: p < 5e-8
   - Suggestive threshold: p < 1e-5
   - Check genomic inflation (λ should be 0.95-1.05)

4. **GRM Usage Decision Tree** (NEW in v0.1.1):

   .. code-block:: python
   
      # Decision tree for GRM usage
      if has_related_samples:
          # Use GRM in analysis
          alpha_df, gwas_df = edge.run_full_analysis(
              ..., grm_matrix=grm_matrix, grm_sample_ids=grm_ids
          )
      else:
          # Use PCs only (faster)
          alpha_df, gwas_df = edge.run_full_analysis(
              ..., covariates=covariates_with_pcs
          )

Performance Optimization
~~~~~~~~~~~~~~~~~~~~~~~~

1. **Use All CPU Cores**:

   .. code-block:: python
   
      edge = EDGEAnalysis(outcome_type='binary', n_jobs=-1)

2. **Memory Management for Large Datasets**:

   .. code-block:: python
   
      # Process by chromosome
      for chrom in range(1, 23):
          # Load, analyze, save, then delete
          geno = load_chromosome(chrom)
          results = analyze(geno)
          save_results(results, f'chr{chrom}')
          del geno
          gc.collect()

3. **Use Approximate PCA for Large Cohorts** (NEW in v0.1.1):

   .. code-block:: python
   
      # For >5000 samples
      pca_df = calculate_pca_plink(
          'genotypes',
          n_pcs=10,
          approx=True,
          approx_samples=5000
      )

4. **Sparse GRM for Large Datasets** (NEW in v0.1.1):

   .. code-block:: python
   
      # Use sparse GRM for >10,000 samples
      grm_prefix = calculate_grm_gcta(
          'genotypes',
          method='grm-sparse'
      )

Interpretation Guidelines
~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Alpha Value Interpretation**:

   - α ≈ 0: Recessive model (only homozygotes affected)
   - α ≈ 0.5: Additive model (heterozygotes intermediate)
   - α ≈ 1: Dominant model (heterozygotes same as homozygotes)
   - α < 0 or α > 1: Over-dominance or under-recessiveness

2. **When EDGE Shows Advantage**:

   - Look for α far from 0.5 (non-additive effects)
   - Check if effect sizes differ between het and hom
   - Validate in independent cohorts
   - Consider biological mechanisms

3. **Follow-Up Analysis**:

   .. code-block:: python
   
      # For significant nonadditive hits
      sig_nonadditive = gwas_df[
          (gwas_df['pval'] < 5e-8) &
          ((gwas_df['alpha_value'] < 0.3) | (gwas_df['alpha_value'] > 0.7))
      ]
      
      # Prioritize for:
      # - Functional annotation
      # - eQTL analysis
      # - Validation experiments

Common Pitfalls to Avoid
~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Not controlling for population structure** → Use PCA + GRM
2. **Ignoring relatedness** → Use PC-AiR or GRM
3. **Not transforming skewed outcomes** → Use rank-based inverse normal
4. **Using different sample IDs formats** → Convert all to strings
5. **Forgetting to filter rare variants** → Use MAF > 0.01
6. **Not checking genomic inflation** → Always calculate λ

Troubleshooting Common Issues
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **"No common samples found"**:

   .. code-block:: python
   
      # Check ID formats
      print(geno_df.index[:5])
      print(pheno_df.index[:5])
      print(pca_df.index[:5])
      
      # Convert all to strings
      geno_df.index = geno_df.index.astype(str)
      pheno_df.index = pheno_df.index.astype(str)

2. **High genomic inflation (λ > 1.05)**:

   - Check for population stratification
   - Add more PCs
   - Consider using GRM

3. **Memory errors**:

   - Process chromosomes separately
   - Use sparse GRM
   - Reduce number of parallel jobs

4. **PLINK2/GCTA not found**:

   .. code-block:: bash
   
      # Install tools
      edge-gwas-install-tools
      
      # Check installation
      edge-gwas-check-tools

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
