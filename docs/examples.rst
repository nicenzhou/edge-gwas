.. _examples:

Example Workflows
=================

This page provides complete, real-world examples of EDGE GWAS analysis.

Example 1: Basic Binary Outcome Analysis
-----------------------------------------

Complete workflow for case-control study.

.. code-block:: python

   import pandas as pd
   from edge_gwas import EDGEAnalysis
   from edge_gwas.utils import (
       load_plink_data,
       prepare_phenotype_data,
       stratified_train_test_split,
       filter_variants_by_maf,
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
       covariate_cols=['age', 'sex'] + [f'PC{i}' for i in range(1, 11)],
       sample_id_col='IID'
   )
   
   # 2. Quality control
   print("Applying QC filters...")
   genotype_df = filter_variants_by_maf(genotype_df, min_maf=0.01)
   
   # 3. Split data
   print("Splitting into train/test sets...")
   train_g, test_g, train_p, test_p = stratified_train_test_split(
       genotype_df,
       phenotype_df,
       outcome_col='disease_status',
       test_size=0.5,
       random_state=42
   )
   
   # 4. Run EDGE analysis
   print("Running EDGE analysis...")
   edge = EDGEAnalysis(outcome_type='binary', n_jobs=8, verbose=True)
   alpha_df, gwas_df = edge.run_full_analysis(
       train_genotype=train_g,
       train_phenotype=train_p,
       test_genotype=test_g,
       test_phenotype=test_p,
       outcome='disease_status',
       covariates=['age', 'sex'] + [f'PC{i}' for i in range(1, 11)],
       variant_info=variant_info,
       output_prefix='chr1_edge_results'
   )
   
   # 5. Check results
   print(f"\nTotal variants tested: {len(gwas_df)}")
   print(f"Significant hits (p < 5e-8): {(gwas_df['pval'] < 5e-8).sum()}")
   
   lambda_gc = calculate_genomic_inflation(gwas_df['pval'])
   print(f"Genomic inflation factor: {lambda_gc:.3f}")
   
   # 6. Visualize
   manhattan_plot(gwas_df, 'chr1_manhattan.png', title='Chromosome 1 EDGE GWAS')
   qq_plot(gwas_df, 'chr1_qq.png')
   
   # 7. Extract top hits
   top_hits = gwas_df[gwas_df['pval'] < 5e-8].sort_values('pval')
   print(f"\nTop hits:\n{top_hits[['variant_id', 'chr', 'pos', 'pval', 'coef']]}")

Example 2: Quantitative Trait Analysis
---------------------------------------

Analysis of continuous phenotype.

.. code-block:: python

   from edge_gwas import EDGEAnalysis
   from edge_gwas.utils import load_plink_data, prepare_phenotype_data
   
   # Load data
   genotype_df, variant_info = load_plink_data('data.bed', 'data.bim', 'data.fam')
   
   # Prepare quantitative phenotype (e.g., BMI)
   phenotype_df = prepare_phenotype_data(
       phenotype_file='phenotypes.txt',
       outcome_col='BMI',
       covariate_cols=['age', 'sex', 'PC1', 'PC2', 'PC3'],
       log_transform_outcome=False  # Set True if needed
   )
   
   # Split data
   from edge_gwas.utils import stratified_train_test_split
   train_g, test_g, train_p, test_p = stratified_train_test_split(
       genotype_df, phenotype_df, 'BMI', test_size=0.5, is_binary=False
   )
   
   # Run analysis with continuous outcome
   edge = EDGEAnalysis(outcome_type='continuous', n_jobs=-1)
   alpha_df, gwas_df = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       outcome='BMI',
       covariates=['age', 'sex', 'PC1', 'PC2', 'PC3']
   )
   
   # Visualize
   from edge_gwas.visualize import manhattan_plot, qq_plot, plot_alpha_distribution
   manhattan_plot(gwas_df, 'bmi_manhattan.png')
   qq_plot(gwas_df, 'bmi_qq.png')
   plot_alpha_distribution(alpha_df, 'bmi_alpha_dist.png')

Example 3: Multi-Chromosome Analysis
-------------------------------------

Analyze multiple chromosomes and combine results.

.. code-block:: python

   import pandas as pd
   from edge_gwas import EDGEAnalysis
   from edge_gwas.utils import load_plink_data, prepare_phenotype_data
   
   # Load phenotype once
   phenotype_df = prepare_phenotype_data(
       'phenotypes.txt', 'disease', 
       ['age', 'sex'] + [f'PC{i}' for i in range(1, 11)]
   )
   
   # Initialize analysis
   edge = EDGEAnalysis(outcome_type='binary', n_jobs=8)
   
   all_alpha = []
   all_gwas = []
   
   # Process each chromosome
   for chrom in range(1, 23):  # Chromosomes 1-22
       print(f"\nProcessing chromosome {chrom}...")
       
       # Load chromosome data
       genotype_df, variant_info = load_plink_data(
           bed_file=f'data/chr{chrom}.bed',
           bim_file=f'data/chr{chrom}.bim',
           fam_file=f'data/chr{chrom}.fam'
       )
       
       # Split data
       from edge_gwas.utils import stratified_train_test_split
       train_g, test_g, train_p, test_p = stratified_train_test_split(
           genotype_df, phenotype_df, 'disease', test_size=0.5, random_state=42
       )
       
       # Run EDGE
       alpha_df, gwas_df = edge.run_full_analysis(
           train_g, train_p, test_g, test_p,
           outcome='disease',
           covariates=['age', 'sex'] + [f'PC{i}' for i in range(1, 11)],
           variant_info=variant_info,
           output_prefix=f'chr{chrom}_edge'
       )
       
       # Store results
       all_alpha.append(alpha_df)
       all_gwas.append(gwas_df)
   
   # Combine results
   combined_alpha = pd.concat(all_alpha, ignore_index=True)
   combined_gwas = pd.concat(all_gwas, ignore_index=True)
   
   # Save combined results
   combined_alpha.to_csv('all_chromosomes_alpha.txt', sep='\t', index=False)
   combined_gwas.to_csv('all_chromosomes_gwas.txt', sep='\t', index=False)
   
   # Visualize genome-wide results
   from edge_gwas.visualize import manhattan_plot, qq_plot
   manhattan_plot(combined_gwas, 'genome_wide_manhattan.png', 
                  title='Genome-Wide EDGE GWAS')
   qq_plot(combined_gwas, 'genome_wide_qq.png')
   
   # Summary
   sig_hits = combined_gwas[combined_gwas['pval'] < 5e-8]
   print(f"\nTotal genome-wide significant hits: {len(sig_hits)}")
   print(f"Chromosomes with hits: {sig_hits['chr'].unique()}")

Example 4: Using Pre-Calculated Alpha Values
---------------------------------------------

Apply previously calculated alpha values to new test data.

.. code-block:: python

   import pandas as pd
   from edge_gwas import EDGEAnalysis
   from edge_gwas.io_handlers import load_alpha_values
   from edge_gwas.utils import load_plink_data, prepare_phenotype_data
   
   # Load pre-calculated alpha values
   alpha_df = load_alpha_values('previous_alpha_values.txt')
   print(f"Loaded {len(alpha_df)} alpha values")
   
   # Load new test data
   test_genotype, variant_info = load_plink_data(
       'new_test_data.bed', 'new_test_data.bim', 'new_test_data.fam'
   )
   
   test_phenotype = prepare_phenotype_data(
       'new_phenotypes.txt', 'disease',
       ['age', 'sex', 'PC1', 'PC2', 'PC3']
   )
   
   # Apply alpha values
   edge = EDGEAnalysis(outcome_type='binary')
   gwas_df = edge.apply_alpha(
       genotype_data=test_genotype,
       phenotype_df=test_phenotype,
       outcome='disease',
       covariates=['age', 'sex', 'PC1', 'PC2', 'PC3'],
       alpha_values=alpha_df
   )
   
   # Save and visualize
   gwas_df.to_csv('new_cohort_gwas.txt', sep='\t', index=False)
   
   from edge_gwas.visualize import manhattan_plot, qq_plot
   manhattan_plot(gwas_df, 'new_cohort_manhattan.png')
   lambda_gc = qq_plot(gwas_df, 'new_cohort_qq.png')
   print(f"Genomic inflation: {lambda_gc:.3f}")

Example 5: Custom Quality Control Pipeline
-------------------------------------------

Detailed QC with multiple filters including new HWE and sample filtering.

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
       stratified_train_test_split
   )
   
   # Load data
   genotype_df, variant_info = load_plink_data('data.bed', 'data.bim', 'data.fam')
   phenotype_df = prepare_phenotype_data(
       'phenotypes.txt', 'disease', ['age', 'sex'] + [f'PC{i}' for i in range(1, 11)]
   )
   
   print(f"Starting variants: {genotype_df.shape[1]}")
   print(f"Starting samples: {genotype_df.shape[0]}")
   
   # QC Step 1: Filter by MAF
   genotype_df = filter_variants_by_maf(
       genotype_df, min_maf=0.01, verbose=True
   )
   print(f"After MAF filter: {genotype_df.shape[1]} variants")
   
   # QC Step 2: Filter by missingness
   genotype_df = filter_variants_by_missing(
       genotype_df, max_missing=0.05, verbose=True
   )
   print(f"After missingness filter: {genotype_df.shape[1]} variants")
   
   # QC Step 3: Filter by Hardy-Weinberg Equilibrium (NEW)
   genotype_df = filter_variants_by_hwe(
       genotype_df, hwe_threshold=1e-6, verbose=True
   )
   print(f"After HWE filter: {genotype_df.shape[1]} variants")
   
   # QC Step 4: Filter samples by call rate (NEW)
   genotype_df, phenotype_df = filter_samples_by_call_rate(
       genotype_df, phenotype_df, min_call_rate=0.95, verbose=True
   )
   print(f"After sample QC: {genotype_df.shape[0]} samples")
   
   # QC Step 5: Check case/control balance (NEW)
   balance = check_case_control_balance(
       phenotype_df, outcome_col='disease', verbose=True
   )
   
   # Continue with analysis
   train_g, test_g, train_p, test_p = stratified_train_test_split(
       genotype_df, phenotype_df, 'disease', test_size=0.5, random_state=42
   )
   
   edge = EDGEAnalysis(outcome_type='binary', n_jobs=8)
   alpha_df, gwas_df = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       outcome='disease',
       covariates=['age', 'sex'] + [f'PC{i}' for i in range(1, 11)]
   )

Example 6: Interpretation of Alpha Values
------------------------------------------

Understanding and interpreting alpha values.

.. code-block:: python

   import pandas as pd
   import numpy as np
   from edge_gwas import EDGEAnalysis
   from edge_gwas.visualize import plot_alpha_distribution
   
   # After running analysis
   edge = EDGEAnalysis(outcome_type='binary')
   alpha_df, gwas_df = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       outcome='disease', covariates=['age', 'sex']
   )
   
   # Categorize alpha values by inheritance pattern
   def categorize_alpha(alpha):
       """Categorize alpha by inheritance model"""
       if alpha < 0:
           return 'Under-recessive'
       elif 0 <= alpha < 0.3:
           return 'Recessive'
       elif 0.3 <= alpha < 0.7:
           return 'Additive'
       elif 0.7 <= alpha < 1.2:
           return 'Dominant'
       else:
           return 'Over-dominant'
   
   alpha_df['inheritance_model'] = alpha_df['alpha_value'].apply(categorize_alpha)
   
   # Summary statistics
   print("\nAlpha value distribution:")
   print(alpha_df['inheritance_model'].value_counts())
   print(f"\nMean alpha: {alpha_df['alpha_value'].mean():.3f}")
   print(f"Median alpha: {alpha_df['alpha_value'].median():.3f}")
   print(f"Std alpha: {alpha_df['alpha_value'].std():.3f}")
   
   # Find variants with extreme alpha values
   recessive = alpha_df[alpha_df['alpha_value'] < 0.3].sort_values('alpha_value')
   dominant = alpha_df[alpha_df['alpha_value'] > 0.7].sort_values('alpha_value', ascending=False)
   
   print(f"\nRecessive variants (alpha < 0.3): {len(recessive)}")
   print(recessive[['variant_id', 'alpha_value', 'eaf']].head(10))
   
   print(f"\nDominant variants (alpha > 0.7): {len(dominant)}")
   print(dominant[['variant_id', 'alpha_value', 'eaf']].head(10))
   
   # Merge with GWAS results to find significant nonadditive effects
   from edge_gwas.utils import merge_alpha_with_gwas
   merged = merge_alpha_with_gwas(gwas_df, alpha_df)
   
   # Significant nonadditive effects
   sig_nonadditive = merged[
       (merged['pval'] < 5e-8) & 
       ((merged['alpha_value'] < 0.3) | (merged['alpha_value'] > 0.7))
   ]
   
   print(f"\nSignificant nonadditive effects: {len(sig_nonadditive)}")
   print(sig_nonadditive[['variant_id', 'chr', 'pos', 'pval', 'alpha_value', 'inheritance_model']])
   
   # Visualize
   plot_alpha_distribution(alpha_df, 'alpha_distribution.png')

Example 7: Generate Comprehensive Report
-----------------------------------------

Create a complete analysis report.

.. code-block:: python

   from edge_gwas import EDGEAnalysis
   from edge_gwas.io_handlers import create_summary_report, save_results
   from edge_gwas.utils import calculate_genomic_inflation
   from edge_gwas.visualize import manhattan_plot, qq_plot, plot_alpha_distribution
   
   # Run analysis
   edge = EDGEAnalysis(outcome_type='binary', n_jobs=8, verbose=True)
   alpha_df, gwas_df = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       outcome='disease',
       covariates=['age', 'sex'] + [f'PC{i}' for i in range(1, 11)],
       output_prefix='final_analysis'
   )
   
   # Calculate statistics
   lambda_gc = calculate_genomic_inflation(gwas_df['pval'])
   sig_count = (gwas_df['pval'] < 5e-8).sum()
   suggestive_count = ((gwas_df['pval'] >= 5e-8) & (gwas_df['pval'] < 1e-5)).sum()
   
   # Generate visualizations
   manhattan_plot(gwas_df, 'final_manhattan.png', title='Final EDGE GWAS Results')
   qq_plot(gwas_df, 'final_qq.png', title='Final QQ Plot')
   plot_alpha_distribution(alpha_df, 'final_alpha_dist.png')
   
   # Create summary report
   report = create_summary_report(
       gwas_df=gwas_df,
       alpha_df=alpha_df,
       significance_threshold=5e-8,
       output_file='analysis_summary.txt'
   )
   
   # Save all results
   output_files = save_results(
       gwas_df=gwas_df,
       alpha_df=alpha_df,
       output_prefix='final_results',
       save_alpha=True
   )
   
   # Print summary
   print("=" * 60)
   print("EDGE GWAS ANALYSIS SUMMARY")
   print("=" * 60)
   print(f"\nTotal variants tested: {len(gwas_df):,}")
   print(f"Genome-wide significant (p < 5e-8): {sig_count}")
   print(f"Suggestive (p < 1e-5): {suggestive_count}")
   print(f"Genomic inflation factor (λ): {lambda_gc:.3f}")
   print(f"\nMean alpha: {alpha_df['alpha_value'].mean():.3f}")
   print(f"Median alpha: {alpha_df['alpha_value'].median():.3f}")
   print(f"\nOutput files:")
   for key, filepath in output_files.items():
       print(f"  {key}: {filepath}")
   print(f"  Manhattan plot: final_manhattan.png")
   print(f"  QQ plot: final_qq.png")
   print(f"  Alpha distribution: final_alpha_dist.png")
   print(f"  Summary report: analysis_summary.txt")
   print("=" * 60)

Example 8: Handling Large Datasets
-----------------------------------

Memory-efficient analysis of large genomic datasets.

.. code-block:: python

   import pandas as pd
   from edge_gwas import EDGEAnalysis
   from edge_gwas.utils import load_plink_data
   
   # Process chromosomes in batches to save memory
   def process_chromosome_batch(chromosomes, phenotype_df, output_prefix):
       """Process multiple chromosomes efficiently"""
       edge = EDGEAnalysis(outcome_type='binary', n_jobs=-1, verbose=True)
       
       results = []
       for chrom in chromosomes:
           print(f"\nProcessing chromosome {chrom}...")
           
           # Load one chromosome at a time
           geno, info = load_plink_data(
               f'data/chr{chrom}.bed',
               f'data/chr{chrom}.bim',
               f'data/chr{chrom}.fam'
           )
           
           # Split
           from edge_gwas.utils import stratified_train_test_split
           train_g, test_g, train_p, test_p = stratified_train_test_split(
               geno, phenotype_df, 'disease', test_size=0.5, random_state=42
           )
           
           # Analyze
           alpha_df, gwas_df = edge.run_full_analysis(
               train_g, train_p, test_g, test_p,
               outcome='disease',
               covariates=['age', 'sex'] + [f'PC{i}' for i in range(1, 11)],
               output_prefix=f'{output_prefix}_chr{chrom}'
           )
           
           results.append((alpha_df, gwas_df))
           
           # Clear memory
           del geno, info, train_g, test_g
           
       return results
   
   # Load phenotype once
   from edge_gwas.utils import prepare_phenotype_data
   phenotype_df = prepare_phenotype_data(
       'phenotypes.txt', 'disease',
       ['age', 'sex'] + [f'PC{i}' for i in range(1, 11)]
   )
   
   # Process in batches
   batch1 = process_chromosome_batch([1, 2, 3, 4, 5], phenotype_df, 'batch1')
   batch2 = process_chromosome_batch([6, 7, 8, 9, 10], phenotype_df, 'batch2')
   batch3 = process_chromosome_batch([11, 12, 13, 14, 15], phenotype_df, 'batch3')
   batch4 = process_chromosome_batch([16, 17, 18, 19, 20, 21, 22], phenotype_df, 'batch4')
   
   # Combine all results
   all_results = batch1 + batch2 + batch3 + batch4
   all_alpha = pd.concat([r[0] for r in all_results], ignore_index=True)
   all_gwas = pd.concat([r[1] for r in all_results], ignore_index=True)
   
   # Save combined results
   all_alpha.to_csv('genome_wide_alpha.txt', sep='\t', index=False)
   all_gwas.to_csv('genome_wide_gwas.txt', sep='\t', index=False)
   
   print(f"\nTotal variants analyzed: {len(all_gwas):,}")

Example 9: Comparing EDGE with Standard Additive GWAS
------------------------------------------------------

Compare EDGE results with traditional additive GWAS using built-in function.

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
       covariates=['age', 'sex', 'PC1', 'PC2']
   )
   
   # Run standard additive GWAS using built-in function (NEW)
   additive_results = additive_gwas(
       genotype_df=test_g,
       phenotype_df=test_p,
       outcome='disease',
       covariates=['age', 'sex', 'PC1', 'PC2'],
       outcome_type='binary'
   )
   
   # Merge results for comparison
   comparison = edge_gwas_results[['variant_id', 'pval', 'coef']].merge(
       additive_results[['variant_id', 'pval', 'coef']],
       on='variant_id',
       suffixes=('_edge', '_additive')
   )
   
   # Find variants where EDGE performs better
   comparison['-log10p_edge'] = -np.log10(comparison['pval_edge'])
   comparison['-log10p_additive'] = -np.log10(comparison['pval_additive'])
   comparison['edge_advantage'] = comparison['-log10p_edge'] - comparison['-log10p_additive']
   
   # Variants where EDGE is more significant
   edge_better = comparison[comparison['edge_advantage'] > 2].sort_values(
       'edge_advantage', ascending=False
   )
   
   print(f"\nVariants where EDGE shows advantage (>100-fold p-value improvement):")
   print(edge_better[['variant_id', 'pval_edge', 'pval_additive', 'edge_advantage']].head(10))
   
   # Plot comparison
   fig, axes = plt.subplots(1, 2, figsize=(14, 6))
   
   # Scatter plot
   axes[0].scatter(comparison['-log10p_additive'], comparison['-log10p_edge'], 
                   alpha=0.5, s=10)
   axes[0].plot([0, comparison['-log10p_additive'].max()], 
                [0, comparison['-log10p_additive'].max()], 
                'r--', label='y=x')
   axes[0].set_xlabel('Additive -log10(p)')
   axes[0].set_ylabel('EDGE -log10(p)')
   axes[0].set_title('EDGE vs Additive GWAS')
   axes[0].legend()
   axes[0].grid(True, alpha=0.3)
   
   # Advantage histogram
   axes[1].hist(comparison['edge_advantage'], bins=50, edgecolor='black', alpha=0.7)
   axes[1].axvline(0, color='red', linestyle='--', label='No difference')
   axes[1].set_xlabel('EDGE Advantage (-log10p difference)')
   axes[1].set_ylabel('Number of variants')
   axes[1].set_title('Distribution of EDGE Advantage')
   axes[1].legend()
   axes[1].grid(True, alpha=0.3, axis='y')
   
   plt.tight_layout()
   plt.savefig('edge_vs_additive_comparison.png', dpi=300)
   plt.close()

Example 10: Cross-Validation Strategy
--------------------------------------

Implement k-fold cross-validation using built-in function.

.. code-block:: python

   import pandas as pd
   from edge_gwas.utils import cross_validated_edge_analysis
   
   # Run cross-validated analysis using built-in function (NEW)
   avg_alpha, meta_gwas, all_alpha, all_gwas = cross_validated_edge_analysis(
       genotype_df=genotype_data,
       phenotype_df=phenotype_df,
       outcome='disease',
       covariates=['age', 'sex'] + [f'PC{i}' for i in range(1, 11)],
       outcome_type='binary',  # NEW: Support for continuous outcomes
       n_folds=5,
       n_jobs=8,
       random_state=42
   )
   
   # Check alpha stability across folds
   print("\nAlpha value stability:")
   print(f"Mean alpha std across variants: {avg_alpha['alpha_std'].mean():.3f}")
   
   unstable_alpha = avg_alpha[avg_alpha['alpha_std'] > 0.3]
   print(f"Variants with unstable alpha (std > 0.3): {len(unstable_alpha)}")
   
   # Save results
   avg_alpha.to_csv('cv_averaged_alpha.txt', sep='\t', index=False)
   meta_gwas.to_csv('cv_meta_gwas.txt', sep='\t', index=False)
   
   # Visualize
   from edge_gwas.visualize import manhattan_plot, qq_plot
   manhattan_plot(meta_gwas, 'cv_manhattan.png', title='Cross-Validated EDGE GWAS')
   qq_plot(meta_gwas, 'cv_qq.png')

Tips and Best Practices
------------------------

Data Preparation
~~~~~~~~~~~~~~~~

1. **Always include principal components** (at least 10 PCs) to control for population structure
2. **Check phenotype distribution** before analysis
3. **Standardize covariates** if scales differ greatly
4. **Verify sample IDs match** between genotype and phenotype files

Quality Control
~~~~~~~~~~~~~~~

1. **MAF filtering**: Use MAF > 0.01 for binary traits, can go lower for quantitative traits
2. **Missingness**: Remove variants with >5% missing data
3. **HWE**: Test in controls, use p > 1e-6 threshold
4. **Sample QC**: Remove samples with call rate < 95%

Analysis Strategy
~~~~~~~~~~~~~~~~~

1. **Train/test split**: Use 50/50 split for balanced power
2. **Stratification**: Always stratify by outcome for binary traits
3. **Multiple testing**: Use genome-wide significance threshold (p < 5e-8)
4. **Genomic control**: Check λ, should be between 0.95-1.05

Interpretation
~~~~~~~~~~~~~~

1. **Alpha values**: Check distribution, look for extreme values
2. **Compare with additive**: Identify variants where EDGE shows advantage
3. **Biological validation**: Follow up significant hits with functional studies
4. **Replication**: Validate findings in independent cohorts

Performance Optimization
~~~~~~~~~~~~~~~~~~~~~~~~

1. **Use all CPU cores**: Set ``n_jobs=-1`` for parallel processing
2. **Process by chromosome**: Analyze chromosomes separately for large datasets
3. **Filter early**: Apply QC filters before splitting data
4. **Monitor memory**: Use batch processing for very large datasets

See Also
--------

* :ref:`api_reference` - Complete API documentation
* :ref:`statistical_model` - Statistical methodology
* :ref:`visualization` - Visualization options
* :ref:`user_guide` - Detailed user guide
