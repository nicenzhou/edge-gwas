.. _user_guide:

User Guide
==========

This comprehensive guide covers all aspects of using edge-gwas v0.1.1.

Overview
--------

edge-gwas is designed for researchers performing large-scale GWAS analysis with flexible genetic encoding to identify nonadditive SNP effects.

**Key Features:**

* Two-stage analysis framework (training/test split)
* Flexible alpha encoding for nonadditive effects detection
* Support for multiple genotype file formats
* Comprehensive quality control tools
* Built-in visualization functions
* Parallel processing for large datasets

Core Concepts
-------------

EDGE Methodology
~~~~~~~~~~~~~~~~

**EDGE (Elastic Data-Driven Encoding)** identifies genetic variants with nonadditive effects by:

1. **Training Stage**: Calculate flexible encoding parameters (alpha) from training data
2. **Test Stage**: Apply alpha values to test data for association testing

This approach detects:

* Recessive effects (α ≈ 0)
* Additive effects (α ≈ 0.5)
* Dominant effects (α ≈ 1)
* Over-dominant effects (α > 1)
* Under-recessive effects (α < 0)

See :ref:`statistical_model` for detailed mathematical framework.

Two-Stage Analysis
~~~~~~~~~~~~~~~~~~

**Why split data?**

The two-stage approach prevents overfitting and provides unbiased estimates:

1. **Training set**: Estimate alpha (encoding parameter) for each variant
2. **Test set**: Perform GWAS using pre-calculated alpha values

Recommended split: **50/50** for balanced power and bias control.

Data Formats
------------

Supported Input Formats
~~~~~~~~~~~~~~~~~~~~~~~

edge-gwas v0.1.1 supports all major genotype file formats:

**Genotype Data:**

* **PLINK (.bed/.bim/.fam)** - Standard PLINK binary format
* **PLINK 2 (.pgen/.pvar/.psam)** - New PLINK 2 format with better compression
* **BGEN (.bgen)** - UK Biobank format with genotype dosages
* **VCF (.vcf/.vcf.gz)** - Variant Call Format (compressed or uncompressed)

**Phenotype Data:**

* **TXT/TSV/CSV** - Tab/comma-delimited text files
* Required columns: Sample ID, outcome, covariates

See :ref:`examples` for format specifications and loading examples.

Output Formats
~~~~~~~~~~~~~~

**Association Results:**

* **TXT/TSV** - Tab-delimited association results
* **CSV** - Comma-delimited results for Excel

**Visualization:**

* **PNG** - High-resolution plots (300 DPI)
* **PDF** - Vector graphics for publications
* **SVG** - Scalable vector graphics

**Alpha Values:**

* **TXT/TSV** - Encoding parameters for each variant

Workflow Components
-------------------

1. Data Loading
~~~~~~~~~~~~~~~

Load and prepare genomic and phenotype data:

.. code-block:: python

   from edge_gwas.utils import load_plink_data, prepare_phenotype_data
   
   # Load genotypes (choose format)
   geno, info = load_plink_data('data.bed', 'data.bim', 'data.fam')
   # OR
   geno, info = load_pgen_data('data.pgen', 'data.pvar', 'data.psam')
   # OR
   geno, info = load_bgen_data('data.bgen', 'data.sample')
   # OR
   geno, info = load_vcf_data('data.vcf.gz', dosage=True)
   
   # Load phenotypes
   pheno = prepare_phenotype_data(
       'phenotypes.txt',
       outcome_col='disease',
       covariate_cols=['age', 'sex', 'PC1', 'PC2', 'PC3']
   )

**Best Practice:** Use binary formats (PLINK, PLINK 2, BGEN) for large datasets.

2. Quality Control
~~~~~~~~~~~~~~~~~~

Filter variants and samples based on standard QC criteria:

.. code-block:: python

   from edge_gwas.utils import (
       filter_variants_by_maf,
       filter_variants_by_missing,
       filter_variants_by_hwe,
       filter_samples_by_call_rate,
       check_case_control_balance
   )
   
   # Variant QC
   geno = filter_variants_by_maf(geno, min_maf=0.01)
   geno = filter_variants_by_missing(geno, max_missing=0.05)
   geno = filter_variants_by_hwe(geno, hwe_threshold=1e-6)
   
   # Sample QC
   geno, pheno = filter_samples_by_call_rate(geno, pheno, min_call_rate=0.95)
   
   # Check case/control balance (for binary outcomes)
   balance = check_case_control_balance(pheno, 'disease')

**Recommended Thresholds:**

* MAF > 0.01 (1%)
* Missingness < 5%
* HWE p-value > 1e-6 (in controls)
* Sample call rate > 95%

3. Data Splitting
~~~~~~~~~~~~~~~~~

Split data into training and test sets:

.. code-block:: python

   from edge_gwas.utils import stratified_train_test_split
   
   train_g, test_g, train_p, test_p = stratified_train_test_split(
       geno, pheno,
       outcome_col='disease',
       test_size=0.5,
       random_state=42,
       is_binary=True  # Stratify for balanced case/control
   )

**Important:** Use stratification for binary outcomes to maintain case/control balance.

4. EDGE Analysis
~~~~~~~~~~~~~~~~

Run the two-stage EDGE analysis:

.. code-block:: python

   from edge_gwas import EDGEAnalysis
   
   # Initialize
   edge = EDGEAnalysis(
       outcome_type='binary',  # or 'continuous'
       n_jobs=8,               # parallel processing
       verbose=True
   )
   
   # Run complete analysis
   alpha_df, gwas_df = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       outcome='disease',
       covariates=['age', 'sex', 'PC1', 'PC2', 'PC3'],
       output_prefix='my_analysis'
   )

**Alternative (Two-Step):**

.. code-block:: python

   # Step 1: Calculate alpha on training data
   alpha_df = edge.calculate_alpha(train_g, train_p, 'disease', covariates)
   
   # Step 2: Apply alpha on test data
   gwas_df = edge.apply_alpha(test_g, test_p, 'disease', covariates, alpha_df)

5. Results Interpretation
~~~~~~~~~~~~~~~~~~~~~~~~~

Interpret GWAS results and alpha values:

.. code-block:: python

   # Significant associations
   sig_hits = gwas_df[gwas_df['pval'] < 5e-8]
   print(f"Significant variants: {len(sig_hits)}")
   
   # Check alpha distribution
   print(f"Mean alpha: {alpha_df['alpha_value'].mean():.3f}")
   
   # Categorize by inheritance pattern
   recessive = alpha_df[alpha_df['alpha_value'] < 0.3]
   additive = alpha_df[(alpha_df['alpha_value'] >= 0.3) & (alpha_df['alpha_value'] < 0.7)]
   dominant = alpha_df[alpha_df['alpha_value'] >= 0.7]
   
   print(f"Recessive: {len(recessive)}")
   print(f"Additive: {len(additive)}")
   print(f"Dominant: {len(dominant)}")

**Interpretation Guide:**

* **p-value < 5×10⁻⁸**: Genome-wide significant
* **p-value < 1×10⁻⁵**: Suggestive
* **α ≈ 0**: Recessive inheritance
* **α ≈ 0.5**: Additive inheritance
* **α ≈ 1**: Dominant inheritance

6. Visualization
~~~~~~~~~~~~~~~~

Create publication-quality plots:

.. code-block:: python

   from edge_gwas.visualize import manhattan_plot, qq_plot, plot_alpha_distribution
   
   # Manhattan plot
   manhattan_plot(gwas_df, 'manhattan.png', title='EDGE GWAS Results')
   
   # QQ plot with lambda
   lambda_gc = qq_plot(gwas_df, 'qq_plot.png')
   print(f"Genomic inflation: {lambda_gc:.3f}")
   
   # Alpha distribution
   plot_alpha_distribution(alpha_df, 'alpha_dist.png')

Best Practices
--------------

Quality Control
~~~~~~~~~~~~~~~

**Pre-Analysis QC:**

1. Filter variants with MAF < 0.01
2. Remove variants with missingness > 5%
3. Test Hardy-Weinberg Equilibrium in controls (p > 1e-6)
4. Remove samples with call rate < 95%
5. Check for batch effects and population stratification

**During Analysis:**

6. Include at least 10 principal components as covariates
7. Monitor convergence warnings
8. Check for skipped SNPs due to convergence issues

**Post-Analysis QC:**

9. Calculate genomic inflation factor (λ should be 0.95-1.05)
10. Verify case/control balance in split datasets
11. Check for systematic bias in QQ plots

Statistical Analysis
~~~~~~~~~~~~~~~~~~~~

**Covariates:**

1. **Always include**: Age, sex, principal components (PC1-PC10)
2. **Consider**: Genotyping batch, study site, other technical factors
3. **Avoid**: Post-treatment variables, mediators

**Multiple Testing:**

4. Use genome-wide significance threshold: **p < 5×10⁻⁸**
5. Report suggestive threshold: **p < 1×10⁻⁵**
6. Consider Bonferroni correction for hypothesis testing
7. Use FDR for exploratory analysis

**Validation:**

8. Replicate findings in independent cohorts
9. Validate top hits with different methods (e.g., additive GWAS)
10. Compare EDGE results with standard additive models
11. Perform functional validation for biological plausibility

Performance Optimization
~~~~~~~~~~~~~~~~~~~~~~~~

**Data Handling:**

1. Use binary file formats (PLINK, PLINK 2, BGEN) for faster I/O
2. Filter data before splitting to reduce memory usage
3. Process chromosomes separately for very large datasets
4. Use chunking for extremely large cohorts

**Computational Efficiency:**

5. Utilize parallel processing: set ``n_jobs=-1`` to use all cores
6. Increase ``max_iter`` if convergence warnings appear
7. Monitor memory usage and adjust batch sizes
8. Use HPC/cloud computing for genome-wide analyses

**Memory Management:**

.. code-block:: python

   # Process chromosomes in batches
   for chrom in range(1, 23):
       geno, info = load_plink_data(f'chr{chrom}.bed', f'chr{chrom}.bim', f'chr{chrom}.fam')
       # ... process ...
       del geno, info  # Free memory

Advanced Topics
---------------

Cross-Validation
~~~~~~~~~~~~~~~~

Assess alpha stability across folds:

.. code-block:: python

   from edge_gwas.utils import cross_validated_edge_analysis
   
   avg_alpha, meta_gwas, all_alpha, all_gwas = cross_validated_edge_analysis(
       geno, pheno,
       outcome='disease',
       covariates=['age', 'sex', 'PC1', 'PC2'],
       outcome_type='binary',
       n_folds=5
   )
   
   # Check alpha stability
   unstable = avg_alpha[avg_alpha['alpha_std'] > 0.3]
   print(f"Unstable variants: {len(unstable)}")

Comparison with Additive GWAS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Compare EDGE with standard additive models:

.. code-block:: python

   from edge_gwas.utils import additive_gwas
   
   # Run standard additive GWAS
   additive_results = additive_gwas(
       test_g, test_p,
       outcome='disease',
       covariates=['age', 'sex'],
       outcome_type='binary'
   )
   
   # Compare results
   merged = gwas_df.merge(additive_results, on='variant_id', suffixes=('_edge', '_add'))
   merged['edge_advantage'] = -np.log10(merged['pval_edge']) - (-np.log10(merged['pval_add']))
   
   # Find variants where EDGE performs better
   edge_better = merged[merged['edge_advantage'] > 2]  # 100-fold improvement

Population Stratification
~~~~~~~~~~~~~~~~~~~~~~~~~

Control for population structure:

.. code-block:: python

   # Include first 10 PCs as covariates
   covariates = ['age', 'sex'] + [f'PC{i}' for i in range(1, 11)]
   
   # Check genomic inflation
   lambda_gc = calculate_genomic_inflation(gwas_df['pval'])
   
   # If λ > 1.05, consider adding more PCs or using mixed models

Meta-Analysis
~~~~~~~~~~~~~

Combine results from multiple studies:

.. code-block:: python

   from scipy.stats import combine_pvalues
   
   # Combine p-values from multiple cohorts using Fisher's method
   cohort1_pvals = gwas_cohort1['pval']
   cohort2_pvals = gwas_cohort2['pval']
   
   _, combined_pval = combine_pvalues([cohort1_pvals, cohort2_pvals], method='fisher')

Rare Variant Analysis
~~~~~~~~~~~~~~~~~~~~~

Special considerations for rare variants (MAF < 1%):

* Use lower MAF threshold (e.g., 0.001)
* Increase sample size requirements
* Consider gene-based or region-based tests
* Validate with sequencing data

Gene-Based Analysis
~~~~~~~~~~~~~~~~~~~

Aggregate signals within genes (future feature):

.. code-block:: python

   # Planned for v0.2.0
   # gene_results = edge.gene_based_test(gwas_df, gene_annotations)

Troubleshooting
---------------

Common Issues
~~~~~~~~~~~~~

**Memory errors**

* Reduce batch size or process chromosomes separately
* Filter more aggressively (increase MAF threshold)
* Use lower-memory file formats
* Close unused variables with ``del``

**Slow performance**

* Enable parallel processing: ``n_jobs=-1``
* Use binary file formats (PLINK, PLINK 2)
* Filter data before analysis
* Use HPC or cloud resources

**Convergence warnings**

* Increase ``max_iter`` parameter
* Check for multicollinearity in covariates
* Verify data quality (no constant variants)
* Use ``get_skipped_snps()`` to identify problematic variants

**Missing values**

* Check data preprocessing
* Verify sample ID matching between genotype and phenotype
* Use ``dropna()`` or imputation before analysis

**High genomic inflation (λ > 1.05)**

* Add more principal components (PC1-PC20)
* Check for population stratification
* Verify case/control balance
* Look for batch effects or technical artifacts

Error Messages
~~~~~~~~~~~~~~

**"No common samples found"**

* Check sample ID formatting (IID column in phenotype file)
* Verify sample IDs match between genotype and phenotype files
* Ensure index is set correctly on both DataFrames

**"Model fitting failed"**

* Check for variants with no variation (all same genotype)
* Verify sufficient sample size
* Check for perfect separation in binary outcomes
* Review covariate values for extreme outliers

**"Alpha value is NaN"**

* Homozygous coefficient too close to zero
* Insufficient sample size for that variant
* Convergence failure - increase ``max_iter``

Getting Help
------------

Resources
~~~~~~~~~

* **Documentation**: https://edge-gwas.readthedocs.io/
* **API Reference**: :ref:`api_reference`
* **Examples**: :ref:`examples`
* **Statistical Model**: :ref:`statistical_model`
* **Installation**: :ref:`installation`

Support Channels
~~~~~~~~~~~~~~~~

* **GitHub Issues**: https://github.com/nicenzhou/edge-gwas/issues
* **GitHub Discussions**: https://github.com/nicenzhou/edge-gwas/discussions
* **Code Questions**: jyzhou@stanford.edu
* **Research Questions**: molly.hall@pennmedicine.upenn.edu
