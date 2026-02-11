.. _troubleshooting:

Troubleshooting Guide
===================

This guide covers potential troubleshooting walkthroughs in edge-gwas.

Common Issues and Solutions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Issue 1: High proportion of SNPs failing to converge (>10%)**

*Symptoms:* Many SNPs in skipped_snps.log

*Solutions:*

.. code-block:: python

   # Try more robust optimization
   edge = EDGEAnalysis(
       outcome_type='continuous',
       ols_method='nm',  # Nelder-Mead is more robust
       max_iter=5000
   )
   
   # Or increase MAF threshold
   genotype_filtered = genotype.loc[:, maf > 0.05]

**Issue 2: Genomic inflation (λ > 1.05)**

*Symptoms:* Systematic p-value inflation, Q-Q plot deviation

*Solutions:*

.. code-block:: python

   # Add more PCs
   covariates_extended = ['age', 'sex'] + [f'pc{i}' for i in range(1, 21)]
   
   # Check for population stratification
   # Verify sample homogeneity

**Issue 3: Unstable alpha estimates**

*Symptoms:* α values vary wildly between nearby SNPs

*Solutions:*

.. code-block:: python

   # Increase training set size
   train_ids, test_ids = train_test_split(
       genotype.index, test_size=0.4, random_state=42
   )
   
   # Filter low-MAF variants
   genotype_filtered = genotype.loc[:, maf > 0.05]

**Issue 4: Transformation produces NaN/Inf values**

*Symptoms:* Warning messages about invalid values

*Solutions:*

.. code-block:: python

   # For log transformation with zeros
   phenotype['outcome_shifted'] = phenotype['outcome'] + 1
   
   # Use rank-based transformation (more robust)
   edge = EDGEAnalysis(
       outcome_type='continuous',
       outcome_transform='rank_inverse_normal'
   )
   
   # Remove extreme outliers before transformation
   from scipy import stats
   z_scores = np.abs(stats.zscore(phenotype['outcome']))
   phenotype_clean = phenotype[z_scores < 5]

**Issue 5: Memory errors with large datasets**

*Symptoms:* Out-of-memory errors

*Solutions:*

.. code-block:: python

   # Use L-BFGS (lower memory)
   edge = EDGEAnalysis(
       outcome_type='continuous',
       ols_method='lbfgs'
   )
   
   # Process chromosomes separately
   for chrom in range(1, 23):
       geno_chr = genotype.filter(regex=f'^chr{chrom}:')
       alpha_chr, gwas_chr = edge.run_full_analysis(
           train_geno_chr, train_p, test_geno_chr, test_p,
           outcome, covariates,
           output_prefix=f'results/chr{chrom}'
       )

**Issue 6: Slow computation**

*Symptoms:* Analysis taking too long

*Solutions:*

.. code-block:: python

   # Use faster optimization method
   edge = EDGEAnalysis(
       outcome_type='continuous',
       ols_method='newton',  # Faster when converging
       max_iter=500
   )
   
   # Pre-filter variants
   genotype_filtered = genotype.loc[:, (maf > 0.01) & (hwe_p > 1e-6)]
   
   # Use parallel processing (set n_jobs in future version)
   # Process chromosomes in parallel

Common Errors and Solutions
----------------------------

Error 1: "ValueError: No common samples found between analysis data and GRM"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Cause:** Sample IDs don't match between genotype data and GRM.

**Solution:**

.. code-block:: python

   # Check sample ID formats
   print("Genotype sample IDs (first 5):")
   print(genotype.index[:5])
   
   print("\nGRM sample IDs (first 5):")
   print(grm_sample_ids['IID'][:5])
   
   # Ensure consistent format (string vs integer)
   genotype.index = genotype.index.astype(str)
   grm_sample_ids['IID'] = grm_sample_ids['IID'].astype(str)
   
   # Or rename to match
   # If GRM uses 'FID_IID' format:
   grm_sample_ids['sample_id'] = grm_sample_ids['FID'] + '_' + grm_sample_ids['IID']

Error 2: "LinAlgError: Matrix is singular"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Cause:** Perfect multicollinearity in covariates or genotypes.

**Solution:**

.. code-block:: python

   # Check covariate correlations
   import seaborn as sns
   
   corr_matrix = phenotype[covariates].corr()
   sns.heatmap(corr_matrix, annot=True)
   plt.savefig('covariate_correlation.png')
   
   # Remove highly correlated covariates
   def remove_collinear(df, threshold=0.95):
       """Remove covariates with correlation > threshold."""
       corr = df.corr().abs()
       upper = corr.where(np.triu(np.ones(corr.shape), k=1).astype(bool))
       to_drop = [col for col in upper.columns if any(upper[col] > threshold)]
       return df.drop(columns=to_drop)
   
   covariates_clean = remove_collinear(phenotype[covariates])

Error 3: "ConvergenceWarning: Maximum iterations reached"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Cause:** Model failed to converge within max_iter iterations.

**Solution:**

.. code-block:: python

   # Increase max iterations
   edge = EDGEAnalysis(
       outcome_type='continuous',
       max_iter=5000,  # Increased from 1000
       ols_method='bfgs'
   )
   
   # Or try more robust method
   edge = EDGEAnalysis(
       outcome_type='continuous',
       max_iter=2000,
       ols_method='nm'  # Nelder-Mead is more robust
   )
   
   # Check which SNPs failed
   skipped = edge.get_skipped_snps()
   
   # Examine skipped SNPs for common characteristics
   skipped_info = variant_info[variant_info['variant_id'].isin(skipped)]
   print("MAF distribution of skipped SNPs:")
   print(skipped_info['maf'].describe())

Error 4: "ValueError: Log transformation requires all positive values"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Cause:** Attempting log transformation on data with zeros or negative values.

**Solution:**

.. code-block:: python

   # Check for zeros/negatives
   print(f"Minimum value: {phenotype['outcome'].min()}")
   print(f"Zeros: {(phenotype['outcome'] == 0).sum()}")
   print(f"Negatives: {(phenotype['outcome'] < 0).sum()}")
   
   # Solution 1: Add constant before log
   phenotype['outcome_shifted'] = phenotype['outcome'] + 1
   
   # Solution 2: Use log1p (log(1+x))
   phenotype['outcome_log1p'] = np.log1p(phenotype['outcome'])
   
   # Solution 3: Use rank-based transformation instead
   edge = EDGEAnalysis(
       outcome_type='continuous',
       outcome_transform='rank_inverse_normal'  # Handles any distribution
   )

Error 5: "KeyError: 'variant_id' not found"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Cause:** Missing required columns in input data.

**Solution:**

.. code-block:: python

   # Check required columns
   print("Genotype columns:", genotype.columns[:5])
   print("Phenotype columns:", phenotype.columns)
   
   # Ensure variant_info has correct structure
   if variant_info is not None:
       required_cols = ['variant_id', 'ref_allele', 'alt_allele']
       missing = [col for col in required_cols if col not in variant_info.columns]
       if missing:
           print(f"Missing columns in variant_info: {missing}")
           
           # Add missing columns with defaults
           if 'ref_allele' not in variant_info.columns:
               variant_info['ref_allele'] = 'REF'
           if 'alt_allele' not in variant_info.columns:
               variant_info['alt_allele'] = 'ALT'

Error 6: "MemoryError: Unable to allocate array"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Cause:** Insufficient memory for large dataset.

**Solution:**

.. code-block:: python

   # Solution 1: Process by chromosome
   for chrom in range(1, 23):
       geno_chr = genotype.filter(regex=f'^chr{chrom}:')
       # Run EDGE on chromosome
   
   # Solution 2: Use chunking (see Q20 above)
   
   # Solution 3: Use L-BFGS (lower memory)
   edge = EDGEAnalysis(
       outcome_type='continuous',
       ols_method='lbfgs'  # More memory-efficient
   )
   
   # Solution 4: Reduce number of covariates
   # Use fewer PCs (e.g., 10 instead of 20)

Error 7: "RuntimeWarning: invalid value encountered in divide"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Cause:** Division by zero or very small numbers in alpha calculation.

**Solution:**

.. code-block:: python

   # This is automatically handled by EDGE
   # SNPs with coef_hom ≈ 0 are skipped
   
   # But you can check which SNPs caused this
   skipped = edge.get_skipped_snps()
   
   # Examine these SNPs
   for snp in skipped[:10]:  # First 10
       snp_data = genotype[snp]
       print(f"\n{snp}:")
       print(f"  Genotype counts: {snp_data.value_counts()}")
       print(f"  MAF: {calculate_maf(snp_data):.4f}")
   
   # Usually indicates very rare variants or monomorphic SNPs

Performance Optimization Tips
------------------------------

Speed Optimization
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # 1. Pre-filter variants before analysis
   # Remove rare variants, HWE failures upfront
   qc_pass = (maf > 0.01) & (hwe_p > 1e-6) & (call_rate > 0.95)
   genotype_filtered = genotype.loc[:, qc_pass]
   
   # 2. Use faster optimization method
   edge = EDGEAnalysis(
       outcome_type='continuous',
       ols_method='newton',  # Faster when it converges
       max_iter=500  # Reduce if most SNPs converge quickly
   )
   
   # 3. Reduce covariate set if possible
   # Start with essential covariates
   covariates_minimal = ['age', 'sex', 'pc1', 'pc2', 'pc3']
   
   # 4. Process chromosomes in parallel
   # See Q16 for parallel processing example
   
   # 5. Use binary file formats
   # Parquet is faster than CSV for large datasets
   genotype.to_parquet('genotypes.parquet')
   genotype = pd.read_parquet('genotypes.parquet')

Memory Optimization
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # 1. Use appropriate data types
   # Convert genotypes to int8 (saves 7/8 memory vs float64)
   genotype = genotype.astype('int8')
   
   # 2. Process in chunks (see Q20)
   
   # 3. Delete large objects when done
   import gc
   
   del train_genotype
   gc.collect()
   
   # 4. Use generators for very large datasets
   def genotype_generator(file_path, chunk_size=1000):
       """Yield genotype chunks from file."""
       for chunk in pd.read_csv(file_path, chunksize=chunk_size):
           yield chunk
   
   # 5. Monitor memory usage
   import psutil
   import os
   
   process = psutil.Process(os.getpid())
   print(f"Memory usage: {process.memory_info().rss / 1e9:.2f} GB")

Accuracy Optimization
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # 1. Use larger training set for stable alpha
   train_ids, test_ids = train_test_split(
       genotype.index, 
       test_size=0.4,  # 60% training
       random_state=42
   )
   
   # 2. Filter low-quality variants
   # Stricter MAF threshold
   genotype_hq = genotype.loc[:, maf > 0.05]
   
   # 3. Use appropriate transformation
   # Check normality before/after
   from scipy import stats
   
   _, p_before = stats.shapiro(phenotype['outcome'].sample(5000))
   
   # Try transformation
   outcome_transformed = transform_outcome(phenotype['outcome'])
   _, p_after = stats.shapiro(outcome_transformed.sample(5000))
   
   print(f"Normality p-value: {p_before:.4f} → {p_after:.4f}")
   
   # 4. Verify convergence
   # Check proportion of successful fits
   success_rate = 1 - len(edge.get_skipped_snps()) / n_variants
   print(f"Convergence rate: {success_rate:.1%}")
   
   if success_rate < 0.9:
       print("Warning: Low convergence rate. Consider:")
       print("- Increasing max_iter")
       print("- Trying different optimization method")
       print("- Filtering rare variants")

Reproducibility Best Practices
-------------------------------

Setting Random Seeds
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import numpy as np
   import random
   
   # Set all random seeds
   RANDOM_SEED = 42
   
   random.seed(RANDOM_SEED)
   np.random.seed(RANDOM_SEED)
   
   # For train/test split
   train_ids, test_ids = train_test_split(
       genotype.index,
       test_size=0.5,
       random_state=RANDOM_SEED
   )

Version Control
~~~~~~~~~~~~~~~

.. code-block:: python

   # Record versions in analysis log
   import edge_gwas
   import pandas as pd
   import numpy as np
   import statsmodels
   
   print("Software versions:")
   print(f"edge-gwas: {edge_gwas.__version__}")
   print(f"pandas: {pd.__version__}")
   print(f"numpy: {np.__version__}")
   print(f"statsmodels: {statsmodels.__version__}")
   
   # Save to file
   with open('software_versions.txt', 'w') as f:
       f.write(f"edge-gwas: {edge_gwas.__version__}\n")
       f.write(f"pandas: {pd.__version__}\n")
       f.write(f"numpy: {np.__version__}\n")
       f.write(f"statsmodels: {statsmodels.__version__}\n")

Parameter Documentation
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Document all analysis parameters
   analysis_params = {
       'outcome_type': 'continuous',
       'outcome_transform': 'rank_inverse_normal',
       'ols_method': 'bfgs',
       'max_iter': 1000,
       'train_test_split': 0.5,
       'random_seed': 42,
       'maf_threshold': 0.01,
       'hwe_threshold': 1e-6,
       'covariates': ['age', 'sex', 'bmi'] + [f'pc{i}' for i in range(1, 11)],
       'n_samples_train': len(train_ids),
       'n_samples_test': len(test_ids),
       'n_variants': len(genotype.columns)
   }
   
   # Save parameters
   import json
   with open('analysis_parameters.json', 'w') as f:
       json.dump(analysis_params, f, indent=2)


Appendix A: Troubleshooting Decision Tree
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Flowchart for Common Issues:**

.. code-block:: text

   START: EDGE Analysis Problem
        |
        ├─> High λ (>1.05)?
        |     YES → Add more PCs → Still high?
        |              YES → Check relatedness → Remove related samples
        |              NO  → RESOLVED
        |     NO ↓
        |
        ├─> Low convergence (<90%)?
        |     YES → Check MAF distribution
        |           ├─> Many rare variants?
        |           |     YES → Filter MAF>0.05 → RESOLVED
        |           |     NO ↓
        |           ├─> Try different method
        |           |     ols_method='nm' → Converges?
        |           |           YES → RESOLVED
        |           |           NO → Increase max_iter → RESOLVED
        |     NO ↓
        |
        ├─> Extreme alpha values?
        |     YES → Check specific SNPs
        |           ├─> Monomorphic?
        |           |     YES → Filter → RESOLVED
        |           ├─> Low MAF?
        |           |     YES → Increase MAF threshold → RESOLVED
        |           ├─> Covariate issue?
        |                 YES → Check multicollinearity → RESOLVED
        |     NO ↓
        |
        ├─> Memory errors?
        |     YES → Process by chromosome → RESOLVED
        |           OR use ols_method='lbfgs' → RESOLVED
        |     NO ↓
        |
        ├─> No significant hits?
        |     YES → Expected for some traits
        |           ├─> Check power
        |           ├─> Compare with additive GWAS
        |           └─> Consider larger sample size
        |
        └─> Other issues → Consult documentation
                         → Post GitHub issue

See Also
--------

**Documentation:**

* `Documentation Home <index.html>`_ - Home
* :ref:`installation` - Installation instructions and requirements
* :ref:`quickstart` - Getting started guide with simple examples
* :ref:`statistical_model` - Statistical methods and mathematical background
* :ref:`examples` - Example analyses and case studies
* :ref:`visualization` - Plotting and visualization guide
* :ref:`api_reference` - Complete API documentation
* :ref:`troubleshooting` - Troubleshooting guide and common issues
* :ref:`faq` - Frequently asked questions
* :ref:`citation` - How to cite EDGE in publications
* :ref:`changelog` - Version history and release notes
* :ref:`futureupdates` - Planned features and roadmap

---

*Last updated: 2026-02-10 for edge-gwas v0.1.2*

*For questions or issues, visit:* https://github.com/nicenzhou/edge-gwas/issues
