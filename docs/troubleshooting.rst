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

Frequently Asked Questions (FAQ)
---------------------------------

General Questions
~~~~~~~~~~~~~~~~~

**Q1: When should I use EDGE instead of standard additive GWAS?**

A: Use EDGE when:

- You suspect nonadditive inheritance patterns
- Previous additive GWAS found few associations
- You have large sample size (N > 2,000 under the case to control ratio at 70/30)
- You're following up candidate genes with known non-additive effects
- You want to explore genetic architecture

Standard additive GWAS is sufficient for initial exploratory analyses or when sample size is limited.

**Q2: What sample size do I need for EDGE analysis?**

A: Minimum recommendations:

- Binary outcomes: 1,400 cases + 600 controls
- Continuous outcomes: 2,000 total samples

Optimal recommendations:

- Binary outcomes: 5,000+ cases and controls
- Continuous outcomes: 10,000+ samples

Larger samples provide more stable α estimates and better power.

**Q3: How do I choose between optimization methods?**

A: Decision guide:

- Start with default (BFGS) - works for >95% of cases
- Use L-BFGS for very large datasets (N > 100,000)
- Try Newton if BFGS has convergence issues
- Use Nelder-Mead for maximum robustness (slower)

**Q4: Should I always use outcome transformations?**

A: No, only when needed:

- Check outcome distribution first (histogram, Q-Q plot)
- Use transformations for non-normal distributions
- RINT is safest choice for unknown distributions
- No transformation needed if outcome is already normal

Technical Questions
~~~~~~~~~~~~~~~~~~~
**Q5: Why do some SNPs fail to converge?**

A: Common reasons:

- **Low MAF**: Rare variants have insufficient heterozygotes
- **Perfect separation**: In logistic regression, genotype perfectly predicts outcome
- **Multicollinearity**: Genotype correlated with covariates
- **Numerical instability**: Extreme coefficient values

Solutions:

.. code-block:: python

   # 1. Filter rare variants
   genotype_filtered = genotype.loc[:, maf > 0.05]
   
   # 2. Try different optimization method
   edge = EDGEAnalysis(
       outcome_type='continuous',
       ols_method='nm',  # More robust
       max_iter=5000
   )
   
   # 3. Check skipped SNPs
   skipped = edge.get_skipped_snps()
   # Review for patterns (e.g., all in same region, all rare)

**Q6: How do I interpret alpha values?**

A: Alpha interpretation guide:

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Alpha
     - Pattern
     - Biological Example
   * - α ≈ 0
     - Recessive
     - Cystic fibrosis (CFTR mutations)
   * - α ≈ 0.5
     - Additive
     - Height, most complex traits
   * - α ≈ 1
     - Dominant
     - Huntington's disease (HTT expansion)
   * - α > 1
     - Over-dominant
     - Sickle cell trait (heterozygote advantage)
   * - α < 0
     - Under-recessive
     - Heterozygote disadvantage (rare)

**Q7: What if lambda (genomic inflation) is high?**

A: Troubleshooting high λ (>1.05):

.. code-block:: python

   # 1. Add more principal components
   covariates_extended = covariates + [f'pc{i}' for i in range(11, 21)]
   
   # 2. Check for stratification
   # Run PCA and examine PC plots by case/control status
   
   # 3. Verify QC filters
   # - Check for related individuals
   # - Verify genotype quality
   # - Check for batch effects
   
   # 4. Consider genomic control correction (post-hoc)
   gwas_results['pval_gc'] = adjust_genomic_control(
       gwas_results['pval'], 
       lambda_gc
   )

**Q8: How do I handle related individuals?**

A: Current options:

1. **Remove related individuals** (recommended for v0.1.1):

   .. code-block:: python
   
      # Identify unrelated set using KING or PLINK
      # kinship threshold: 0.0884 (3rd degree relatives)
      unrelated_ids = identify_unrelated(kinship_matrix, threshold=0.0884)
      genotype_unrelated = genotype.loc[unrelated_ids]

2. **Wait for mixed model support** (planned for v0.2.0):

   - Will account for relatedness using GRM
   - Allows use of full dataset including families
   - More powerful for related samples

**Q9: Can I use EDGE with imputed genotypes?**

A: Yes, with considerations:

.. code-block:: python

   # 1. Filter by imputation quality
   high_quality = imputation_info > 0.8  # INFO score
   genotype_imputed = genotype.loc[:, high_quality]
   
   # 2. Use dosages (0-2 continuous) instead of hard calls
   # EDGE works with dosages
   
   # 3. Be aware of uncertainty
   # Imputed variants may have less stable alpha estimates
   # Prefer directly genotyped SNPs for alpha calculation

**Q10: How do I replicate EDGE findings in an independent cohort?**

A: Replication workflow:

.. code-block:: python

   # Step 1: Extract alpha values from discovery cohort
   discovery_alpha = pd.read_csv('discovery_alpha_values.csv')
   alpha_dict = dict(zip(
       discovery_alpha['variant_id'],
       discovery_alpha['alpha_value']
   ))
   
   # Step 2: Apply same alpha values to replication cohort
   edge_replication = EDGEAnalysis(outcome_type='continuous')
   
   replication_results = edge_replication.apply_alpha(
       genotype_data=replication_geno,
       phenotype_df=replication_pheno,
       outcome='trait',
       covariates=['age', 'sex', 'pc1', 'pc2'],
       alpha_values=discovery_alpha  # Use discovery alphas
   )
   
   # Step 3: Check for consistency
   sig_variants = discovery_gwas[discovery_gwas['pval'] < 5e-8]['variant_id']
   
   replication_sig = replication_results[
       replication_results['variant_id'].isin(sig_variants)
   ]
   
   # Count replicated loci (p < 0.05, same direction)
   replicated = (
       (replication_sig['pval'] < 0.05) &
       (np.sign(replication_sig['coef']) == np.sign(discovery_gwas['coef']))
   ).sum()
   
   print(f"Replication rate: {replicated / len(sig_variants):.1%}")

Interpretation Questions
~~~~~~~~~~~~~~~~~~~~~~~~

**Q11: Can EDGE detect epistasis (SNP-SNP interactions)?**

A: No, EDGE models single-SNP effects:

- EDGE captures nonadditive effects at **single loci** (dominance)
- Does not model **interactions between loci** (epistasis)
- For epistasis, use dedicated interaction testing methods
- EDGE and epistasis testing are complementary approaches

**Q12: How do I know if EDGE found a truly novel locus?**

A: Validation checklist:

.. code-block:: python

   # 1. Check it's not in LD with known loci
   def check_ld_independence(novel_snp, known_snps, ld_matrix):
       """Check if novel SNP is independent of known loci."""
       max_ld_r2 = 0
       for known_snp in known_snps:
           r2 = ld_matrix[novel_snp, known_snp]
           if r2 > max_ld_r2:
               max_ld_r2 = r2
       return max_ld_r2 < 0.1  # Independent if r² < 0.1
   
   # 2. Verify it's not a QC artifact
   # - Check MAF, HWE, call rate
   # - Look for clustered missingness
   # - Verify in multiple genotyping platforms
   
   # 3. Examine alpha biological plausibility
   # - Does inheritance pattern make sense?
   # - Are nearby genes relevant to phenotype?
   
   # 4. Replicate in independent cohort
   # - Must replicate with similar alpha value
   # - Check effect direction consistency

**Q13: What does it mean if alpha is negative?**

A: Negative alpha indicates **under-recessive** pattern:

- Heterozygotes have **opposite effect direction** from homozygotes
- Example: If homozygotes increase trait, heterozygotes decrease it
- Rare but possible biological scenarios:
  
  - Antagonistic pleiotropy
  - Complex regulatory mechanisms
  - Possible false positive (check carefully!)

Recommendation: Scrutinize negative alpha SNPs carefully for QC issues.

**Q14: Should I filter SNPs based on alpha values?**

A: Generally no, but consider:

**Don't filter** alpha values for:

- Primary GWAS analysis (test all variants)
- Exploratory analysis
- Novel discovery

**Consider filtering** for:

- Extremely unstable alphas (|α| > 5): likely convergence issues
- Follow-up functional studies: focus on biologically plausible patterns
- Replication: use only well-estimated alphas (SE < 0.5)

.. code-block:: python

   # Flag potentially unstable alphas
   alpha_df['stable'] = (
       (alpha_df['alpha_value'].abs() < 5) &
       (alpha_df['std_err_het'] < 0.5) &
       (alpha_df['std_err_hom'] < 0.5)
   )
   
   print(f"Stable alpha estimates: {alpha_df['stable'].sum()} / {len(alpha_df)}")

**Q15: How do I compare EDGE results across different traits?**

A: Cross-trait comparison:

.. code-block:: python

   # Load results for multiple traits
   traits = ['bmi', 'height', 'triglycerides']
   alpha_by_trait = {}
   
   for trait in traits:
       alpha_df = pd.read_csv(f'results/{trait}_alpha_values.csv')
       alpha_by_trait[trait] = alpha_df
   
   # Find variants tested in all traits
   common_variants = set(alpha_by_trait[traits[0]]['variant_id'])
   for trait in traits[1:]:
       common_variants &= set(alpha_by_trait[trait]['variant_id'])
   
   # Compare alpha values
   comparison = pd.DataFrame({
       trait: alpha_by_trait[trait].set_index('variant_id')['alpha_value']
       for trait in traits
   }).loc[list(common_variants)]
   
   # Calculate correlation
   print("Alpha correlation across traits:")
   print(comparison.corr())
   
   # Identify trait-specific non-additive effects
   # (different alpha values across traits)
   comparison['alpha_range'] = comparison.max(axis=1) - comparison.min(axis=1)
   trait_specific = comparison.nlargest(20, 'alpha_range')
   print("\nVariants with trait-specific inheritance patterns:")
   print(trait_specific)

Practical Implementation Questions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Q16: Can I run EDGE on multiple chromosomes in parallel?**

A: Yes, recommended for large datasets:

.. code-block:: python

   from joblib import Parallel, delayed
   
   def analyze_chromosome(chrom, train_g, train_p, test_g, test_p):
       """Run EDGE on a single chromosome."""
       # Filter to chromosome
       chrom_snps = [col for col in train_g.columns if col.startswith(f'chr{chrom}:')]
       
       edge = EDGEAnalysis(outcome_type='continuous')
       alpha_df, gwas_df = edge.run_full_analysis(
           train_genotype=train_g[chrom_snps],
           train_phenotype=train_p,
           test_genotype=test_g[chrom_snps],
           test_phenotype=test_p,
           outcome='trait',
           covariates=['age', 'sex'],
           output_prefix=f'results/chr{chrom}'
       )
       
       return alpha_df, gwas_df
   
   # Run in parallel
   results = Parallel(n_jobs=22)(
       delayed(analyze_chromosome)(chrom, train_g, train_p, test_g, test_p)
       for chrom in range(1, 23)
   )
   
   # Combine results
   all_alpha = pd.concat([r[0] for r in results], ignore_index=True)
   all_gwas = pd.concat([r[1] for r in results], ignore_index=True)

**Q17: How much disk space do I need?**

A: Storage requirements:

.. code-block:: python

   # Estimates for typical GWAS
   
   # Input data
   # - Genotypes: ~4 bytes per genotype
   #   Example: 500K SNPs × 50K samples = 100 GB
   
   # - Phenotypes: negligible (~1 MB)
   
   # Output data
   # - Alpha values: ~500 bytes per SNP
   #   Example: 500K SNPs = 250 MB
   
   # - GWAS results: ~1 KB per SNP
   #   Example: 500K SNPs = 500 MB
   
   # Total: ~100 GB input + 1 GB output
   
   # Recommendations:
   # - 200 GB free space for safety
   # - SSD for faster I/O
   # - Compress intermediate files

**Q18: Can I use EDGE with X chromosome variants?**

A: Limited support in v0.1.1:

.. code-block:: python

   # Current workaround for X chromosome
   
   # 1. Analyze males and females separately
   males = phenotype[phenotype['sex'] == 'male'].index
   females = phenotype[phenotype['sex'] == 'female'].index
   
   # 2. For males: recode genotypes (0/1 → 0/2)
   #    Males are hemizygous, so 0=ref, 1=alt (recode as 0 and 2)
   genotype_X_males = genotype_X.loc[males].replace({1: 2})
   
   # 3. Run EDGE separately
   edge_females = EDGEAnalysis(outcome_type='continuous')
   alpha_f, gwas_f = edge_females.run_full_analysis(
       train_genotype_X.loc[females_train],
       phenotype.loc[females_train],
       test_genotype_X.loc[females_test],
       phenotype.loc[females_test],
       outcome, covariates
   )
   
   # Note: Proper X chromosome support coming in v0.2.0

**Q19: How do I save and load EDGE results for later use?**

A: Best practices for saving results:

.. code-block:: python

   # Save results
   alpha_df.to_csv('alpha_values.csv', index=False)
   gwas_results.to_csv('gwas_results.csv', index=False)
   
   # For large datasets, use compressed format
   alpha_df.to_csv('alpha_values.csv.gz', index=False, compression='gzip')
   
   # Or use Parquet for faster loading
   alpha_df.to_parquet('alpha_values.parquet')
   gwas_results.to_parquet('gwas_results.parquet')
   
   # Load results later
   alpha_df = pd.read_csv('alpha_values.csv')
   # or
   alpha_df = pd.read_parquet('alpha_values.parquet')
   
   # Apply saved alphas to new data
   edge = EDGEAnalysis(outcome_type='continuous')
   new_gwas = edge.apply_alpha(
       new_genotype,
       new_phenotype,
       outcome,
       covariates,
       alpha_values=alpha_df  # Use saved alphas
   )

**Q20: What are the memory requirements?**

A: Memory usage guide:

.. code-block:: python

   # Rough estimates
   
   # Genotype data: 8 bytes per value (pandas float64)
   # Example: 100K SNPs × 10K samples = 8 GB
   
   # For large datasets:
   # 1. Process by chromosome
   # 2. Use chunking
   
   def analyze_in_chunks(genotype, phenotype, chunk_size=10000):
       """Process SNPs in chunks to reduce memory."""
       all_results = []
       
       n_snps = len(genotype.columns)
       n_chunks = (n_snps + chunk_size - 1) // chunk_size
       
       for i in range(n_chunks):
           start_idx = i * chunk_size
           end_idx = min((i + 1) * chunk_size, n_snps)
           
           chunk_snps = genotype.columns[start_idx:end_idx]
           geno_chunk = genotype[chunk_snps]
           
           edge = EDGEAnalysis(outcome_type='continuous')
           alpha_chunk, gwas_chunk = edge.run_full_analysis(
               geno_chunk.loc[train_ids],
               phenotype.loc[train_ids],
               geno_chunk.loc[test_ids],
               phenotype.loc[test_ids],
               outcome='trait',
               covariates=['age', 'sex']
           )
           
           all_results.append((alpha_chunk, gwas_chunk))
           
           # Clear memory
           del geno_chunk
           import gc
           gc.collect()
       
       # Combine results
       alpha_combined = pd.concat([r[0] for r in all_results], ignore_index=True)
       gwas_combined = pd.concat([r[1] for r in all_results], ignore_index=True)
       
       return alpha_combined, gwas_combined

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







Appendix E: Optimization Method Details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**BFGS Algorithm:**

The Broyden-Fletcher-Goldfarb-Shannon algorithm approximates the Hessian matrix using gradient information:

.. math::

   H_{k+1} = H_k + \frac{y_k y_k^T}{y_k^T s_k} - \frac{H_k s_k s_k^T H_k}{s_k^T H_k s_k}

where:

- :math:`s_k = x_{k+1} - x_k` (step)
- :math:`y_k = \nabla f(x_{k+1}) - \nabla f(x_k)` (gradient change)
- :math:`H_k` is the approximate Hessian at iteration k

**Convergence criteria:**

.. math::

   \|\nabla f(x_k)\| < \epsilon \quad \text{or} \quad |f(x_k) - f(x_{k-1})| < \epsilon

Default tolerance: :math:`\epsilon = 10^{-6}`

**L-BFGS Algorithm:**

Limited-memory BFGS stores only m recent vector pairs (typically m=10):

Memory requirement: :math:`O(m \cdot n)` instead of :math:`O(n^2)` for full BFGS

Trade-off: Slightly slower convergence but much lower memory usage

**Newton-Raphson:**

Uses exact Hessian matrix:

.. math::

   x_{k+1} = x_k - H^{-1}(x_k) \nabla f(x_k)

where :math:`H(x_k)` is the true Hessian matrix.

**Convergence rate:** Quadratic (fastest when close to optimum)

**Cost:** Expensive Hessian computation at each iteration

**Nelder-Mead Simplex:**

Derivative-free method using simplex of n+1 points:

Operations: reflection, expansion, contraction, shrinkage

**Pros:** Robust, no gradient needed
**Cons:** Slow, may not converge for high dimensions

Appendix F: Sample Size Calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Power calculation for EDGE:**

.. code-block:: python

   from scipy.stats import norm, nct
   import numpy as np
   
   def edge_power_calculation(n, maf, alpha_true, beta, sig_level=5e-8):
       """
       Calculate statistical power for EDGE analysis.
       
       Args:
           n: Sample size
           maf: Minor allele frequency
           alpha_true: True alpha value
           beta: Effect size
           sig_level: Significance threshold
       
       Returns:
           Statistical power (0-1)
       """
       # Genotype frequencies under HWE
       p_00 = (1 - maf) ** 2
       p_01 = 2 * maf * (1 - maf)
       p_11 = maf ** 2
       
       # EDGE encoding variance
       edge_00 = 1.0
       edge_01 = alpha_true
       edge_11 = 0.0
       
       mean_edge = p_00 * edge_00 + p_01 * edge_01 + p_11 * edge_11
       var_edge = (p_00 * edge_00**2 + p_01 * edge_01**2 + 
                   p_11 * edge_11**2 - mean_edge**2)
       
       # Non-centrality parameter
       ncp = beta * np.sqrt(n * var_edge)
       
       # Critical value
       z_crit = norm.ppf(1 - sig_level / 2)
       
       # Power
       power = 1 - norm.cdf(z_crit - ncp) + norm.cdf(-z_crit - ncp)
       
       return power
   
   # Example: Power curves
   import matplotlib.pyplot as plt
   
   sample_sizes = np.arange(1000, 50000, 1000)
   alpha_values = [0, 0.25, 0.5, 0.75, 1.0]
   
   fig, ax = plt.subplots(figsize=(10, 6))
   
   for alpha in alpha_values:
       powers = [edge_power_calculation(n, maf=0.3, alpha_true=alpha, beta=0.1)
                 for n in sample_sizes]
       ax.plot(sample_sizes, powers, label=f'α={alpha}')
   
   ax.axhline(0.8, color='gray', linestyle='--', label='80% power')
   ax.set_xlabel('Sample Size')
   ax.set_ylabel('Statistical Power')
   ax.set_title('EDGE Power by Inheritance Pattern\n(MAF=0.3, β=0.1)')
   ax.legend()
   ax.grid(True, alpha=0.3)
   plt.savefig('power_analysis.png', dpi=300)

**Sample size requirement:**

For 80% power at genome-wide significance (p < 5×10⁻⁸):

.. code-block:: python

   def required_sample_size(maf, alpha_true, beta, power=0.8, sig_level=5e-8):
       """
       Calculate required sample size for desired power.
       
       Returns:
           Required sample size
       """
       from scipy.optimize import fsolve
       
       def power_equation(n):
           return edge_power_calculation(n, maf, alpha_true, beta, sig_level) - power
       
       n_required = fsolve(power_equation, x0=10000)[0]
       return int(np.ceil(n_required))
   
   # Examples
   scenarios = [
       {'maf': 0.3, 'alpha': 0.2, 'beta': 0.1, 'label': 'Recessive'},
       {'maf': 0.3, 'alpha': 0.5, 'beta': 0.1, 'label': 'Additive'},
       {'maf': 0.3, 'alpha': 0.8, 'beta': 0.1, 'label': 'Dominant'},
   ]
   
   print("Sample size requirements for 80% power:")
   print("=" * 60)
   for scenario in scenarios:
       n_req = required_sample_size(
           scenario['maf'], 
           scenario['alpha'], 
           scenario['beta']
       )
       print(f"{scenario['label']:12s}: N = {n_req:,}")

Appendix G: Useful Helper Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Calculate Minor Allele Frequency:**

.. code-block:: python

   def calculate_maf(genotype_series):
       """
       Calculate minor allele frequency from genotype data.
       
       Args:
           genotype_series: Pandas Series with genotypes (0, 1, 2)
       
       Returns:
           MAF value
       """
       valid_geno = genotype_series.dropna()
       
       # Count alleles
       n_alt = (2 * (valid_geno == 2).sum() + 
                (valid_geno == 1).sum())
       n_total = 2 * len(valid_geno)
       
       # Frequency
       freq_alt = n_alt / n_total
       
       # MAF is minimum of alt and ref frequency
       maf = min(freq_alt, 1 - freq_alt)
       
       return maf

**Calculate Hardy-Weinberg Equilibrium p-value:**

.. code-block:: python

   def calculate_hwe(genotype_series, case_control=None):
       """
       Test for Hardy-Weinberg Equilibrium.
       
       Args:
           genotype_series: Pandas Series with genotypes
           case_control: Optional Series indicating cases (1) and controls (0)
                        If provided, tests HWE in controls only
       
       Returns:
           HWE p-value
       """
       from scipy.stats import chi2
       
       # Use controls only if case/control provided
       if case_control is not None:
           geno = genotype_series[case_control == 0]
       else:
           geno = genotype_series
       
       # Count genotypes
       n_00 = (geno == 0).sum()
       n_01 = (geno == 1).sum()
       n_11 = (geno == 2).sum()
       n_total = n_00 + n_01 + n_11
       
       # Allele frequency
       p_alt = (2 * n_11 + n_01) / (2 * n_total)
       p_ref = 1 - p_alt
       
       # Expected counts under HWE
       exp_00 = n_total * p_ref ** 2
       exp_01 = n_total * 2 * p_ref * p_alt
       exp_11 = n_total * p_alt ** 2
       
       # Chi-square test (1 df)
       chi2_stat = ((n_00 - exp_00) ** 2 / exp_00 +
                    (n_01 - exp_01) ** 2 / exp_01 +
                    (n_11 - exp_11) ** 2 / exp_11)
       
       p_value = 1 - chi2.cdf(chi2_stat, df=1)
       
       return p_value

**Calculate Genomic Inflation Factor:**

.. code-block:: python

   def calculate_lambda(pvalues):
       """
       Calculate genomic inflation factor (lambda).
       
       Args:
           pvalues: Array of p-values
       
       Returns:
           Lambda value
       """
       from scipy.stats import chi2
       
       # Remove NaN
       pvals_clean = pvalues[~np.isnan(pvalues)]
       
       # Convert to chi-square statistics
       chi2_stats = chi2.ppf(1 - pvals_clean, df=1)
       
       # Lambda is ratio of observed to expected median
       lambda_gc = np.median(chi2_stats) / chi2.ppf(0.5, df=1)
       
       return lambda_gc

**Apply Genomic Control:**

.. code-block:: python

   def apply_genomic_control(pvalues, lambda_gc):
       """
       Adjust p-values using genomic control.
       
       Args:
           pvalues: Original p-values
           lambda_gc: Genomic inflation factor
       
       Returns:
           Adjusted p-values
       """
       from scipy.stats import chi2
       
       # Convert to chi-square
       chi2_stats = chi2.ppf(1 - pvalues, df=1)
       
       # Adjust by lambda
       chi2_adjusted = chi2_stats / lambda_gc
       
       # Convert back to p-values
       pvalues_adjusted = 1 - chi2.cdf(chi2_adjusted, df=1)
       
       return pvalues_adjusted

**Identify Unrelated Samples:**

.. code-block:: python

   def identify_unrelated_samples(kinship_matrix, threshold=0.0884):
       """
       Identify maximal set of unrelated individuals.
       
       Args:
           kinship_matrix: N×N kinship matrix
           threshold: Kinship threshold (0.0884 = 3rd degree)
       
       Returns:
           List of unrelated sample indices
       """
       import networkx as nx
       
       n = kinship_matrix.shape[0]
       
       # Create graph of related pairs
       G = nx.Graph()
       G.add_nodes_from(range(n))
       
       for i in range(n):
           for j in range(i + 1, n):
               if kinship_matrix[i, j] > threshold:
                   G.add_edge(i, j)
       
       # Find maximum independent set (NP-hard, use greedy approximation)
       unrelated = nx.maximal_independent_set(G)
       
       return sorted(unrelated)

**Convert Between Genotype Formats:**

.. code-block:: python

   def convert_genotype_format(genotype, from_format, to_format):
       """
       Convert between different genotype encoding formats.
       
       Args:
           genotype: Genotype data
           from_format: '012', 'nucleotide', 'dosage'
           to_format: '012', 'nucleotide', 'dosage'
       
       Returns:
           Converted genotype
       """
       if from_format == '012' and to_format == 'nucleotide':
           # Requires allele information
           # Example: 0->AA, 1->AG, 2->GG
           raise NotImplementedError("Requires allele information")
       
       elif  from_format == 'dosage' and to_format == '012':
           # Round dosages to nearest integer
           return np.round(genotype).astype(int)
       
       elif from_format == '012' and to_format == 'dosage':
           # Already compatible (0, 1, 2 are valid dosages)
           return genotype.astype(float)
       
       else:
           raise ValueError(f"Conversion from {from_format} to {to_format} not implemented")

**Annotate Variants with Gene Information:**

.. code-block:: python

   def annotate_variants_with_genes(variants_df, gtf_file, window_kb=0):
       """
       Annotate variants with nearest gene information.
       
       Args:
           variants_df: DataFrame with 'chr', 'pos' columns
           gtf_file: Path to GTF annotation file
           window_kb: Window around gene (in kb)
       
       Returns:
           DataFrame with gene annotations
       """
       import pyranges as pr
       
       # Load GTF
       genes = pr.read_gtf(gtf_file)
       genes = genes[genes.Feature == "gene"]
       
       # Convert variants to ranges
       variants_gr = pr.PyRanges(
           chromosomes=variants_df['chr'].values,
           starts=variants_df['pos'].values,
           ends=variants_df['pos'].values + 1
       )
       
       # Find nearest genes
       nearest = variants_gr.nearest(
           genes,
           suffix="_gene",
           how="upstream"
       )
       
       # Add gene information to variants
       variants_annotated = variants_df.copy()
       variants_annotated['gene_name'] = nearest.Name.values
       variants_annotated['gene_distance'] = nearest.Distance.values
       
       return variants_annotated

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
