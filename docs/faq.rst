.. _faq:

Frequently Asked Questions (FAQ)
=================

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
