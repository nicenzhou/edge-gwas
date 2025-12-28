.. _quickstart:

Quick Start Guide
=================

This guide demonstrates a complete EDGE GWAS analysis in minutes using v0.1.1 features.

Overview
--------

edge-gwas performs GWAS analysis with flexible genetic encoding to detect nonadditive SNP effects.

**Key Features:**

* Two-stage analysis (training/test split)
* Flexible encoding for nonadditive effects
* Population structure control (PCA, GRM)
* Multiple file formats (PLINK, PGEN, BGEN, VCF)
* Outcome transformations
* Built-in visualization

EDGE Methodology
----------------

**Two-stage approach:**

1. **Training**: Calculate encoding parameters (α) for each variant
2. **Testing**: Apply α values to perform GWAS

**Alpha (α) Interpretation:**

* **α ≈ 0**: Recessive (only homozygotes affected)
* **α ≈ 0.5**: Additive (heterozygotes intermediate)
* **α ≈ 1**: Dominant (heterozygotes = homozygotes)

Complete Workflow Example
--------------------------

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

Step-by-Step Guide
------------------

1. Load Data
~~~~~~~~~~~~

**Genotype data (multiple formats):**

.. code-block:: python

   # PLINK format
   geno, info = load_plink_data('data.bed', 'data.bim', 'data.fam')
   
   # PLINK2 format
   from edge_gwas.utils import load_pgen_data
   geno, info = load_pgen_data('data.pgen', 'data.pvar', 'data.psam')
   
   # BGEN format
   from edge_gwas.utils import load_bgen_data
   geno, info = load_bgen_data('data.bgen', 'data.sample')
   
   # VCF format
   from edge_gwas.utils import load_vcf_data
   geno, info = load_vcf_data('data.vcf.gz')

**Phenotype data:**

.. code-block:: python

   pheno = prepare_phenotype_data(
       'pheno.txt',
       outcome_col='disease',
       covariate_cols=['age', 'sex'],
       sample_id_col='IID'
   )

**Phenotype file format (tab-delimited):**

.. code-block:: text

   IID       disease  age  sex
   SAMPLE1   1        45   0
   SAMPLE2   0        38   1
   SAMPLE3   1        52   0

2. Quality Control
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from edge_gwas.utils import filter_genotype_data
   
   # All-in-one QC
   geno_qc, pheno_qc = filter_genotype_data(
       geno, pheno,
       min_maf=0.01,                    # MAF > 1%
       max_missing_per_variant=0.05,   # < 5% missing
       min_call_rate_per_sample=0.95   # > 95% call rate
   )

**Or apply filters individually:**

.. code-block:: python

   from edge_gwas.utils import (
       filter_variants_by_maf,
       filter_variants_by_hwe,
       filter_samples_by_call_rate
   )
   
   geno = filter_variants_by_maf(geno, min_maf=0.01)
   geno = filter_variants_by_hwe(geno, hwe_threshold=1e-6)
   geno, pheno = filter_samples_by_call_rate(geno, pheno, min_call_rate=0.95)

**QC Thresholds:**

* MAF > 1%
* Missingness < 5%
* HWE p > 1e-6
* Sample call rate > 95%

3. Population Structure Control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Option A: PCA (unrelated samples):**

.. code-block:: python

   pca_df = calculate_pca_plink('data', n_pcs=10)
   pheno = attach_pcs_to_phenotype(pheno, pca_df, n_pcs=10)

**Option B: GRM (related samples):**

.. code-block:: python

   grm_prefix = calculate_grm_gcta('data')
   grm_matrix, grm_ids = load_grm_gcta(grm_prefix)

**Option C: PC-AiR (related samples):**

.. code-block:: python

   from edge_gwas.utils import calculate_pca_pcair
   pca_df = calculate_pca_pcair('data', n_pcs=10)
   pheno = attach_pcs_to_phenotype(pheno, pca_df, n_pcs=10)

4. Split Data
~~~~~~~~~~~~~

.. code-block:: python

   train_g, test_g, train_p, test_p = stratified_train_test_split(
       geno, pheno,
       outcome_col='disease',
       test_size=0.5,          # 50/50 split
       random_state=42         # Reproducibility
   )

5. Run EDGE Analysis
~~~~~~~~~~~~~~~~~~~~

**Binary outcome:**

.. code-block:: python

   edge = EDGEAnalysis(outcome_type='binary', n_jobs=-1)
   
   covariates = ['age', 'sex'] + get_pc_covariate_list(10)
   
   alpha_df, gwas_df = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       outcome='disease',
       covariates=covariates,
       variant_info=info
   )

**With GRM:**

.. code-block:: python

   alpha_df, gwas_df = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       outcome='disease',
       covariates=covariates,
       variant_info=info,
       grm_matrix=grm_matrix,
       grm_sample_ids=grm_ids
   )

**Continuous outcome with transformation:**

.. code-block:: python

   edge = EDGEAnalysis(
       outcome_type='continuous',
       outcome_transform='rank_inverse_normal'
   )
   
   alpha_df, gwas_df = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       outcome='bmi',
       covariates=covariates
   )

6. Visualize Results
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from edge_gwas.visualize import manhattan_plot, qq_plot, plot_alpha_distribution
   
   manhattan_plot(gwas_df, 'manhattan.png')
   lambda_gc = qq_plot(gwas_df, 'qq.png')
   plot_alpha_distribution(alpha_df, 'alpha_dist.png')

7. Interpret Results
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Significant hits
   sig = gwas_df[gwas_df['pval'] < 5e-8]
   print(f"Genome-wide significant: {len(sig)}")
   
   # Genomic inflation
   from edge_gwas.utils import calculate_genomic_inflation
   lambda_gc = calculate_genomic_inflation(gwas_df['pval'])
   print(f"Lambda GC: {lambda_gc:.3f}")  # Should be 0.95-1.05
   
   # Alpha statistics
   print(f"Mean alpha: {alpha_df['alpha_value'].mean():.3f}")
   print(f"Median alpha: {alpha_df['alpha_value'].median():.3f}")

**Significance Thresholds:**

* **p < 5×10⁻⁸**: Genome-wide significant
* **p < 1×10⁻⁵**: Suggestive
* **λ = 0.95-1.05**: Good population structure control

Output Files
------------

When using ``output_prefix='results/edge'``:

* ``results/edge_alpha.txt``: Alpha values for each variant
* ``results/edge_gwas.txt``: GWAS results with p-values
* ``results/edge_skipped.txt``: Variants that failed convergence

Best Practices
--------------

Quality Control
~~~~~~~~~~~~~~~

* MAF > 1%
* Missingness < 5%
* HWE p > 1e-6 (controls)
* Sample call rate > 95%

Covariates
~~~~~~~~~~

**Always include:**

* Age, sex, batch
* 10-20 principal components

**For biobank data:**

* Study center, date
* 20+ PCs if diverse ancestry

Data Splitting
~~~~~~~~~~~~~~

* 50/50 split recommended
* Stratify binary outcomes
* Set random seed

Performance
~~~~~~~~~~~

* Use ``n_jobs=-1`` for all cores
* Process chromosomes separately for large data
* Use binary formats (PLINK, PGEN, BGEN)

Advanced Features
-----------------

Cross-Validation
~~~~~~~~~~~~~~~~

.. code-block:: python

   from edge_gwas.utils import cross_validated_edge_analysis
   
   avg_alpha, meta_gwas, all_alpha, all_gwas = cross_validated_edge_analysis(
       geno, pheno, 'disease', covariates,
       outcome_type='binary', n_folds=5
   )

Compare with Additive GWAS
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from edge_gwas.utils import additive_gwas
   
   add_results = additive_gwas(test_g, test_p, 'disease', covariates, 'binary')
   
   comparison = gwas_df.merge(add_results, on='variant_id', 
                              suffixes=('_edge', '_add'))

Outcome Transformations
~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Transformation
     - When to Use
   * - ``'log'``
     - Right-skewed data (Y > 0)
   * - ``'rank_inverse_normal'``
     - Robust to outliers (recommended)
   * - ``'inverse_normal'``
     - Parametric normalization

Troubleshooting
---------------

**High λ (>1.05):** Add more PCs or use GRM

**No common samples:** Ensure IDs are strings

**Memory errors:** Process chromosomes separately, reduce ``n_jobs``

**Convergence warnings:** Increase ``max_iter``, filter rare variants

Quick Reference
---------------

.. code-block:: python

   # Minimal workflow
   from edge_gwas import EDGEAnalysis
   from edge_gwas.utils import *
   
   # Load
   geno, info = load_plink_data('data.bed', 'data.bim', 'data.fam')
   pheno = prepare_phenotype_data('pheno.txt', 'disease', ['age', 'sex'])
   
   # QC
   geno, pheno = filter_genotype_data(geno, pheno, min_maf=0.01)
   
   # PCA
   pca_df = calculate_pca_plink('data', 10)
   pheno = attach_pcs_to_phenotype(pheno, pca_df, 10)
   
   # Split
   train_g, test_g, train_p, test_p = stratified_train_test_split(
       geno, pheno, 'disease', test_size=0.5
   )
   
   # EDGE
   edge = EDGEAnalysis(outcome_type='binary')
   alpha_df, gwas_df = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       'disease', ['age', 'sex'] + get_pc_covariate_list(10)
   )
   
   # Visualize
   from edge_gwas.visualize import manhattan_plot, qq_plot
   manhattan_plot(gwas_df, 'manhattan.png')
   qq_plot(gwas_df, 'qq.png')

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

*Last updated: 2025-12-28 for edge-gwas v0.1.1*

*For questions or issues, visit:* https://github.com/nicenzhou/edge-gwas/issues
