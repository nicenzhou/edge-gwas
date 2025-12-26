.. _user_guide:

User Guide
==========

This guide covers essential aspects of using edge-gwas v0.1.1.

Overview
--------

edge-gwas performs GWAS analysis with flexible genetic encoding to detect nonadditive SNP effects.

**Key Features:**

* Two-stage analysis framework (training/test split)
* Flexible encoding for nonadditive effects detection
* Population structure control (PCA, GRM)
* Multiple file format support (PLINK, PGEN, BGEN, VCF)
* Outcome transformations for continuous traits
* Built-in visualization and QC tools

EDGE Methodology
----------------

**EDGE (Efficient Detection of Genetic Effects)** uses a two-stage approach:

1. **Training**: Calculate encoding parameters (α) for each variant
2. **Testing**: Apply α values to perform GWAS

**Alpha (α) Interpretation:**

* **α ≈ 0**: Recessive (only homozygotes affected)
* **α ≈ 0.5**: Additive (heterozygotes intermediate)
* **α ≈ 1**: Dominant (heterozygotes = homozygotes)

This flexible encoding detects effects that standard additive models miss.

Complete Workflow
-----------------

1. Load Data
~~~~~~~~~~~~

.. code-block:: python

   from edge_gwas.utils import load_plink_data, prepare_phenotype_data
   
   # Load genotypes (multiple formats supported)
   geno, info = load_plink_data('data.bed', 'data.bim', 'data.fam')
   
   # Load phenotypes
   pheno = prepare_phenotype_data(
       'pheno.txt',
       outcome_col='disease',
       covariate_cols=['age', 'sex']
   )

2. Quality Control
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from edge_gwas.utils import (
       filter_variants_by_maf,
       filter_variants_by_hwe,
       filter_samples_by_call_rate
   )
   
   # Apply QC filters
   geno = filter_variants_by_maf(geno, min_maf=0.01)
   geno = filter_variants_by_hwe(geno, hwe_threshold=1e-6)
   geno, pheno = filter_samples_by_call_rate(geno, pheno, min_call_rate=0.95)

**Recommended Thresholds:**

* MAF > 1%
* Missingness < 5%
* HWE p > 1e-6
* Sample call rate > 95%

3. Population Structure Control (NEW in v0.1.1)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Option A: PCA (for unrelated samples):**

.. code-block:: python

   from edge_gwas.utils import calculate_pca_plink, attach_pcs_to_phenotype
   
   pca_df = calculate_pca_plink('data', n_pcs=10)
   pheno = attach_pcs_to_phenotype(pheno, pca_df, n_pcs=10)

**Option B: GRM (for related samples):**

.. code-block:: python

   from edge_gwas.utils import calculate_grm_gcta, load_grm_gcta
   
   grm_prefix = calculate_grm_gcta('data')
   grm_matrix, grm_ids = load_grm_gcta(grm_prefix)

4. Split Data
~~~~~~~~~~~~~

.. code-block:: python

   from edge_gwas.utils import stratified_train_test_split
   
   train_g, test_g, train_p, test_p = stratified_train_test_split(
       geno, pheno, 'disease', test_size=0.5, random_state=42
   )

5. Run EDGE Analysis
~~~~~~~~~~~~~~~~~~~~

**Basic Analysis:**

.. code-block:: python

   from edge_gwas import EDGEAnalysis
   from edge_gwas.utils import get_pc_covariate_list
   
   edge = EDGEAnalysis(outcome_type='binary', n_jobs=-1)
   
   covariates = ['age', 'sex'] + get_pc_covariate_list(10)
   
   alpha_df, gwas_df = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       outcome='disease',
       covariates=covariates
   )

**With GRM (NEW in v0.1.1):**

.. code-block:: python

   alpha_df, gwas_df = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       outcome='disease',
       covariates=covariates,
       grm_matrix=grm_matrix,
       grm_sample_ids=grm_ids
   )

**Continuous Outcomes with Transformation (NEW in v0.1.1):**

.. code-block:: python

   edge = EDGEAnalysis(
       outcome_type='continuous',
       outcome_transform='rank_inverse_normal'
   )
   
   alpha_df, gwas_df = edge.run_full_analysis(...)

6. Interpret Results
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Significant hits
   sig = gwas_df[gwas_df['pval'] < 5e-8]
   print(f"Significant variants: {len(sig)}")
   
   # Check genomic inflation
   from edge_gwas.utils import calculate_genomic_inflation
   lambda_gc = calculate_genomic_inflation(gwas_df['pval'])
   print(f"Lambda GC: {lambda_gc:.3f}")  # Should be 0.95-1.05
   
   # Alpha distribution
   print(f"Mean alpha: {alpha_df['alpha_value'].mean():.3f}")

**Significance Thresholds:**

* **p < 5×10⁻⁸**: Genome-wide significant
* **p < 1×10⁻⁵**: Suggestive
* **λ = 0.95-1.05**: Good control for population structure

7. Visualize
~~~~~~~~~~~~

.. code-block:: python

   from edge_gwas.visualize import manhattan_plot, qq_plot
   
   manhattan_plot(gwas_df, 'manhattan.png')
   lambda_gc = qq_plot(gwas_df, 'qq.png')

File Formats
------------

**Supported Genotype Formats:**

* **PLINK**: ``.bed/.bim/.fam`` → ``load_plink_data()``
* **PLINK2**: ``.pgen/.pvar/.psam`` → ``load_pgen_data()``
* **BGEN**: ``.bgen`` → ``load_bgen_data()``
* **VCF**: ``.vcf/.vcf.gz`` → ``load_vcf_data()``

**Phenotype Format:**

Tab-delimited text with columns:

* Sample ID (IID)
* Outcome variable
* Covariates

Example:

.. code-block:: text

   IID       disease  age  sex
   SAMPLE1   1        45   0
   SAMPLE2   0        38   1
   SAMPLE3   1        52   0

Best Practices
--------------

Quality Control
~~~~~~~~~~~~~~~

1. **Variant QC**: MAF > 0.01, missingness < 5%, HWE p > 1e-6
2. **Sample QC**: Call rate > 95%, check for relatedness
3. **Population Structure**: Always include ≥10 PCs as covariates
4. **Genomic Inflation**: λ should be 0.95-1.05

Covariates
~~~~~~~~~~

**Always include:**

* Age, sex, genotyping batch
* 10-20 principal components

**For biobank data:**

* Study center, assessment date
* 20+ PCs if diverse ancestry

Data Splitting
~~~~~~~~~~~~~~

* **50/50 split** recommended for balanced power
* **Stratify** binary outcomes to maintain case/control balance
* **Same random seed** for reproducibility

Performance
~~~~~~~~~~~

* **Parallel processing**: Set ``n_jobs=-1`` to use all cores
* **Large datasets**: Process chromosomes separately
* **Memory**: Use binary formats (PLINK, PGEN, BGEN)

Advanced Features (v0.1.1)
--------------------------

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
   
   # Compare p-values
   comparison = gwas_df.merge(add_results, on='variant_id', suffixes=('_edge', '_add'))

Population Structure Options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - Scenario
     - Method
     - Command
   * - Unrelated samples
     - Standard PCA
     - ``calculate_pca_plink()``
   * - Related samples
     - PC-AiR
     - ``calculate_pca_pcair()``
   * - Biobank data
     - GRM + PC-AiR
     - ``calculate_grm_gcta()`` + ``calculate_pca_pcair()``

Outcome Transformations
~~~~~~~~~~~~~~~~~~~~~~~

For non-normal continuous traits:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Transformation
     - When to Use
   * - ``'log'``
     - Right-skewed data (e.g., gene expression)
   * - ``'log10'``
     - Right-skewed with large range
   * - ``'inverse_normal'``
     - Parametric normalization
   * - ``'rank_inverse_normal'``
     - Robust to outliers (recommended)

Troubleshooting
---------------

Common Issues
~~~~~~~~~~~~~

**High genomic inflation (λ > 1.05)**

* Add more PCs (try 15-20 instead of 10)
* Use GRM for related samples
* Check for batch effects

**"No common samples found"**

* Ensure all IDs are strings: ``df.index = df.index.astype(str)``
* Check ID column names match

**Memory errors**

* Process chromosomes separately
* Reduce ``n_jobs``
* Use sparse GRM: ``method='grm-sparse'``

**Convergence warnings**

* Increase ``max_iter``
* Check for rare variants (filter MAF > 0.05)
* View skipped SNPs: ``edge.get_skipped_snps()``

Quick Reference
---------------

**Import commonly used functions:**

.. code-block:: python

   from edge_gwas import EDGEAnalysis
   from edge_gwas.utils import (
       load_plink_data, prepare_phenotype_data,
       filter_variants_by_maf, filter_variants_by_hwe,
       calculate_pca_plink, attach_pcs_to_phenotype,
       get_pc_covariate_list, stratified_train_test_split,
       calculate_grm_gcta, load_grm_gcta
   )
   from edge_gwas.visualize import manhattan_plot, qq_plot

**Typical workflow:**

.. code-block:: python

   # 1. Load → 2. QC → 3. PCA → 4. Split → 5. EDGE → 6. Visualize
   geno, info = load_plink_data(...)
   geno = filter_variants_by_maf(geno, 0.01)
   pca_df = calculate_pca_plink('data', 10)
   train_g, test_g, train_p, test_p = stratified_train_test_split(...)
   alpha_df, gwas_df = edge.run_full_analysis(...)
   manhattan_plot(gwas_df, 'results.png')

See Also
--------

**Documentation:**

* :ref:`index` - Index page of the documentation
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
