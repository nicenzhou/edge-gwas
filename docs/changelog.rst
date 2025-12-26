.. _changelog:

Changelog
=========

All notable changes to edge-gwas will be documented in this file.

The format is based on `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_,
and this project adheres to `Semantic Versioning <https://semver.org/spec/v2.0.0.html>`_.

Version 0.1.1 (2025-12-25) - Current Release
--------------------------------------------

**Status:** ⚠️ Public Testing Phase

Added
^^^^^

**Population Structure Control:**

* ``load_grm_gcta()`` - Load GRM calculated by GCTA
* ``calculate_grm_gcta()`` - Calculate GRM using GCTA
* GRM integration in ``calculate_alpha()`` and ``apply_alpha()`` for mixed model analysis
* ``identify_related_samples()`` - Identify related sample pairs from GRM
* ``filter_related_samples()`` - Remove related samples

**Principal Component Analysis:**

* ``calculate_pca_plink()`` - PCA using PLINK2 with LD pruning
* ``calculate_pca_pcair()`` - PC-AiR for related samples
* ``calculate_pca_sklearn()`` - Basic PCA using scikit-learn
* ``attach_pcs_to_phenotype()`` - Merge PCs with phenotype data
* ``get_pc_covariate_list()`` - Helper for generating PC covariate names

**Outcome Transformations:**

* ``outcome_transform`` parameter in ``EDGEAnalysis`` for continuous outcomes
* Support for 'log', 'log10', 'inverse_normal', and 'rank_inverse_normal' transformations

**Additional Data Formats:**

* ``load_pgen_data()`` - Load PLINK2 .pgen/.pvar/.psam files
* ``load_bgen_data()`` - Load BGEN format with dosages
* ``load_vcf_data()`` - Load VCF/VCF.gz files

**Tool Installation:**

* ``edge-gwas-install-tools`` - Interactive installer for PLINK2, GCTA, and R packages
* ``edge-gwas-check-tools`` - Verify external tool installations
* Support for Linux and macOS (including ARM64)

**Quality Control:**

* ``filter_samples_by_call_rate()`` - Filter samples by genotype call rate
* ``calculate_hwe_pvalues()`` - Calculate Hardy-Weinberg Equilibrium p-values
* ``filter_variants_by_hwe()`` - Filter variants by HWE p-value
* ``check_case_control_balance()`` - Check case/control ratio

**Analysis Methods:**

* ``cross_validated_edge_analysis()`` - K-fold cross-validation for EDGE
* ``additive_gwas()`` - Standard additive model for comparison

Changed
^^^^^^^

**Data Handling:**

* **Breaking Change:** Replaced Koalas DataFrames with pandas DataFrames throughout codebase
* Updated ``load_participant_data()`` to use pandas instead of ``participant.retrieve_fields().to_koalas()``
* All PCA functions now return DataFrames with 'IID' column and IID as index
* Consistent sample ID handling (all IDs converted to strings for matching)

**Core Analysis:**

* ``EDGEAnalysis`` class now supports GRM as optional input in all methods
* ``calculate_alpha()`` and ``apply_alpha()`` accept ``grm_matrix`` and ``grm_sample_ids`` parameters
* Linear mixed model implementation for continuous outcomes with GRM
* Logistic mixed model for binary outcomes with GRM

**Dependencies:**

* Removed dependency on ``databricks.koalas``
* All operations now use native pandas (>= 1.2.0)
* Added optional dependencies: ``pgenlib``, ``bgen-reader``, ``cyvcf2``
* External tool requirements: PLINK2, GCTA, R (all optional)

Migration from v0.1.0
^^^^^^^^^^^^^^^^^^^^^

**Breaking Changes:**

* Code using ``.to_koalas()`` must be updated to use pandas DataFrames
* Remove any ``import databricks.koalas`` statements
* PCA functions now return DataFrames indexed by IID

**Migration Example:**

.. code-block:: python

   # Old code (v0.1.0)
   import databricks.koalas as ks
   participant_df = participant.retrieve_fields(fields=fields, engine=dxdata.connect()).to_koalas()
   
   # New code (v0.1.1)
   import pandas as pd
   participant_df = participant.retrieve_fields(fields=fields, engine=dxdata.connect())
   
   # New features in v0.1.1
   from edge_gwas.utils import calculate_grm_gcta, load_grm_gcta, calculate_pca_plink
   
   # Calculate GRM for population structure control
   grm_prefix = calculate_grm_gcta('genotypes', maf_threshold=0.01)
   grm_matrix, grm_ids = load_grm_gcta(grm_prefix)
   
   # Calculate PCA
   pca_df = calculate_pca_plink('genotypes', n_pcs=10)
   
   # Use outcome transformation and GRM in analysis
   edge = EDGEAnalysis(outcome_type='continuous', outcome_transform='rank_inverse_normal')
   alpha_df, gwas_df = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       outcome='trait', covariates=['age', 'sex', 'PC1', 'PC2', 'PC3'],
       grm_matrix=grm_matrix, grm_sample_ids=grm_ids
   )

Bug Fixes
^^^^^^^^^

* Resolved compatibility issues with recent pandas versions
* Improved memory efficiency in data loading operations
* Fixed GRM sample alignment when samples differ between genotype and GRM
* Fixed ID type inconsistencies (all IDs now converted to strings)
* Improved convergence in mixed model fitting

Version 0.1.0 (2025-12-24)
--------------------------

**Status:** Superseded by v0.1.1

Initial Release
~~~~~~~~~~~~~~~

This is the first packaged release of EDGE-GWAS, transitioning from standalone scripts (v0.0.0) 
to a proper Python package.

Added
^^^^^

**Core Functionality:**

* ``EDGEAnalysis`` class for two-stage GWAS analysis
* ``calculate_alpha()`` - Calculate encoding parameters from training data
* ``apply_alpha()`` - Apply encoding to test data for GWAS
* ``run_full_analysis()`` - Complete end-to-end workflow
* Support for binary outcomes (logistic regression)
* Support for quantitative outcomes (linear regression)

**Data Handling:**

* ``load_plink_data()`` - Load PLINK binary format (.bed/.bim/.fam)
* ``prepare_phenotype_data()`` - Load and prepare phenotype files
* ``stratified_train_test_split()`` - Train/test data splitting with stratification
* ``filter_variants_by_maf()`` - Minor allele frequency filtering
* ``filter_variants_by_missing()`` - Missingness filtering

**Statistical Functions:**

* ``calculate_genomic_inflation()`` - Calculate genomic inflation factor (λ)
* Two-stage analysis framework (calculate alpha → apply alpha)
* Parallel processing support with configurable CPU cores
* Convergence monitoring and skipped SNP tracking

**Visualization:**

* ``manhattan_plot()`` - Create Manhattan plots with customization
* ``qq_plot()`` - Create QQ plots with lambda calculation
* ``plot_alpha_distribution()`` - Visualize alpha value distribution

**Documentation:**

* Complete Sphinx documentation with Read the Docs integration
* Installation guide, quick start guide, and API reference

Dependencies
^^^^^^^^^^^^

* numpy >= 1.19.0
* pandas >= 1.2.0
* scipy >= 1.6.0
* statsmodels >= 0.12.0
* scikit-learn >= 0.24.0
* matplotlib >= 3.3.0
* pandas-plink >= 2.0.0
* databricks-koalas (deprecated in v0.1.1)

---

Version 0.0.0 (2024-04-02) - Legacy
------------------------------------

**Status:** 🔒 Deprecated - No Maintenance

Original EDGE implementation as standalone Python scripts.

.. warning::
   Version 0.0.0 is **deprecated** and no longer maintained.
   Users should migrate to v0.1.1 or later.

---

Version History Summary
-----------------------

.. list-table::
   :header-rows: 1
   :widths: 10 15 15 60

   * - Version
     - Date
     - Status
     - Key Features
   * - 0.1.1
     - 2025-12-25
     - **Current**
     - GRM support, PCA methods, outcome transformations, pandas migration
   * - 0.1.0
     - 2025-12-24
     - Superseded
     - First packaged release, core EDGE functionality
   * - 0.0.0
     - 2024-04-02
     - Deprecated
     - Original standalone scripts

---

See Also
--------

* :ref:`installation` - Installation instructions
* :ref:`quickstart` - Getting started guide
* GitHub Releases - https://github.com/YOUR-USERNAME/edge-gwas/releases
* Issue Tracker - https://github.com/YOUR-USERNAME/edge-gwas/issues
