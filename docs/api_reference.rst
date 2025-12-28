.. _api_reference:

API Reference
=============

**UNDER CONSTRUCTION**

Major API documentation for edge-gwas package.

Quick Function Finder
----------------------

**"I want to..."**

* **Get started with test data**:

  * Download test files: ``download_test_files()``

* **Load genetic data**:

  * PLINK binary (.bed/.bim/.fam): ``load_plink_data()``
  * PLINK2 (.pgen/.pvar/.psam): ``load_pgen_data()``
  * VCF: ``load_vcf_data()``
  * BGEN: ``load_bgen_data()``
  * Phenotype: ``prepare_phenotype_data()``

* **Quality control**:

  * Comprehensive QC: ``filter_genotype_data()``
  * Filter by MAF: ``filter_variants_by_maf()``
  * Filter by missingness: ``filter_variants_by_missing()``
  * Filter by HWE: ``filter_variants_by_hwe()``
  * Filter samples: ``filter_samples_by_call_rate()``
  * Check case/control balance: ``check_case_control_balance()``
  * Validate data: ``validate_and_align_data()``

* **Control for population structure**:

  * Calculate GRM: ``calculate_grm_gcta()``
  * Load GRM: ``load_grm_gcta()``
  * Calculate PCs (standard): ``calculate_pca_plink()``
  * Calculate PCs (related samples): ``calculate_pca_pcair()``
  * Add PCs to phenotype: ``attach_pcs_to_phenotype()``
  * Find related samples: ``identify_related_samples()``
  * Remove related samples: ``filter_related_samples()``

* **Prepare data for analysis**:

  * Split train/test: ``stratified_train_test_split()``
  * Impute missing covariates: ``impute_covariates()``

* **Run EDGE analysis**:

  * Calculate alpha: ``calculate_alpha()``
  * Apply alpha: ``apply_alpha()``
  * Full workflow: ``run_full_analysis()``
  * Cross-validation: ``cross_validated_edge_analysis()``
  * Check failed SNPs: ``get_skipped_snps()``

* **Compare with standard GWAS**:

  * Run additive model: ``additive_gwas()`` / ``standard_gwas()``

* **Visualize results**:

  * Manhattan plot: ``manhattan_plot()``
  * QQ plot: ``qq_plot()``
  * Alpha distribution: ``plot_alpha_distribution()``

* **Save and format results**:

  * Save results: ``save_results()``
  * Load alpha values: ``load_alpha_values()``
  * Format for LocusZoom: ``format_gwas_output_for_locuszoom()``
  * Save for LocusZoom: ``save_for_locuszoom()``
  * Validate LocusZoom format: ``validate_locuszoom_format()``
  * Create summary report: ``create_summary_report()``

* **Statistical utilities**:

  * Calculate genomic inflation: ``calculate_genomic_inflation()``
  * Merge alpha with GWAS: ``merge_alpha_with_gwas()``

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
