.. _api_reference:

=============
API Reference
=============

Quick Function Finder
======================

**"I want to..."**

* **Get started with test data**:

  * Download test files: :func:`download_test_files`

* **Load genetic data**:

  * PLINK binary (.bed/.bim/.fam): :func:`load_plink_data`
  * PLINK2 (.pgen/.pvar/.psam): :func:`load_pgen_data`
  * VCF: :func:`load_vcf_data`
  * BGEN: :func:`load_bgen_data`
  * Phenotype: :func:`prepare_phenotype_data`

* **Validate and check data quality**:

  * Validate genotype encoding: :func:`validate_genotype_df`
  * Fix encoding issues: :func:`validate_and_fix_encoding`
  * Validate phenotype data: :func:`validate_phenotype_df`
  * Validate and align datasets: :func:`validate_and_align_data`

* **Quality control**:

  * Comprehensive QC: :func:`filter_genotype_data`
  * Filter by MAF: :func:`filter_variants_by_maf`
  * Filter by missingness: :func:`filter_variants_by_missing`
  * Filter by HWE: :func:`filter_variants_by_hwe`
  * Calculate HWE p-values: :func:`calculate_hwe_pvalues`
  * Filter samples: :func:`filter_samples_by_call_rate`
  * Check case/control balance: :func:`check_case_control_balance`

* **Control for population structure**:

  * Calculate GRM: :func:`calculate_grm_gcta`
  * Load GRM: :func:`load_grm_gcta`
  * Calculate PCs (basic): :func:`calculate_pca_sklearn`
  * Calculate PCs (with PLINK2): :func:`calculate_pca_plink`
  * Calculate PCs (for related samples): :func:`calculate_pca_pcair`
  * Add PCs to phenotype: :func:`attach_pcs_to_phenotype`
  * Get PC covariate names: :func:`get_pc_covariate_list`
  * Find related samples: :func:`identify_related_samples`
  * Remove related samples: :func:`filter_related_samples`

* **Prepare data for analysis**:

  * Split train/test: :func:`stratified_train_test_split`
  * Impute missing covariates: :func:`impute_covariates`

* **Run EDGE analysis**:

  * Initialize EDGE: :class:`EDGEAnalysis`
  * Calculate alpha: :meth:`EDGEAnalysis.calculate_alpha`
  * Apply alpha: :meth:`EDGEAnalysis.apply_alpha`
  * Full workflow: :meth:`EDGEAnalysis.run_full_analysis`
  * Cross-validation: :func:`cross_validated_edge_analysis`
  * Check failed SNPs: :meth:`EDGEAnalysis.get_skipped_snps`

* **Compare with standard GWAS**:

  * Run additive model: :func:`standard_gwas` / :func:`additive_gwas`

* **Visualize results**:

  * Manhattan plot: :func:`manhattan_plot`
  * QQ plot: :func:`qq_plot`
  * Alpha distribution: :func:`plot_alpha_distribution`

* **Save and format results**:

  * Save results: :func:`save_results`
  * Load alpha values: :func:`load_alpha_values`
  * Format for publication: :func:`format_gwas_output_for_locuszoom`
  * Save for LocusZoom: :func:`save_for_locuszoom`
  * Validate LocusZoom format: :func:`validate_locuszoom_format`
  * Create summary report: :func:`create_summary_report`

----
.. _core_analysis:

Core Analysis
=============

.. class:: EDGEAnalysis(outcome_type='binary', outcome_transform=None, ols_method='bfgs', n_jobs=-1, max_iter=1000, verbose=True)

   Main class for EDGE GWAS analysis.

   :param outcome_type: Type of outcome - 'binary' for logistic regression or 'continuous' for linear regression
   :type outcome_type: str
   :param outcome_transform: Transformation for continuous outcomes. Options: None, 'log', 'log10', 'inverse_normal', 'rank_inverse_normal'
   :type outcome_transform: str, optional
   :param ols_method: Optimization method for OLS regression. Options: 'newton', 'bfgs', 'lbfgs', 'nm', 'cg', 'ncg', 'powell', 'basinhopping'
   :type ols_method: str
   :param n_jobs: Number of parallel jobs (-1 uses all available cores)
   :type n_jobs: int
   :param max_iter: Maximum iterations for model convergence
   :type max_iter: int
   :param verbose: Print progress information
   :type verbose: bool

   .. method:: calculate_alpha(genotype_data, phenotype_df, outcome, covariates, variant_info=None, grm_matrix=None, grm_sample_ids=None, mean_centered=False, use_fast_approximation=True)

      Calculate EDGE alpha values from training data.

      :param genotype_data: Genotype data with samples as index and variants as columns (0/1/2 encoding)
      :type genotype_data: pd.DataFrame
      :param phenotype_df: Phenotype data with sample IDs as index
      :type phenotype_df: pd.DataFrame
      :param outcome: Name of outcome variable in phenotype_df
      :type outcome: str
      :param covariates: List of covariate names in phenotype_df
      :type covariates: list
      :param variant_info: Optional variant information with variant_id as index
      :type variant_info: pd.DataFrame, optional
      :param grm_matrix: Optional GRM matrix from GCTA (for population structure control)
      :type grm_matrix: np.ndarray, optional
      :param grm_sample_ids: DataFrame with FID, IID, and sample_id corresponding to GRM rows
      :type grm_sample_ids: pd.DataFrame, optional
      :param mean_centered: If True, use mean-centered codominant model without intercept
      :type mean_centered: bool
      :param use_fast_approximation: If True, use faster approximation for GRM-based binary models
      :type use_fast_approximation: bool
      :returns: DataFrame with alpha values for each variant
      :rtype: pd.DataFrame

   .. method:: apply_alpha(genotype_data, phenotype_df, outcome, covariates, alpha_values=None, grm_matrix=None, grm_sample_ids=None, variant_info=None, use_fast_approximation=True)

      Apply EDGE alpha values to test data and perform GWAS.

      :param genotype_data: Genotype data with samples as index and variants as columns (0/1/2 encoding)
      :type genotype_data: pd.DataFrame
      :param phenotype_df: Phenotype data with sample IDs as index
      :type phenotype_df: pd.DataFrame
      :param outcome: Name of outcome variable in phenotype_df
      :type outcome: str
      :param covariates: List of covariate names in phenotype_df
      :type covariates: list
      :param alpha_values: DataFrame with alpha values (from calculate_alpha). If None, uses self.alpha_values
      :type alpha_values: pd.DataFrame, optional
      :param grm_matrix: Optional GRM matrix from GCTA
      :type grm_matrix: np.ndarray, optional
      :param grm_sample_ids: DataFrame with FID, IID, and sample_id corresponding to GRM rows
      :type grm_sample_ids: pd.DataFrame, optional
      :param variant_info: Optional variant information DataFrame
      :type variant_info: pd.DataFrame, optional
      :param use_fast_approximation: If True, use faster approximation for GRM-based binary models
      :type use_fast_approximation: bool
      :returns: DataFrame with GWAS results
      :rtype: pd.DataFrame

   .. method:: run_full_analysis(train_genotype, train_phenotype, test_genotype, test_phenotype, outcome, covariates, variant_info=None, grm_matrix=None, grm_sample_ids=None, mean_centered=False, use_fast_approximation=True, output_prefix=None)

      Run complete EDGE analysis: calculate alpha on training data, apply alpha on test data.

      :param train_genotype: Training genotype data
      :type train_genotype: pd.DataFrame
      :param train_phenotype: Training phenotype data
      :type train_phenotype: pd.DataFrame
      :param test_genotype: Test genotype data
      :type test_genotype: pd.DataFrame
      :param test_phenotype: Test phenotype data
      :type test_phenotype: pd.DataFrame
      :param outcome: Name of outcome variable
      :type outcome: str
      :param covariates: List of covariate names
      :type covariates: list
      :param variant_info: Optional variant information DataFrame
      :type variant_info: pd.DataFrame, optional
      :param grm_matrix: Optional GRM matrix from GCTA
      :type grm_matrix: np.ndarray, optional
      :param grm_sample_ids: Optional sample IDs for GRM
      :type grm_sample_ids: pd.DataFrame, optional
      :param mean_centered: If True, use mean-centered model without intercept
      :type mean_centered: bool
      :param use_fast_approximation: If True, use faster approximation for GRM-based binary models
      :type use_fast_approximation: bool
      :param output_prefix: Optional prefix for output files
      :type output_prefix: str, optional
      :returns: Tuple of (alpha_df, gwas_df)
      :rtype: tuple

   .. method:: get_skipped_snps()

      Get list of SNPs that were skipped due to convergence issues.

      :returns: List of skipped SNP IDs
      :rtype: list

----
.. _data_loading:

Data Loading
============

.. function:: load_plink_data(bed_file, bim_file, fam_file, minor_allele_as_alt=True, verbose=True)

   Load PLINK binary format data (.bed/.bim/.fam).

   :param bed_file: Path to .bed file
   :type bed_file: str
   :param bim_file: Path to .bim file
   :type bim_file: str
   :param fam_file: Path to .fam file
   :type fam_file: str
   :param minor_allele_as_alt: If True, ensure minor allele is coded as ALT (2)
   :type minor_allele_as_alt: bool
   :param verbose: Print loading information
   :type verbose: bool
   :returns: Tuple of (genotype_df, variant_info_df)
   :rtype: tuple

.. function:: load_pgen_data(pgen_file, pvar_file, psam_file, minor_allele_as_alt=True, verbose=True)

   Load PLINK 2 binary format data (.pgen/.pvar/.psam).

   :param pgen_file: Path to .pgen file
   :type pgen_file: str
   :param pvar_file: Path to .pvar file
   :type pvar_file: str
   :param psam_file: Path to .psam file
   :type psam_file: str
   :param minor_allele_as_alt: If True, ensure minor allele is coded as ALT (2)
   :type minor_allele_as_alt: bool
   :param verbose: Print loading information
   :type verbose: bool
   :returns: Tuple of (genotype_df, variant_info_df)
   :rtype: tuple

   .. note::
      Requires pgenlib package: ``pip install pgenlib``

.. function:: load_vcf_data(vcf_file, dosage=True, minor_allele_as_alt=True, verbose=True)

   Load VCF format data.

   :param vcf_file: Path to .vcf or .vcf.gz file
   :type vcf_file: str
   :param dosage: If True, use dosages (DS field); if False, use hard calls (GT field)
   :type dosage: bool
   :param minor_allele_as_alt: If True, ensure minor allele is coded as ALT (2)
   :type minor_allele_as_alt: bool
   :param verbose: Print loading information
   :type verbose: bool
   :returns: Tuple of (genotype_df, variant_info_df)
   :rtype: tuple

   .. note::
      Requires cyvcf2 package: ``pip install cyvcf2``

.. function:: load_bgen_data(bgen_file, sample_file=None, minor_allele_as_alt=True, verbose=True)

   Load BGEN format data.

   :param bgen_file: Path to .bgen file
   :type bgen_file: str
   :param sample_file: Path to .sample file (optional, can be embedded in BGEN)
   :type sample_file: str, optional
   :param minor_allele_as_alt: If True, ensure minor allele is coded as ALT (2)
   :type minor_allele_as_alt: bool
   :param verbose: Print loading information
   :type verbose: bool
   :returns: Tuple of (genotype_df, variant_info_df) - genotypes are dosages
   :rtype: tuple

   .. note::
      Requires bgen_reader package: ``pip install bgen-reader``

.. function:: prepare_phenotype_data(phenotype_file, outcome_col, covariate_cols, sample_id_col='IID', sep='\\t', log_transform_outcome=False)

   Load and prepare phenotype data.

   :param phenotype_file: Path to phenotype file
   :type phenotype_file: str
   :param outcome_col: Name of outcome column
   :type outcome_col: str
   :param covariate_cols: List of covariate column names
   :type covariate_cols: list
   :param sample_id_col: Name of sample ID column (will become index)
   :type sample_id_col: str
   :param sep: File separator
   :type sep: str
   :param log_transform_outcome: Apply log10(x+1) transformation to outcome
   :type log_transform_outcome: bool
   :returns: DataFrame with sample IDs as index, outcome and covariates as columns
   :rtype: pd.DataFrame

.. function:: download_test_files(output_dir='tests', version='v0.1.2', overwrite=False, verbose=True)

   Download test files from GitHub repository.

   :param output_dir: Directory to save test files
   :type output_dir: str
   :param version: GitHub release version tag
   :type version: str
   :param overwrite: If True, overwrite existing files
   :type overwrite: bool
   :param verbose: Print download progress
   :type verbose: bool
   :returns: Dictionary with download results (downloaded, skipped, failed)
   :rtype: dict

----
.. _data_validation:

Data Validation
===============

.. function:: validate_genotype_df(genotype_df, variant_info_df=None, name='genotype_df', check_encoding=True, verbose=True, return_details=False)

   Validate genotype DataFrame format and encoding.

   :param genotype_df: Genotype DataFrame (samples x variants)
   :type genotype_df: pd.DataFrame
   :param variant_info_df: Optional variant information DataFrame
   :type variant_info_df: pd.DataFrame, optional
   :param name: Name for error messages
   :type name: str
   :param check_encoding: If True, validate encoding (requires variant_info_df)
   :type check_encoding: bool
   :param verbose: Print validation results
   :type verbose: bool
   :param return_details: If True, return (passed, report_df)
   :type return_details: bool
   :returns: None (raises errors if invalid) OR bool (validation passed) OR Tuple[bool, pd.DataFrame] (if return_details=True)
   :rtype: None, bool, or tuple

.. function:: validate_and_fix_encoding(genotype_df, variant_info_df, verbose=True)

   Validate and automatically fix genotype encoding.

   :param genotype_df: Genotype DataFrame
   :type genotype_df: pd.DataFrame
   :param variant_info_df: Variant info DataFrame
   :type variant_info_df: pd.DataFrame
   :param verbose: Print progress
   :type verbose: bool
   :returns: Tuple of (fixed_genotype_df, fixed_variant_info_df, report_df)
   :rtype: tuple

.. function:: validate_phenotype_df(phenotype_df, outcome_col, covariate_cols, name='phenotype_df')

   Validate phenotype DataFrame format.

   :param phenotype_df: Phenotype DataFrame to validate
   :type phenotype_df: pd.DataFrame
   :param outcome_col: Name of outcome column
   :type outcome_col: str
   :param covariate_cols: List of covariate column names
   :type covariate_cols: list
   :param name: Name of the DataFrame for error messages
   :type name: str
   :raises TypeError: If not a pandas DataFrame
   :raises ValueError: If required columns are missing or DataFrame is invalid

.. function:: validate_and_align_data(genotype_df, phenotype_df, outcome_col=None, covariate_cols=None, geno_id_col=None, pheno_id_col=None, keep_only_common=True, verbose=True)

   Validate and align genotype and phenotype data by sample IDs.

   :param genotype_df: Genotype DataFrame (samples x variants)
   :type genotype_df: pd.DataFrame
   :param phenotype_df: Phenotype DataFrame
   :type phenotype_df: pd.DataFrame
   :param outcome_col: Name of outcome column (optional, for validation)
   :type outcome_col: str, optional
   :param covariate_cols: List of covariate columns (optional, for validation)
   :type covariate_cols: list, optional
   :param geno_id_col: Column name for sample IDs in genotype_df (None = use index)
   :type geno_id_col: str, optional
   :param pheno_id_col: Column name for sample IDs in phenotype_df (None = use index)
   :type pheno_id_col: str, optional
   :param keep_only_common: If True, keep only samples present in both datasets
   :type keep_only_common: bool
   :param verbose: Print validation information
   :type verbose: bool
   :returns: Tuple of (aligned_genotype_df, aligned_phenotype_df)
   :rtype: tuple
   :raises ValueError: If no common samples found or if keep_only_common=False and samples don't match

----
.. _quality_control:

Quality Control
===============

.. function:: filter_genotype_data(genotype_df, phenotype_df=None, min_maf=None, max_missing_per_variant=None, min_call_rate_per_sample=None, verbose=True)

   Comprehensive genotype data filtering with multiple QC criteria.

   :param genotype_df: Genotype DataFrame (samples x variants)
   :type genotype_df: pd.DataFrame
   :param phenotype_df: Optional phenotype DataFrame (required if filtering samples)
   :type phenotype_df: pd.DataFrame, optional
   :param min_maf: Minimum minor allele frequency (e.g., 0.01 for 1%). If None, no MAF filtering
   :type min_maf: float, optional
   :param max_missing_per_variant: Maximum missing rate per variant (e.g., 0.1 for 10%). If None, no filtering
   :type max_missing_per_variant: float, optional
   :param min_call_rate_per_sample: Minimum call rate per sample (e.g., 0.95 for 95%). If None, no filtering
   :type min_call_rate_per_sample: float, optional
   :param verbose: Print filtering information
   :type verbose: bool
   :returns: filtered_genotype_df OR (filtered_genotype_df, filtered_phenotype_df)
   :rtype: pd.DataFrame or tuple

.. function:: filter_variants_by_maf(genotype_df, min_maf=0.01, verbose=True)

   Filter variants by minor allele frequency.

   :param genotype_df: Genotype DataFrame (works with both hard calls and dosages)
   :type genotype_df: pd.DataFrame
   :param min_maf: Minimum minor allele frequency
   :type min_maf: float
   :param verbose: Print filtering information
   :type verbose: bool
   :returns: Filtered genotype DataFrame
   :rtype: pd.DataFrame

.. function:: filter_variants_by_missing(genotype_df, max_missing=0.1, verbose=True)

   Filter variants by missing genotype rate.

   :param genotype_df: Genotype DataFrame
   :type genotype_df: pd.DataFrame
   :param max_missing: Maximum proportion of missing genotypes allowed (0-1)
   :type max_missing: float
   :param verbose: Print filtering information
   :type verbose: bool
   :returns: Filtered genotype DataFrame
   :rtype: pd.DataFrame

.. function:: filter_samples_by_call_rate(genotype_df, phenotype_df, min_call_rate=0.95, verbose=True)

   Filter samples by genotype call rate.

   :param genotype_df: Genotype DataFrame (samples as index)
   :type genotype_df: pd.DataFrame
   :param phenotype_df: Phenotype DataFrame (sample IDs as index)
   :type phenotype_df: pd.DataFrame
   :param min_call_rate: Minimum call rate (proportion of non-missing genotypes, 0-1)
   :type min_call_rate: float
   :param verbose: Print filtering information
   :type verbose: bool
   :returns: Tuple of (filtered_genotype_df, filtered_phenotype_df)
   :rtype: tuple

.. function:: calculate_hwe_pvalues(genotype_df, verbose=True)

   Calculate Hardy-Weinberg Equilibrium p-values for each variant.

   :param genotype_df: Genotype DataFrame
   :type genotype_df: pd.DataFrame
   :param verbose: Print calculation information
   :type verbose: bool
   :returns: Series of HWE p-values for each variant
   :rtype: pd.Series

.. function:: filter_variants_by_hwe(genotype_df, hwe_threshold=1e-6, verbose=True)

   Filter variants by Hardy-Weinberg Equilibrium p-value.

   :param genotype_df: Genotype DataFrame
   :type genotype_df: pd.DataFrame
   :param hwe_threshold: Minimum HWE p-value threshold
   :type hwe_threshold: float
   :param verbose: Print filtering information
   :type verbose: bool
   :returns: Filtered genotype DataFrame
   :rtype: pd.DataFrame

.. function:: check_case_control_balance(phenotype_df, outcome_col, verbose=True)

   Check case/control balance in binary outcome.

   :param phenotype_df: Phenotype DataFrame
   :type phenotype_df: pd.DataFrame
   :param outcome_col: Name of outcome column
   :type outcome_col: str
   :param verbose: Print balance information
   :type verbose: bool
   :returns: Dictionary with case_count, control_count, and ratio
   :rtype: dict

----
.. _population_structure:

Population Structure Control
=============================

.. function:: calculate_grm_gcta(plink_prefix, output_prefix=None, maf_threshold=0.01, method='grm', max_threads=1, verbose=True)

   Calculate genetic relationship matrix (GRM) using GCTA.

   :param plink_prefix: Prefix for PLINK binary files (.bed/.bim/.fam)
   :type plink_prefix: str
   :param output_prefix: Prefix for output GRM files (default: temp directory)
   :type output_prefix: str, optional
   :param maf_threshold: MAF threshold for variant filtering
   :type maf_threshold: float
   :param method: GRM calculation method ('grm' for full, 'grm-sparse' for sparse)
   :type method: str
   :param max_threads: Maximum number of threads to use
   :type max_threads: int
   :param verbose: Print progress information
   :type verbose: bool
   :returns: Path to output GRM prefix
   :rtype: str

   .. note::
      Requires GCTA to be installed and available in PATH.
      Download from: https://yanglab.westlake.edu.cn/software/gcta/

.. function:: load_grm_gcta(grm_prefix, verbose=True)

   Load GRM calculated by GCTA.

   :param grm_prefix: Prefix for GRM files (without .grm.bin extension)
   :type grm_prefix: str
   :param verbose: Print loading information
   :type verbose: bool
   :returns: Tuple of (grm_matrix, sample_ids_df)
   :rtype: tuple
   :raises FileNotFoundError: If GRM files are not found

.. function:: calculate_pca_sklearn(genotype_df, n_pcs=10, verbose=True)

   Calculate principal components using scikit-learn (basic PCA without relatedness correction).

   :param genotype_df: Genotype DataFrame (samples x variants)
   :type genotype_df: pd.DataFrame
   :param n_pcs: Number of principal components to calculate
   :type n_pcs: int
   :param verbose: Print progress information
   :type verbose: bool
   :returns: DataFrame with 'IID' as index and PC1, PC2, ..., PCn columns
   :rtype: pd.DataFrame

   .. note::
      This is a basic PCA without correction for relatedness.
      For more robust PCA accounting for relatedness, use :func:`calculate_pca_plink`.

.. function:: calculate_pca_plink(file_prefix, n_pcs=10, file_format='bfile', output_prefix=None, maf_threshold=0.01, ld_window=50, ld_step=5, ld_r2=0.2, approx=False, verbose=True)

   Calculate principal components using PLINK2.

   :param file_prefix: Prefix for input files
   :type file_prefix: str
   :param n_pcs: Number of principal components to calculate
   :type n_pcs: int
   :param file_format: Input file format ('bfile', 'pfile', 'vcf', 'bgen')
   :type file_format: str
   :param output_prefix: Prefix for output files (default: temp directory)
   :type output_prefix: str, optional
   :param maf_threshold: MAF threshold for variant filtering (None to skip)
   :type maf_threshold: float, optional
   :param ld_window: Window size for LD pruning in variant count (None to skip)
   :type ld_window: int, optional
   :param ld_step: Step size for LD pruning in variant count (None to skip)
   :type ld_step: int, optional
   :param ld_r2: R² threshold for LD pruning (None to skip)
   :type ld_r2: float, optional
   :param approx: Use approximate PCA for large cohorts
   :type approx: bool
   :param verbose: Print progress information
   :type verbose: bool
   :returns: DataFrame with IID as index and PC1, PC2, ..., PCn columns
   :rtype: pd.DataFrame

.. function:: calculate_pca_pcair(plink_prefix, n_pcs=10, kinship_matrix=None, divergence_matrix=None, output_prefix=None, kin_threshold=0.0884, div_threshold=-0.0884, maf_threshold=0.01, verbose=True)

   Calculate PC-AiR (Principal Components - Analysis in Related samples).

   :param plink_prefix: Prefix for PLINK binary files
   :type plink_prefix: str
   :param n_pcs: Number of principal components to calculate
   :type n_pcs: int
   :param kinship_matrix: Path to kinship matrix (if None, calculates using GCTA)
   :type kinship_matrix: str, optional
   :param divergence_matrix: Path to divergence matrix (optional)
   :type divergence_matrix: str, optional
   :param output_prefix: Prefix for output files (default: temp directory)
   :type output_prefix: str, optional
   :param kin_threshold: Kinship threshold for defining relatives
   :type kin_threshold: float
   :param div_threshold: Divergence threshold
   :type div_threshold: float
   :param maf_threshold: MAF threshold for variant filtering
   :type maf_threshold: float
   :param verbose: Print progress information
   :type verbose: bool
   :returns: DataFrame with IID as index and PC1, PC2, ..., PCn columns
   :rtype: pd.DataFrame

   .. note::
      Requires R with GENESIS, SNPRelate, and gdsfmt packages installed.

.. function:: attach_pcs_to_phenotype(phenotype_df, pca_df, n_pcs=10, pc_prefix='PC', sample_id_col=None, drop_na=False, verbose=True)

   Attach principal components to phenotype DataFrame.

   :param phenotype_df: Phenotype DataFrame (IID as index or column)
   :type phenotype_df: pd.DataFrame
   :param pca_df: PCA DataFrame with IID as index and PC columns
   :type pca_df: pd.DataFrame
   :param n_pcs: Number of PCs to attach (will use PC1 to PCn)
   :type n_pcs: int
   :param pc_prefix: Prefix for PC column names
   :type pc_prefix: str
   :param sample_id_col: Column name in phenotype_df to use for matching. If None, uses index
   :type sample_id_col: str, optional
   :param drop_na: If True, remove samples with missing PCs after merging
   :type drop_na: bool
   :param verbose: Print information about merging
   :type verbose: bool
   :returns: Phenotype DataFrame with PC columns added
   :rtype: pd.DataFrame
   :raises ValueError: If requested PCs are not available in pca_df

.. function:: get_pc_covariate_list(n_pcs, pc_prefix='PC')

   Generate list of PC covariate names for use in EDGE analysis.

   :param n_pcs: Number of PCs
   :type n_pcs: int
   :param pc_prefix: Prefix for PC column names
   :type pc_prefix: str
   :returns: List of PC column names ['PC1', 'PC2', ..., 'PCn']
   :rtype: list

.. function:: identify_related_samples(grm_matrix, sample_ids, threshold=0.0884, verbose=True)

   Identify pairs of related samples based on GRM threshold.

   :param grm_matrix: n_samples x n_samples GRM matrix
   :type grm_matrix: np.ndarray
   :param sample_ids: DataFrame with sample IDs (from load_grm_gcta)
   :type sample_ids: pd.DataFrame
   :param threshold: Relatedness threshold. Common values: 0.354 (1st degree), 0.177 (2nd degree), 0.0884 (3rd degree)
   :type threshold: float
   :param verbose: Print summary statistics
   :type verbose: bool
   :returns: DataFrame with columns IID1, IID2, kinship (sorted by kinship descending)
   :rtype: pd.DataFrame

.. function:: filter_related_samples(phenotype_df, grm_matrix, sample_ids, threshold=0.0884, method='greedy', sample_id_col=None, verbose=True)

   Filter out related samples to create an unrelated subset.

   :param phenotype_df: Phenotype DataFrame
   :type phenotype_df: pd.DataFrame
   :param grm_matrix: n_samples x n_samples GRM matrix
   :type grm_matrix: np.ndarray
   :param sample_ids: DataFrame with sample IDs (from load_grm_gcta)
   :type sample_ids: pd.DataFrame
   :param threshold: Relatedness threshold
   :type threshold: float
   :param method: Method for selecting unrelated samples ('greedy' or 'random')
   :type method: str
   :param sample_id_col: Column name in phenotype_df for sample IDs. If None, uses index
   :type sample_id_col: str, optional
   :param verbose: Print filtering information
   :type verbose: bool
   :returns: Filtered phenotype DataFrame with unrelated samples only
   :rtype: pd.DataFrame

----
.. _data_preparation:

Data Preparation
================

.. function:: stratified_train_test_split(genotype_df, phenotype_df, outcome_col, test_size=0.5, random_state=42, is_binary=True, geno_id_col=None, pheno_id_col=None)

   Split data into training and test sets with stratification.

   :param genotype_df: Genotype DataFrame (samples x variants)
   :type genotype_df: pd.DataFrame
   :param phenotype_df: Phenotype DataFrame
   :type phenotype_df: pd.DataFrame
   :param outcome_col: Name of outcome column for stratification
   :type outcome_col: str
   :param test_size: Proportion of samples in test set
   :type test_size: float
   :param random_state: Random seed for reproducibility
   :type random_state: int
   :param is_binary: Whether outcome is binary (enables stratification)
   :type is_binary: bool
   :param geno_id_col: Column name in genotype_df for sample IDs. If None, uses index
   :type geno_id_col: str, optional
   :param pheno_id_col: Column name in phenotype_df for sample IDs. If None, uses index
   :type pheno_id_col: str, optional
   :returns: Tuple of (train_geno, test_geno, train_pheno, test_pheno)
   :rtype: tuple
   :raises ValueError: If no common samples found or stratification fails

.. function:: impute_covariates(phenotype_df, covariate_cols, method='median', drop_na=False, verbose=True)

   Impute missing values in covariates.

   :param phenotype_df: Phenotype DataFrame with covariates
   :type phenotype_df: pd.DataFrame
   :param covariate_cols: List of covariate column names to impute
   :type covariate_cols: list
   :param method: Imputation method - 'drop', 'mean', 'median', 'mode', 'knn', 'missforest', 'mice'
   :type method: str
   :param drop_na: If True, drop rows with missing outcome after imputation
   :type drop_na: bool
   :param verbose: Print imputation information
   :type verbose: bool
   :returns: DataFrame with imputed covariates
   :rtype: pd.DataFrame

   .. note::
      For 'missforest' and 'mice', install: ``pip install missingpy``

----
.. _standard_gwas:

Standard GWAS
=============

.. function:: standard_gwas(genotype_df, phenotype_df, outcome, covariates, outcome_type='binary')

   Perform standard additive GWAS for comparison with EDGE.

   :param genotype_df: Genotype DataFrame
   :type genotype_df: pd.DataFrame
   :param phenotype_df: Phenotype DataFrame
   :type phenotype_df: pd.DataFrame
   :param outcome: Name of outcome column
   :type outcome: str
   :param covariates: List of covariate column names
   :type covariates: list
   :param outcome_type: 'binary' for logistic regression, 'continuous' for linear regression
   :type outcome_type: str
   :returns: DataFrame with variant_id, coef, pval, std_err
   :rtype: pd.DataFrame

.. function:: additive_gwas(genotype_df, phenotype_df, outcome, covariates, outcome_type='binary')

   Alias for :func:`standard_gwas`. Perform standard additive GWAS.

   :param genotype_df: Genotype DataFrame
   :type genotype_df: pd.DataFrame
   :param phenotype_df: Phenotype DataFrame
   :type phenotype_df: pd.DataFrame
   :param outcome: Name of outcome column
   :type outcome: str
   :param covariates: List of covariate column names
   :type covariates: list
   :param outcome_type: 'binary' for logistic regression, 'continuous' for linear regression
   :type outcome_type: str
   :returns: DataFrame with variant_id, coef, pval, std_err
   :rtype: pd.DataFrame

.. function:: cross_validated_edge_analysis(genotype_df, phenotype_df, outcome, covariates, outcome_type='binary', n_folds=5, n_jobs=8, random_state=42)

   Perform k-fold cross-validation for EDGE analysis.

   :param genotype_df: Genotype DataFrame
   :type genotype_df: pd.DataFrame
   :param phenotype_df: Phenotype DataFrame
   :type phenotype_df: pd.DataFrame
   :param outcome: Name of outcome column
   :type outcome: str
   :param covariates: List of covariate column names
   :type covariates: list
   :param outcome_type: 'binary' or 'continuous'
   :type outcome_type: str
   :param n_folds: Number of cross-validation folds
   :type n_folds: int
   :param n_jobs: Number of parallel jobs for EDGE analysis
   :type n_jobs: int
   :param random_state: Random seed for reproducibility
   :type random_state: int
   :returns: Tuple of (avg_alpha, meta_gwas_df, combined_alpha, combined_gwas)
   :rtype: tuple

----
.. _visualization:

Visualization
=============

.. function:: manhattan_plot(gwas_df, output='manhattan.png', title='EDGE GWAS Manhattan Plot', sig_threshold=5e-8, figsize=(14, 6), colors=None)

   Create Manhattan plot from EDGE GWAS results.

   :param gwas_df: DataFrame or list of DataFrames with columns 'chrom', 'pos', 'pval'
   :type gwas_df: pd.DataFrame or list
   :param output: Output filename for the plot
   :type output: str
   :param title: Plot title
   :type title: str
   :param sig_threshold: Genome-wide significance threshold
   :type sig_threshold: float
   :param figsize: Figure size as (width, height)
   :type figsize: tuple
   :param colors: List of two colors for alternating chromosomes
   :type colors: list, optional

.. function:: qq_plot(gwas_df, output='qq_plot.png', title='EDGE GWAS QQ Plot', figsize=(8, 8))

   Create QQ plot from EDGE GWAS results and calculate genomic inflation factor.

   :param gwas_df: DataFrame or list of DataFrames with column 'pval'
   :type gwas_df: pd.DataFrame or list
   :param output: Output filename for the plot
   :type output: str
   :param title: Plot title
   :type title: str
   :param figsize: Figure size as (width, height)
   :type figsize: tuple
   :returns: Genomic inflation factor (lambda_gc)
   :rtype: float

.. function:: plot_alpha_distribution(alpha_df, output='alpha_distribution.png', bins=50, figsize=(10, 6), xlim=None)

   Plot distribution of alpha values.

   :param alpha_df: DataFrame with 'alpha_value' column
   :type alpha_df: pd.DataFrame
   :param output: Output filename
   :type output: str
   :param bins: Number of histogram bins
   :type bins: int
   :param figsize: Figure size as (width, height)
   :type figsize: tuple
   :param xlim: Optional tuple (min, max) for x-axis limits. If None, uses full range
   :type xlim: tuple, optional

----
.. _input_output:

Input/Output
============

.. function:: save_results(gwas_df, alpha_df=None, output_prefix='edge_gwas', save_alpha=True)

   Save EDGE GWAS results to files.

   :param gwas_df: GWAS results DataFrame
   :type gwas_df: pd.DataFrame
   :param alpha_df: Alpha values DataFrame
   :type alpha_df: pd.DataFrame, optional
   :param output_prefix: Prefix for output files
   :type output_prefix: str
   :param save_alpha: Whether to save alpha values
   :type save_alpha: bool
   :returns: Dictionary with output file paths
   :rtype: dict

.. function:: load_alpha_values(alpha_file)

   Load pre-calculated alpha values.

   :param alpha_file: Path to alpha values file
   :type alpha_file: str
   :returns: DataFrame with alpha values
   :rtype: pd.DataFrame

.. function:: format_gwas_output_for_locuszoom(gwas_df, include_alpha=True, sort_by='pval', format_for_locuszoom=False)

   Format GWAS output for publication/reporting or LocusZoom upload.

   :param gwas_df: GWAS results DataFrame
   :type gwas_df: pd.DataFrame
   :param include_alpha: Include alpha-related columns
   :type include_alpha: bool
   :param sort_by: Column to sort by
   :type sort_by: str
   :param format_for_locuszoom: If True, format for LocusZoom upload with required columns
   :type format_for_locuszoom: bool
   :returns: Formatted DataFrame
   :rtype: pd.DataFrame

   .. note::
      For LocusZoom format, the output will be tab-delimited with columns:
      chrom, pos, ref, alt, pval, beta, se, eaf, and optionally alpha_value.
      The file should be sorted by chrom and pos, compressed with bgzip,
      and indexed with tabix for optimal LocusZoom performance.

.. function:: save_for_locuszoom(gwas_df, output_file, include_alpha=True, compress=True)

   Save GWAS results in LocusZoom-compatible format.

   :param gwas_df: GWAS results DataFrame
   :type gwas_df: pd.DataFrame
   :param output_file: Output file path (will add .gz if compress=True)
   :type output_file: str
   :param include_alpha: Include alpha_value column
   :type include_alpha: bool
   :param compress: If True, compress with gzip
   :type compress: bool

   .. note::
      For best performance with LocusZoom:
      
      1. Compress with bgzip: ``bgzip output_file.tsv``
      2. Index with tabix: ``tabix -s 1 -b 2 -e 2 output_file.tsv.gz``

.. function:: validate_locuszoom_format(gwas_df)

   Validate that GWAS results meet LocusZoom format requirements.

   :param gwas_df: GWAS results DataFrame
   :type gwas_df: pd.DataFrame
   :returns: Dictionary with validation results (valid, errors, warnings, info)
   :rtype: dict

.. function:: create_summary_report(gwas_df, alpha_df=None, significance_threshold=5e-8, output_file=None)

   Create a summary report of EDGE GWAS analysis.

   :param gwas_df: GWAS results DataFrame
   :type gwas_df: pd.DataFrame
   :param alpha_df: Alpha values DataFrame
   :type alpha_df: pd.DataFrame, optional
   :param significance_threshold: P-value threshold for significance
   :type significance_threshold: float
   :param output_file: Optional file to save report
   :type output_file: str, optional
   :returns: Summary report as string
   :rtype: str

----

.. _function_index:

Function Index
==============

**Core Analysis** (see :ref:`core_analysis`)
   * :class:`EDGEAnalysis`
   * :meth:`EDGEAnalysis.calculate_alpha`
   * :meth:`EDGEAnalysis.apply_alpha`
   * :meth:`EDGEAnalysis.run_full_analysis`
   * :meth:`EDGEAnalysis.get_skipped_snps`

**Data Loading** (see :ref:`data_loading`)
   * :func:`load_plink_data`
   * :func:`load_pgen_data`
   * :func:`load_vcf_data`
   * :func:`load_bgen_data`
   * :func:`prepare_phenotype_data`
   * :func:`download_test_files`

**Data Validation** (see :ref:`data_validation`)
   * :func:`validate_genotype_df`
   * :func:`validate_and_fix_encoding`
   * :func:`validate_phenotype_df`
   * :func:`validate_and_align_data`

**Quality Control** (see :ref:`quality_control`)
   * :func:`filter_genotype_data`
   * :func:`filter_variants_by_maf`
   * :func:`filter_variants_by_missing`
   * :func:`filter_samples_by_call_rate`
   * :func:`calculate_hwe_pvalues`
   * :func:`filter_variants_by_hwe`
   * :func:`check_case_control_balance`

**Population Structure Control** (see :ref:`population_structure`)
   * :func:`calculate_grm_gcta`
   * :func:`load_grm_gcta`
   * :func:`calculate_pca_sklearn`
   * :func:`calculate_pca_plink`
   * :func:`calculate_pca_pcair`
   * :func:`attach_pcs_to_phenotype`
   * :func:`get_pc_covariate_list`
   * :func:`identify_related_samples`
   * :func:`filter_related_samples`

**Data Preparation** (see :ref:`data_preparation`)
   * :func:`stratified_train_test_split`
   * :func:`impute_covariates`

**Standard GWAS** (see :ref:`standard_gwas`)
   * :func:`standard_gwas`
   * :func:`additive_gwas`
   * :func:`cross_validated_edge_analysis`

**Visualization** (see :ref:`visualization`)
   * :func:`manhattan_plot`
   * :func:`qq_plot`
   * :func:`plot_alpha_distribution`

**Input/Output** (see :ref:`input_output`)
   * :func:`save_results`
   * :func:`load_alpha_values`
   * :func:`format_gwas_output_for_locuszoom`
   * :func:`save_for_locuszoom`
   * :func:`validate_locuszoom_format`
   * :func:`create_summary_report`


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
