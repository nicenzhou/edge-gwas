.. _api_reference:

API Reference
=============

Complete API documentation for edge-gwas package.

Core Module
-----------

EDGEAnalysis Class
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from edge_gwas import EDGEAnalysis
   
   edge = EDGEAnalysis(
       outcome_type='binary',  # 'binary' or 'continuous'
       n_jobs=-1,              # Number of parallel jobs
       max_iter=1000,          # Max iterations for convergence
       verbose=True            # Print progress
   )

**Parameters:**

* ``outcome_type`` (str): Type of outcome variable

  * ``'binary'``: For case-control studies (logistic regression)
  * ``'continuous'``: For quantitative traits (linear regression)

* ``n_jobs`` (int): Number of CPU cores for parallel processing

  * ``-1``: Use all available cores
  * ``1``: Single-threaded (no parallelization)
  * ``n``: Use n cores

* ``max_iter`` (int): Maximum iterations for model convergence
* ``verbose`` (bool): Print progress messages

Methods
^^^^^^^

calculate_alpha()
"""""""""""""""""

Calculate encoding parameters from training data.

.. code-block:: python

   alpha_df = edge.calculate_alpha(
       genotype_data,
       phenotype_df,
       outcome,
       covariates,
       variant_info=None
   )

**Parameters:**

* ``genotype_data`` (pandas.DataFrame): Genotype matrix (samples × variants)
* ``phenotype_df`` (pandas.DataFrame): Phenotype data with outcome and covariates
* ``outcome`` (str): Name of outcome column
* ``covariates`` (list): List of covariate column names
* ``variant_info`` (pandas.DataFrame, optional): Variant information (chr, pos, etc.)

**Returns:**

* ``pandas.DataFrame``: Alpha values with columns:

  * ``variant_id``: SNP identifier
  * ``alpha_value``: Encoding parameter
  * ``eaf``: Effect allele frequency
  * ``coef_het``: Heterozygous coefficient
  * ``coef_hom``: Homozygous coefficient

apply_alpha()
"""""""""""""

Apply alpha values to test data for GWAS.

.. code-block:: python

   gwas_df = edge.apply_alpha(
       genotype_data,
       phenotype_df,
       outcome,
       covariates,
       alpha_values
   )

**Parameters:**

* ``genotype_data`` (pandas.DataFrame): Test genotype matrix
* ``phenotype_df`` (pandas.DataFrame): Test phenotype data
* ``outcome`` (str): Name of outcome column
* ``covariates`` (list): List of covariate column names
* ``alpha_values`` (pandas.DataFrame): Alpha values from training

**Returns:**

* ``pandas.DataFrame``: GWAS results with columns:

  * ``variant_id``: SNP identifier
  * ``chr``: Chromosome
  * ``pos``: Position
  * ``pval``: P-value
  * ``coef``: Effect coefficient
  * ``std_err``: Standard error
  * ``stat``: Test statistic
  * ``alpha_value``: Applied alpha

run_full_analysis()
"""""""""""""""""""

Complete two-stage EDGE analysis.

.. code-block:: python

   alpha_df, gwas_df = edge.run_full_analysis(
       train_genotype,
       train_phenotype,
       test_genotype,
       test_phenotype,
       outcome,
       covariates,
       variant_info=None,
       output_prefix=None
   )

**Parameters:**

* ``train_genotype`` (pandas.DataFrame): Training genotype data
* ``train_phenotype`` (pandas.DataFrame): Training phenotype data
* ``test_genotype`` (pandas.DataFrame): Test genotype data
* ``test_phenotype`` (pandas.DataFrame): Test phenotype data
* ``outcome`` (str): Outcome column name
* ``covariates`` (list): Covariate column names
* ``variant_info`` (pandas.DataFrame, optional): Variant information
* ``output_prefix`` (str, optional): Prefix for output files

**Returns:**

* ``tuple``: (alpha_df, gwas_df)

get_skipped_snps()
""""""""""""""""""

Get list of SNPs that failed convergence.

.. code-block:: python

   skipped = edge.get_skipped_snps()

**Returns:**

* ``list``: List of variant IDs that were skipped

Utilities Module
----------------

Data Loading Functions
~~~~~~~~~~~~~~~~~~~~~~

load_plink_data()
"""""""""""""""""

Load genotype data from PLINK binary files.

.. code-block:: python

   from edge_gwas.utils import load_plink_data
   
   genotype_df, variant_info = load_plink_data(
       bed_file='data.bed',
       bim_file='data.bim',
       fam_file='data.fam',
       verbose=True
   )

**Parameters:**

* ``bed_file`` (str): Path to .bed file
* ``bim_file`` (str): Path to .bim file
* ``fam_file`` (str): Path to .fam file
* ``verbose`` (bool): Print loading information

**Returns:**

* ``tuple``: (genotype_df, variant_info)

  * ``genotype_df``: Genotype matrix (samples × variants)
  * ``variant_info``: Variant information (chr, pos, ref, alt)

prepare_phenotype_data()
""""""""""""""""""""""""

Load and prepare phenotype data.

.. code-block:: python

   from edge_gwas.utils import prepare_phenotype_data
   
   phenotype_df = prepare_phenotype_data(
       phenotype_file='pheno.txt',
       outcome_col='disease',
       covariate_cols=['age', 'sex', 'PC1', 'PC2'],
       sample_id_col='IID',
       sep='\t',
       log_transform_outcome=False
   )

**Parameters:**

* ``phenotype_file`` (str): Path to phenotype file
* ``outcome_col`` (str): Name of outcome column
* ``covariate_cols`` (list): List of covariate columns
* ``sample_id_col`` (str): Sample ID column name
* ``sep`` (str): File delimiter
* ``log_transform_outcome`` (bool): Apply log transformation to outcome

**Returns:**

* ``pandas.DataFrame``: Formatted phenotype data

Data Processing Functions
~~~~~~~~~~~~~~~~~~~~~~~~~

stratified_train_test_split()
"""""""""""""""""""""""""""""

Split data into training and test sets.

.. code-block:: python

   from edge_gwas.utils import stratified_train_test_split
   
   train_g, test_g, train_p, test_p = stratified_train_test_split(
       genotype_df,
       phenotype_df,
       outcome_col='disease',
       test_size=0.5,
       random_state=42,
       is_binary=True
   )

**Parameters:**

* ``genotype_df`` (pandas.DataFrame): Genotype data
* ``phenotype_df`` (pandas.DataFrame): Phenotype data
* ``outcome_col`` (str): Outcome column for stratification
* ``test_size`` (float): Proportion for test set (0-1)
* ``random_state`` (int): Random seed for reproducibility
* ``is_binary`` (bool): Whether to stratify by outcome

**Returns:**

* ``tuple``: (train_geno, test_geno, train_pheno, test_pheno)

filter_variants_by_maf()
""""""""""""""""""""""""

Filter variants by minor allele frequency.

.. code-block:: python

   from edge_gwas.utils import filter_variants_by_maf
   
   filtered_geno = filter_variants_by_maf(
       genotype_df,
       min_maf=0.01,
       verbose=True
   )

**Parameters:**

* ``genotype_df`` (pandas.DataFrame): Genotype data
* ``min_maf`` (float): Minimum MAF threshold
* ``verbose`` (bool): Print filtering information

**Returns:**

* ``pandas.DataFrame``: Filtered genotype data

filter_variants_by_missing()
""""""""""""""""""""""""""""

Filter variants by missingness rate.

.. code-block:: python

   from edge_gwas.utils import filter_variants_by_missing
   
   filtered_geno = filter_variants_by_missing(
       genotype_df,
       max_missing=0.05,
       verbose=True
   )

**Parameters:**

* ``genotype_df`` (pandas.DataFrame): Genotype data
* ``max_missing`` (float): Maximum missingness rate (0-1)
* ``verbose`` (bool): Print filtering information

**Returns:**

* ``pandas.DataFrame``: Filtered genotype data

Statistical Functions
~~~~~~~~~~~~~~~~~~~~~

calculate_genomic_inflation()
"""""""""""""""""""""""""""""

Calculate genomic inflation factor (λ).

.. code-block:: python

   from edge_gwas.utils import calculate_genomic_inflation
   
   lambda_gc = calculate_genomic_inflation(pvalues)

**Parameters:**

* ``pvalues`` (pandas.Series or array): P-values from GWAS

**Returns:**

* ``float``: Genomic inflation factor (λ)

merge_alpha_with_gwas()
"""""""""""""""""""""""

Merge GWAS results with alpha values.

.. code-block:: python

   from edge_gwas.utils import merge_alpha_with_gwas
   
   merged_df = merge_alpha_with_gwas(gwas_df, alpha_df)

**Parameters:**

* ``gwas_df`` (pandas.DataFrame): GWAS results
* ``alpha_df`` (pandas.DataFrame): Alpha values

**Returns:**

* ``pandas.DataFrame``: Merged results with both GWAS and alpha information

Visualization Module
--------------------

manhattan_plot()
~~~~~~~~~~~~~~~~

Create Manhattan plot of GWAS results.

.. code-block:: python

   from edge_gwas.visualize import manhattan_plot
   
   manhattan_plot(
       gwas_df,
       output='manhattan.png',
       title='EDGE GWAS Manhattan Plot',
       sig_threshold=5e-8,
       suggestive_threshold=1e-5,
       figsize=(14, 6),
       colors=['#1f77b4', '#ff7f0e']
   )

**Parameters:**

* ``gwas_df`` (pandas.DataFrame): GWAS results with columns: chr, pos, pval
* ``output`` (str): Output filename
* ``title`` (str): Plot title
* ``sig_threshold`` (float): Genome-wide significance threshold
* ``suggestive_threshold`` (float): Suggestive significance threshold
* ``figsize`` (tuple): Figure size (width, height) in inches
* ``colors`` (list): Colors for alternating chromosomes

**Returns:**

* None (saves plot to file)

qq_plot()
~~~~~~~~~

Create QQ plot and calculate genomic inflation factor.

.. code-block:: python

   from edge_gwas.visualize import qq_plot
   
   lambda_gc = qq_plot(
       gwas_df,
       output='qq_plot.png',
       title='QQ Plot',
       figsize=(8, 8)
   )

**Parameters:**

* ``gwas_df`` (pandas.DataFrame): GWAS results with pval column
* ``output`` (str): Output filename
* ``title`` (str): Plot title
* ``figsize`` (tuple): Figure size (width, height) in inches

**Returns:**

* ``float``: Genomic inflation factor (λ)

plot_alpha_distribution()
~~~~~~~~~~~~~~~~~~~~~~~~~

Plot distribution of alpha values.

.. code-block:: python

   from edge_gwas.visualize import plot_alpha_distribution
   
   plot_alpha_distribution(
       alpha_df,
       output='alpha_distribution.png',
       bins=50,
       figsize=(10, 6),
       title='Alpha Value Distribution'
   )

**Parameters:**

* ``alpha_df`` (pandas.DataFrame): Alpha values with alpha_value column
* ``output`` (str): Output filename
* ``bins`` (int): Number of histogram bins
* ``figsize`` (tuple): Figure size (width, height)
* ``title`` (str): Plot title

**Returns:**

* None (saves plot to file)

I/O Handlers Module
-------------------

save_results()
~~~~~~~~~~~~~~

Save GWAS and alpha results to files.

.. code-block:: python

   from edge_gwas.io_handlers import save_results
   
   output_files = save_results(
       gwas_df,
       alpha_df=None,
       output_prefix='edge_gwas',
       save_alpha=True
   )

**Parameters:**

* ``gwas_df`` (pandas.DataFrame): GWAS results
* ``alpha_df`` (pandas.DataFrame, optional): Alpha values
* ``output_prefix`` (str): Prefix for output files
* ``save_alpha`` (bool): Whether to save alpha values

**Returns:**

* ``dict``: Dictionary with keys 'gwas' and 'alpha' containing file paths

load_alpha_values()
~~~~~~~~~~~~~~~~~~~

Load pre-calculated alpha values from file.

.. code-block:: python

   from edge_gwas.io_handlers import load_alpha_values
   
   alpha_df = load_alpha_values('alpha_values.txt')

**Parameters:**

* ``filename`` (str): Path to alpha values file

**Returns:**

* ``pandas.DataFrame``: Alpha values

create_summary_report()
~~~~~~~~~~~~~~~~~~~~~~~

Generate text summary of GWAS results.

.. code-block:: python

   from edge_gwas.io_handlers import create_summary_report
   
   report = create_summary_report(
       gwas_df,
       alpha_df=None,
       significance_threshold=5e-8,
       output_file='summary.txt'
   )

**Parameters:**

* ``gwas_df`` (pandas.DataFrame): GWAS results
* ``alpha_df`` (pandas.DataFrame, optional): Alpha values
* ``significance_threshold`` (float): P-value threshold for significance
* ``output_file`` (str, optional): Output file path

**Returns:**

* ``str``: Summary report text

Data Structures
---------------

Genotype DataFrame Format
~~~~~~~~~~~~~~~~~~~~~~~~~

**Expected format:**

* Index: Sample IDs
* Columns: Variant IDs
* Values: Genotype encoding (0, 1, 2, or NaN)

  * 0 = Reference homozygote
  * 1 = Heterozygote
  * 2 = Alternate homozygote
  * NaN = Missing

**Example:**

.. code-block:: python

   #              rs123  rs456  rs789
   # Sample1        0      1      2
   # Sample2        1      2      0
   # Sample3        2      0      1

Phenotype DataFrame Format
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Expected format:**

* Index: Sample IDs (matching genotype data)
* Columns: Outcome and covariate columns

**Example:**

.. code-block:: python

   #          disease  age  sex  PC1   PC2
   # Sample1        1   45    0  0.1  -0.2
   # Sample2        0   38    1 -0.3   0.4
   # Sample3        1   52    0  0.2   0.1

Variant Info DataFrame Format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Expected format:**

* Columns: variant_id, chr, pos, ref, alt, etc.

**Example:**

.. code-block:: python

   #    variant_id  chr      pos ref alt
   # 0       rs123    1  1234567   A   G
   # 1       rs456    1  2345678   C   T
   # 2       rs789    2  3456789   G   A

Alpha DataFrame Format
~~~~~~~~~~~~~~~~~~~~~~

**Output format from calculate_alpha():**

.. code-block:: python

   #    variant_id  alpha_value    eaf  coef_het  coef_hom
   # 0       rs123         0.45  0.235      0.23      0.51
   # 1       rs456         0.52  0.412      0.31      0.60
   # 2       rs789         0.38  0.189      0.19      0.50

GWAS DataFrame Format
~~~~~~~~~~~~~~~~~~~~~

**Output format from apply_alpha():**

.. code-block:: python

   #    variant_id  chr      pos       pval    coef  std_err    stat  alpha_value
   # 0       rs123    1  1234567  1.23e-08  0.0451   0.0089    5.07         0.45
   # 1       rs456    1  2345678  3.45e-06  0.0312   0.0067    4.66         0.52
   # 2       rs789    2  3456789  5.67e-10  0.0623   0.0095    6.56         0.38

Constants and Defaults
----------------------

.. code-block:: python

   # Significance thresholds
   GWAS_THRESHOLD = 5e-8          # Genome-wide significance
   SUGGESTIVE_THRESHOLD = 1e-5    # Suggestive threshold
   
   # Quality control defaults
   DEFAULT_MAF = 0.01             # Minimum MAF
   DEFAULT_MISSING = 0.05         # Maximum missingness
   DEFAULT_HWE = 1e-6             # HWE p-value threshold
   
   # Analysis defaults
   DEFAULT_TEST_SIZE = 0.5        # Train/test split ratio
   DEFAULT_N_JOBS = -1            # Use all CPU cores
   DEFAULT_MAX_ITER = 1000        # Maximum iterations

Exceptions
----------

.. code-block:: python

   class EdgeGWASError(Exception):
       """Base exception for edge-gwas."""
       pass
   
   class DataFormatError(EdgeGWASError):
       """Raised when input data format is invalid."""
       pass
   
   class ConvergenceError(EdgeGWASError):
       """Raised when model fails to converge."""
       pass
   
   class MissingDataError(EdgeGWASError):
       """Raised when required data is missing."""
       pass

Version Information
-------------------

.. code-block:: python

   import edge_gwas
   
   # Get version string
   print(edge_gwas.__version__)  # '0.1.0'
   
   # Get author information
   print(edge_gwas.__author__)   # 'Jiayan Zhou, Molly Ann Hall'
   
   # Get license
   print(edge_gwas.__license__)  # 'GPL-3.0'

See Also
--------

* :ref:`quickstart` - Quick start guide
* :ref:`examples` - Example workflows
* :ref:`statistical_model` - Statistical methodology
* :ref:`user_guide` - Detailed user guide
