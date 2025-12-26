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
       outcome_type='binary',           # 'binary' or 'continuous'
       outcome_transform=None,          # Transformation for continuous outcomes
       ols_method='bfgs',              # NEW in v0.1.1: Optimization method for OLS
       n_jobs=-1,                       # Number of parallel jobs
       max_iter=1000,                   # Max iterations for convergence
       verbose=True                     # Print progress
   )

**Parameters:**

* ``outcome_type`` (str): Type of outcome variable

  * ``'binary'``: For case-control studies (logistic regression)
  * ``'continuous'``: For quantitative traits (linear regression)

* ``outcome_transform`` (str, optional): **NEW in v0.1.1** - Transformation for continuous outcomes

  * ``None``: No transformation (default)
  * ``'log'``: Natural log transformation
  * ``'log10'``: Log base 10 transformation
  * ``'inverse_normal'``: Inverse normal transformation (parametric)
  * ``'rank_inverse_normal'``: Rank-based inverse normal transformation (robust to outliers)

* ``ols_method`` (str, optional): **NEW in v0.1.1** - Optimization method for OLS regression (continuous outcomes only)

  * ``'bfgs'``: Broyden-Fletcher-Goldfarb-Shannon algorithm (default)
    
    - Best for: General purpose, balanced speed and accuracy
    - Convergence: Fast, uses approximate Hessian
    - Memory: Moderate (O(p²) where p = number of parameters)
    - Recommended: Most analyses

  * ``'newton'``: Newton-Raphson algorithm
    
    - Best for: High precision requirements
    - Convergence: Very fast when close to optimum (quadratic)
    - Memory: High (requires exact Hessian computation)
    - Recommended: Small datasets, need high accuracy

  * ``'lbfgs'``: Limited-memory BFGS
    
    - Best for: Large datasets (N > 100,000)
    - Convergence: Slightly slower than BFGS
    - Memory: Low (O(mp) where m ≈ 10)
    - Recommended: Biobank-scale data, memory-constrained systems

  * ``'nm'``: Nelder-Mead simplex algorithm
    
    - Best for: Robustness, convergence issues with gradient methods
    - Convergence: Slow, derivative-free
    - Memory: Low
    - Recommended: Troubleshooting convergence failures

  * ``'cg'``: Conjugate Gradient
    
    - Best for: Sparse problems, many covariates
    - Convergence: Moderate
    - Memory: Low
    - Recommended: Large covariate sets

  * ``'ncg'``: Newton Conjugate Gradient
    
    - Best for: Large-scale problems with >100 covariates
    - Convergence: Good, uses Hessian-vector products
    - Memory: Moderate
    - Recommended: High-dimensional analyses

  * ``'powell'``: Powell's conjugate direction method
    
    - Best for: Derivative-free optimization
    - Convergence: Slow, robust
    - Memory: Low
    - Recommended: Non-smooth objective functions (rare in GWAS)

  * ``'basinhopping'``: Basin-hopping global optimization
    
    - Best for: Multiple local minima (very rare in GWAS)
    - Convergence: Very slow, explores parameter space thoroughly
    - Memory: Moderate
    - Recommended: Only when suspecting multiple solutions

  **Note:** Optimization method only affects continuous outcomes with linear regression. 
  Binary outcomes always use BFGS for logistic regression.

* ``n_jobs`` (int): Number of CPU cores for parallel processing

  * ``-1``: Use all available cores
  * ``1``: Single-threaded (no parallelization)
  * ``n``: Use n cores

* ``max_iter`` (int): Maximum iterations for model convergence

  * Default: 1000 (usually sufficient for BFGS, Newton)
  * Increase to 2000-5000 for Nelder-Mead or if convergence warnings occur
  * Reduce to 500 if most variants converge quickly

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
       variant_info=None,
       grm_matrix=None,        # NEW in v0.1.1
       grm_sample_ids=None     # NEW in v0.1.1
   )

**Parameters:**

* ``genotype_data`` (pandas.DataFrame): Genotype matrix (samples × variants)
* ``phenotype_df`` (pandas.DataFrame): Phenotype data with outcome and covariates
* ``outcome`` (str): Name of outcome column
* ``covariates`` (list): List of covariate column names
* ``variant_info`` (pandas.DataFrame, optional): Variant information (chr, pos, etc.)
* ``grm_matrix`` (numpy.ndarray, optional): **NEW in v0.1.1** - GRM matrix for population structure control
* ``grm_sample_ids`` (pandas.DataFrame, optional): **NEW in v0.1.1** - Sample IDs corresponding to GRM

**Returns:**

* ``pandas.DataFrame``: Alpha values with columns:

  * ``variant_id``: SNP identifier
  * ``alpha_value``: Encoding parameter
  * ``eaf``: Effect allele frequency
  * ``coef_het``: Heterozygous coefficient
  * ``coef_hom``: Homozygous coefficient

**Example with GRM:**

.. code-block:: python

   from edge_gwas.utils import load_grm_gcta
   
   # Load GRM for population structure control
   grm_matrix, grm_ids = load_grm_gcta('grm_prefix')
   
   # Calculate alpha with GRM
   alpha_df = edge.calculate_alpha(
       genotype_data=train_geno,
       phenotype_df=train_pheno,
       outcome='trait',
       covariates=['age', 'sex'],
       grm_matrix=grm_matrix,
       grm_sample_ids=grm_ids
   )

apply_alpha()
"""""""""""""

Apply alpha values to test data for GWAS.

.. code-block:: python

   gwas_df = edge.apply_alpha(
       genotype_data,
       phenotype_df,
       outcome,
       covariates,
       alpha_values,
       grm_matrix=None,        # NEW in v0.1.1
       grm_sample_ids=None     # NEW in v0.1.1
   )

**Parameters:**

* ``genotype_data`` (pandas.DataFrame): Test genotype matrix
* ``phenotype_df`` (pandas.DataFrame): Test phenotype data
* ``outcome`` (str): Name of outcome column
* ``covariates`` (list): List of covariate column names
* ``alpha_values`` (pandas.DataFrame): Alpha values from training
* ``grm_matrix`` (numpy.ndarray, optional): **NEW in v0.1.1** - GRM matrix for mixed model
* ``grm_sample_ids`` (pandas.DataFrame, optional): **NEW in v0.1.1** - Sample IDs for GRM

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
       grm_matrix=None,        # NEW in v0.1.1
       grm_sample_ids=None,    # NEW in v0.1.1
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
* ``grm_matrix`` (numpy.ndarray, optional): **NEW in v0.1.1** - GRM for both training and testing
* ``grm_sample_ids`` (pandas.DataFrame, optional): **NEW in v0.1.1** - Sample IDs for GRM
* ``output_prefix`` (str, optional): Prefix for output files

**Returns:**

* ``tuple``: (alpha_df, gwas_df)

**Example with outcome transformation and GRM:**

.. code-block:: python

   from edge_gwas import EDGEAnalysis
   from edge_gwas.utils import load_grm_gcta, calculate_pca_plink, attach_pcs_to_phenotype
   
   # Calculate PCA
   pca_df = calculate_pca_plink('genotypes', n_pcs=10)
   train_pheno = attach_pcs_to_phenotype(train_pheno, pca_df, n_pcs=10)
   test_pheno = attach_pcs_to_phenotype(test_pheno, pca_df, n_pcs=10)
   
   # Load GRM
   grm_matrix, grm_ids = load_grm_gcta('grm_prefix')
   
   # Initialize EDGE with outcome transformation
   edge = EDGEAnalysis(
       outcome_type='continuous',
       outcome_transform='rank_inverse_normal'
   )
   
   # Run full analysis
   alpha_df, gwas_df = edge.run_full_analysis(
       train_geno, train_pheno,
       test_geno, test_pheno,
       outcome='trait',
       covariates=['age', 'sex'] + [f'PC{i}' for i in range(1, 11)],
       grm_matrix=grm_matrix,
       grm_sample_ids=grm_ids,
       output_prefix='results/edge'
   )

get_skipped_snps()
""""""""""""""""""

Get list of SNPs that failed convergence.

.. code-block:: python

   skipped = edge.get_skipped_snps()

**Returns:**

* ``list``: List of variant IDs that were skipped

Utilities Module
----------------

Population Structure Control (NEW in v0.1.1)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

calculate_grm_gcta()
""""""""""""""""""""

Calculate genetic relationship matrix (GRM) using GCTA.

.. code-block:: python

   from edge_gwas.utils import calculate_grm_gcta
   
   grm_prefix = calculate_grm_gcta(
       plink_prefix='genotypes',
       output_prefix='output/grm',
       maf_threshold=0.01,
       method='grm',
       max_threads=8,
       verbose=True
   )

**Parameters:**

* ``plink_prefix`` (str): Prefix for PLINK binary files (.bed/.bim/.fam)
* ``output_prefix`` (str, optional): Prefix for output GRM files (default: temp directory)
* ``maf_threshold`` (float): MAF threshold for variant filtering (default: 0.01)
* ``method`` (str): GRM calculation method

  * ``'grm'``: Full GRM (default)
  * ``'grm-sparse'``: Sparse GRM (faster for large cohorts)

* ``max_threads`` (int): Maximum number of threads to use (default: 1)
* ``verbose`` (bool): Print progress information

**Returns:**

* ``str``: Path to output GRM prefix

**Output files:**

* ``prefix.grm.bin``: GRM values (lower triangle, binary format)
* ``prefix.grm.N.bin``: Number of SNPs used for each pair
* ``prefix.grm.id``: Sample IDs (FID and IID)

**Note:**

Requires GCTA to be installed. Install with:

.. code-block:: bash

   edge-gwas-install-tools

load_grm_gcta()
"""""""""""""""

Load GRM calculated by GCTA.

.. code-block:: python

   from edge_gwas.utils import load_grm_gcta
   
   grm_matrix, sample_ids = load_grm_gcta(
       grm_prefix='output/grm',
       verbose=True
   )

**Parameters:**

* ``grm_prefix`` (str): Prefix for GRM files (without .grm.bin extension)
* ``verbose`` (bool): Print loading information

**Returns:**

* ``tuple``: (grm_matrix, sample_ids_df)

  * ``grm_matrix`` (numpy.ndarray): n_samples × n_samples symmetric GRM matrix
  * ``sample_ids_df`` (pandas.DataFrame): DataFrame with FID and IID columns

**Example:**

.. code-block:: python

   # Calculate and load GRM
   grm_prefix = calculate_grm_gcta('genotypes', output_prefix='grm/output')
   grm_matrix, grm_ids = load_grm_gcta('grm/output')
   
   print(f"GRM shape: {grm_matrix.shape}")
   print(f"Mean diagonal: {np.diag(grm_matrix).mean():.3f}")

identify_related_samples()
"""""""""""""""""""""""""""

Identify pairs of related samples based on GRM threshold.

.. code-block:: python

   from edge_gwas.utils import identify_related_samples
   
   related_pairs = identify_related_samples(
       grm_matrix,
       sample_ids,
       threshold=0.0884,  # ~3rd degree relatives
       verbose=True
   )

**Parameters:**

* ``grm_matrix`` (numpy.ndarray): GRM matrix
* ``sample_ids`` (pandas.DataFrame): Sample IDs from load_grm_gcta()
* ``threshold`` (float): Relatedness threshold

  * ``0.354``: 1st degree (parent-offspring, full siblings)
  * ``0.177``: 2nd degree (half-siblings, grandparent-grandchild)
  * ``0.0884``: 3rd degree (first cousins) - default

* ``verbose`` (bool): Print summary statistics

**Returns:**

* ``pandas.DataFrame``: Related pairs with columns:

  * ``IID1``: First sample ID
  * ``IID2``: Second sample ID
  * ``kinship``: Kinship coefficient

**Example:**

.. code-block:: python

   grm_matrix, sample_ids = load_grm_gcta('grm_prefix')
   
   # Find 2nd degree or closer relatives
   related = identify_related_samples(grm_matrix, sample_ids, threshold=0.177)
   print(f"Found {len(related)} related pairs")

filter_related_samples()
"""""""""""""""""""""""""

Filter out related samples to create an unrelated subset.

.. code-block:: python

   from edge_gwas.utils import filter_related_samples
   
   unrelated_pheno = filter_related_samples(
       phenotype_df,
       grm_matrix,
       sample_ids,
       threshold=0.0884,
       method='greedy',
       sample_id_col=None,
       verbose=True
   )

**Parameters:**

* ``phenotype_df`` (pandas.DataFrame): Phenotype DataFrame
* ``grm_matrix`` (numpy.ndarray): GRM matrix
* ``sample_ids`` (pandas.DataFrame): Sample IDs from load_grm_gcta()
* ``threshold`` (float): Relatedness threshold
* ``method`` (str): Method for selecting unrelated samples

  * ``'greedy'``: Iteratively remove sample with most relatives (default)
  * ``'random'``: Randomly remove one from each related pair

* ``sample_id_col`` (str, optional): Column name in phenotype_df for sample IDs (if None, uses index)
* ``verbose`` (bool): Print filtering information

**Returns:**

* ``pandas.DataFrame``: Filtered phenotype DataFrame with unrelated samples only

**Example:**

.. code-block:: python

   # Remove related samples for QC
   unrelated_pheno = filter_related_samples(
       phenotype_df, grm_matrix, sample_ids, 
       threshold=0.0884, method='greedy'
   )
   print(f"Unrelated samples: {len(unrelated_pheno)}")

Principal Component Analysis (NEW in v0.1.1)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

calculate_pca_plink()
"""""""""""""""""""""

Calculate principal components using PLINK2.

.. code-block:: python

   from edge_gwas.utils import calculate_pca_plink
   
   pca_df = calculate_pca_plink(
       plink_prefix='genotypes',
       n_pcs=10,
       output_prefix=None,
       maf_threshold=0.01,
       ld_window=50,
       ld_step=5,
       ld_r2=0.2,
       approx=False,
       approx_samples=5000,
       verbose=True
   )

**Parameters:**

* ``plink_prefix`` (str): Prefix for PLINK binary files
* ``n_pcs`` (int): Number of principal components to calculate (default: 10)
* ``output_prefix`` (str, optional): Prefix for output files (default: temp directory)
* ``maf_threshold`` (float): MAF threshold for variant filtering (default: 0.01)
* ``ld_window`` (int): Window size for LD pruning in kb (default: 50)
* ``ld_step`` (int): Step size for LD pruning (default: 5)
* ``ld_r2`` (float): R² threshold for LD pruning (default: 0.2)
* ``approx`` (bool): Use approximate PCA for large cohorts (default: False)
* ``approx_samples`` (int): Number of samples for approximate PCA (default: 5000)
* ``verbose`` (bool): Print progress information

**Returns:**

* ``pandas.DataFrame``: PCA results with IID as index and PC1, PC2, ..., PCn as columns

**Note:**

Requires PLINK2 to be installed. Install with:

.. code-block:: bash

   edge-gwas-install-tools

**Example:**

.. code-block:: python

   # Standard PCA for <5000 samples
   pca_df = calculate_pca_plink('genotypes', n_pcs=10)
   
   # Approximate PCA for large cohorts (>5000 samples)
   pca_df = calculate_pca_plink(
       'genotypes', 
       n_pcs=20, 
       approx=True, 
       approx_samples=5000
   )
   
   print(f"Calculated {pca_df.shape[1]} PCs for {len(pca_df)} samples")

calculate_pca_pcair()
"""""""""""""""""""""

Calculate PC-AiR (Principal Components - Analysis in Related samples).

.. code-block:: python

   from edge_gwas.utils import calculate_pca_pcair
   
   pca_df = calculate_pca_pcair(
       plink_prefix='genotypes',
       n_pcs=10,
       kinship_matrix=None,
       divergence_matrix=None,
       output_prefix=None,
       kin_threshold=0.0884,
       div_threshold=-0.0884,
       maf_threshold=0.01,
       verbose=True
   )

**Parameters:**

* ``plink_prefix`` (str): Prefix for PLINK binary files
* ``n_pcs`` (int): Number of principal components to calculate (default: 10)
* ``kinship_matrix`` (str, optional): Path to GCTA GRM prefix (if None, will compute automatically)
* ``divergence_matrix`` (str, optional): Path to divergence matrix (optional)
* ``output_prefix`` (str, optional): Prefix for output files
* ``kin_threshold`` (float): Kinship threshold for defining relatedness (default: 0.0884 ~ 2nd degree)
* ``div_threshold`` (float): Divergence threshold (default: -0.0884)
* ``maf_threshold`` (float): MAF threshold for GRM calculation if kinship_matrix is None
* ``verbose`` (bool): Print progress information

**Returns:**

* ``pandas.DataFrame``: PC-AiR results with IID as index and PC1, PC2, ..., PCn as columns

**Note:**

Requires R with GENESIS, SNPRelate, and gdsfmt packages. Install with:

.. code-block:: bash

   edge-gwas-install-tools

**Reference:**

Conomos et al. (2015) Genetic Epidemiology 39(4):276-293

**Example:**

.. code-block:: python

   # PC-AiR for related samples
   pca_df = calculate_pca_pcair(
       plink_prefix='genotypes',
       n_pcs=10,
       kin_threshold=0.0884
   )
   
   # PC-AiR with pre-computed GRM
   pca_df = calculate_pca_pcair(
       plink_prefix='genotypes',
       n_pcs=10,
       kinship_matrix='grm_prefix'
   )

calculate_pca_sklearn()
"""""""""""""""""""""""

Calculate principal components using scikit-learn (basic PCA without relatedness correction).

.. code-block:: python

   from edge_gwas.utils import calculate_pca_sklearn
   
   pca_df = calculate_pca_sklearn(
       genotype_df,
       n_pcs=10,
       verbose=True
   )

**Parameters:**

* ``genotype_df`` (pandas.DataFrame): Genotype matrix (samples × variants)
* ``n_pcs`` (int): Number of principal components to calculate (default: 10)
* ``verbose`` (bool): Print progress information

**Returns:**

* ``pandas.DataFrame``: PCA results with IID as index and PC1, PC2, ..., PCn as columns

**Note:**

This is a basic PCA without correction for relatedness. For robust PCA, use ``calculate_pca_plink()`` or ``calculate_pca_pcair()``.

**Example:**

.. code-block:: python

   # Quick PCA for exploratory analysis
   pca_df = calculate_pca_sklearn(genotype_df, n_pcs=10)
   
   # Check variance explained
   print("PCA complete")

attach_pcs_to_phenotype()
""""""""""""""""""""""""""

Attach principal components to phenotype DataFrame.

.. code-block:: python

   from edge_gwas.utils import attach_pcs_to_phenotype
   
   pheno_with_pcs = attach_pcs_to_phenotype(
       phenotype_df,
       pca_df,
       n_pcs=10,
       pc_prefix='PC',
       sample_id_col=None,
       verbose=True
   )

**Parameters:**

* ``phenotype_df`` (pandas.DataFrame): Phenotype DataFrame
* ``pca_df`` (pandas.DataFrame): PCA DataFrame with IID as index
* ``n_pcs`` (int): Number of PCs to attach (default: 10)
* ``pc_prefix`` (str): Prefix for PC column names (default: 'PC')
* ``sample_id_col`` (str, optional): Column name in phenotype_df for sample IDs (if None, uses index)
* ``verbose`` (bool): Print information about merging

**Returns:**

* ``pandas.DataFrame``: Phenotype DataFrame with PC columns added

**Example:**

.. code-block:: python

   # Calculate PCA
   pca_df = calculate_pca_plink('genotypes', n_pcs=10)
   
   # Attach to phenotype (when phenotype has IID as index)
   pheno_with_pcs = attach_pcs_to_phenotype(pheno_df, pca_df, n_pcs=10)
   
   # Attach to phenotype (when phenotype has IID as column)
   pheno_with_pcs = attach_pcs_to_phenotype(
       pheno_df, pca_df, n_pcs=10, sample_id_col='IID'
   )

get_pc_covariate_list()
"""""""""""""""""""""""

Generate list of PC covariate names for use in EDGE analysis.

.. code-block:: python

   from edge_gwas.utils import get_pc_covariate_list
   
   pc_list = get_pc_covariate_list(n_pcs=10, pc_prefix='PC')

**Parameters:**

* ``n_pcs`` (int): Number of PCs
* ``pc_prefix`` (str): Prefix for PC column names (default: 'PC')

**Returns:**

* ``list``: List of PC column names ['PC1', 'PC2', ..., 'PCn']

**Example:**

.. code-block:: python

   # Generate PC covariate list
   pc_list = get_pc_covariate_list(10)
   # Returns: ['PC1', 'PC2', 'PC3', ..., 'PC10']
   
   # Use in EDGE analysis
   covariates = ['age', 'sex'] + get_pc_covariate_list(10)
   alpha_df, gwas_df = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       outcome='disease',
       covariates=covariates
   )

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

load_pgen_data()
""""""""""""""""

**NEW in v0.1.1** - Load PLINK 2 binary format data (.pgen/.pvar/.psam).

.. code-block:: python

   from edge_gwas.utils import load_pgen_data
   
   genotype_df, variant_info = load_pgen_data(
       pgen_file='data.pgen',
       pvar_file='data.pvar',
       psam_file='data.psam',
       verbose=True
   )

**Parameters:**

* ``pgen_file`` (str): Path to .pgen file
* ``pvar_file`` (str): Path to .pvar file  
* ``psam_file`` (str): Path to .psam file
* ``verbose`` (bool): Print loading information

**Returns:**

* ``tuple``: (genotype_df, variant_info)

**Note:**

Requires ``pgenlib`` package: 

.. code-block:: bash

   pip install pgenlib

load_bgen_data()
""""""""""""""""

**NEW in v0.1.1** - Load BGEN format data with dosages.

.. code-block:: python

   from edge_gwas.utils import load_bgen_data
   
   genotype_df, variant_info = load_bgen_data(
       bgen_file='data.bgen',
       sample_file='data.sample',  # Optional
       verbose=True
   )

**Parameters:**

* ``bgen_file`` (str): Path to .bgen file
* ``sample_file`` (str, optional): Path to .sample file (optional if embedded in BGEN)
* ``verbose`` (bool): Print loading information

**Returns:**

* ``tuple``: (genotype_df, variant_info)
  
  * Genotypes are dosages (0-2 continuous values)

**Note:**

Requires ``bgen-reader`` package: 

.. code-block:: bash

   pip install bgen-reader

load_vcf_data()
"""""""""""""""

**NEW in v0.1.1** - Load VCF/VCF.GZ format data.

.. code-block:: python

   from edge_gwas.utils import load_vcf_data
   
   genotype_df, variant_info = load_vcf_data(
       vcf_file='data.vcf.gz',
       dosage=True,    # Use DS field if available
       verbose=True
   )

**Parameters:**

* ``vcf_file`` (str): Path to .vcf or .vcf.gz file
* ``dosage`` (bool): If True, use dosages (DS field); if False, use hard calls (GT field)
* ``verbose`` (bool): Print loading information

**Returns:**

* ``tuple``: (genotype_df, variant_info)

**Note:**

Requires ``cyvcf2`` package: 

.. code-block:: bash

   pip install cyvcf2

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
* ``log_transform_outcome`` (bool): Apply log transformation to outcome (deprecated, use ``outcome_transform`` in ``EDGEAnalysis`` instead)

**Returns:**

* ``pandas.DataFrame``: Formatted phenotype data

Data Processing Functions
~~~~~~~~~~~~~~~~~~~~~~~~~

stratified_train_test_split()
"""""""""""""""""""""""""""""

Split data into training and test sets with stratification and flexible sample ID matching.

.. code-block:: python

   from edge_gwas.utils import stratified_train_test_split
   
   train_g, test_g, train_p, test_p = stratified_train_test_split(
       genotype_df,
       phenotype_df,
       outcome_col='disease',
       test_size=0.5,
       random_state=42,
       is_binary=True,
       geno_id_col=None,
       pheno_id_col=None
   )

**Parameters:**

* ``genotype_df`` (pandas.DataFrame): Genotype data
* ``phenotype_df`` (pandas.DataFrame): Phenotype data
* ``outcome_col`` (str): Outcome column for stratification
* ``test_size`` (float): Proportion for test set (0-1)
* ``random_state`` (int): Random seed for reproducibility
* ``is_binary`` (bool): Whether to stratify by outcome
* ``geno_id_col`` (str, optional): **NEW in v0.1.1** - Column or index name for sample IDs in genotype_df. If None, uses the DataFrame index (default: None)
* ``pheno_id_col`` (str, optional): **NEW in v0.1.1** - Column or index name for sample IDs in phenotype_df. If None, uses the DataFrame index (default: None)

**Returns:**

* ``tuple``: (train_geno, test_geno, train_pheno, test_pheno)

**Examples:**

.. code-block:: python

   # Example 1: Basic usage with matching indices (recommended)
   pheno = pheno.set_index('IID')  # Set IID as index
   train_g, test_g, train_p, test_p = stratified_train_test_split(
       geno, pheno, 'disease', test_size=0.5
   )
   
   # Example 2: Specify different ID columns for genotype and phenotype
   train_g, test_g, train_p, test_p = stratified_train_test_split(
       geno, pheno, 'disease',
       geno_id_col='sample_id',  # Use 'sample_id' column/index from genotype
       pheno_id_col='IID',        # Use 'IID' column from phenotype
       test_size=0.5
   )
   
   # Example 3: Genotype uses index, phenotype uses column
   train_g, test_g, train_p, test_p = stratified_train_test_split(
       geno, pheno, 'disease',
       geno_id_col=None,          # Use genotype index
       pheno_id_col='IID',        # Use phenotype 'IID' column
       test_size=0.3,
       random_state=123
   )

**Note:**

The function automatically handles type mismatches between sample IDs (e.g., strings vs integers) by converting both to strings for comparison, then mapping back to original types for proper indexing. This ensures samples are matched correctly even if genotype IDs are strings ('11') and phenotype IDs are integers (11).

**Tip:**

For simplest usage, ensure both DataFrames have matching sample IDs as their index:

.. code-block:: python

   # Prepare data with matching indices
   pheno = pheno.set_index('IID')
   
   # Then use without specifying ID columns
   train_g, test_g, train_p, test_p = stratified_train_test_split(
       geno, pheno, 'disease'
   )

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

* ``genotype_df`` (pandas.DataFrame): Genotype data (works with both hard calls and dosages)
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

filter_samples_by_call_rate()
"""""""""""""""""""""""""""""

**NEW in v0.1.1** - Filter samples by genotype call rate.

.. code-block:: python

   from edge_gwas.utils import filter_samples_by_call_rate
   
   filtered_geno, filtered_pheno = filter_samples_by_call_rate(
       genotype_df,
       phenotype_df,
       min_call_rate=0.95,
       verbose=True
   )

**Parameters:**

* ``genotype_df`` (pandas.DataFrame): Genotype data
* ``phenotype_df`` (pandas.DataFrame): Phenotype data
* ``min_call_rate`` (float): Minimum call rate (proportion of non-missing genotypes)
* ``verbose`` (bool): Print filtering information

**Returns:**

* ``tuple``: (filtered_genotype_df, filtered_phenotype_df)

**Example:**

.. code-block:: python

   # Remove samples with >5% missing genotypes
   geno, pheno = filter_samples_by_call_rate(geno, pheno, min_call_rate=0.95)

filter_variants_by_hwe()
""""""""""""""""""""""""

**NEW in v0.1.1** - Filter variants by Hardy-Weinberg Equilibrium p-value.

.. code-block:: python

   from edge_gwas.utils import filter_variants_by_hwe
   
   filtered_geno = filter_variants_by_hwe(
       genotype_df,
       hwe_threshold=1e-6,
       verbose=True
   )

**Parameters:**

* ``genotype_df`` (pandas.DataFrame): Genotype data
* ``hwe_threshold`` (float): Minimum HWE p-value threshold
* ``verbose`` (bool): Print filtering information

**Returns:**

* ``pandas.DataFrame``: Filtered genotype data

**Example:**

.. code-block:: python

   # Filter variants deviating from HWE
   geno = filter_variants_by_hwe(geno, hwe_threshold=1e-6)

check_case_control_balance()
""""""""""""""""""""""""""""

**NEW in v0.1.1** - Check case/control balance in binary outcome.

.. code-block:: python

   from edge_gwas.utils import check_case_control_balance
   
   balance = check_case_control_balance(
       phenotype_df,
       outcome_col='disease',
       verbose=True
   )

**Parameters:**

* ``phenotype_df`` (pandas.DataFrame): Phenotype data
* ``outcome_col`` (str): Name of outcome column
* ``verbose`` (bool): Print balance information

**Returns:**

* ``dict``: Dictionary with keys:
  
  * ``case_count`` (int): Number of cases
  * ``control_count`` (int): Number of controls  
  * ``ratio`` (float): Case/control ratio

**Example:**

.. code-block:: python

   balance = check_case_control_balance(pheno, 'disease')
   print(f"Cases: {balance['case_count']}, Controls: {balance['control_count']}")
   print(f"Ratio: {balance['ratio']:.2f}")

calculate_hwe_pvalues()
"""""""""""""""""""""""

**NEW in v0.1.1** - Calculate Hardy-Weinberg Equilibrium p-values for each variant.

.. code-block:: python

   from edge_gwas.utils import calculate_hwe_pvalues
   
   hwe_pvals = calculate_hwe_pvalues(
       genotype_df,
       verbose=True
   )

**Parameters:**

* ``genotype_df`` (pandas.DataFrame): Genotype data
* ``verbose`` (bool): Print calculation information

**Returns:**

* ``pandas.Series``: HWE p-values for each variant

**Example:**

.. code-block:: python

   hwe_pvals = calculate_hwe_pvalues(genotype_df)
   
   # Find variants deviating from HWE
   hwe_violations = hwe_pvals[hwe_pvals < 1e-6]
   print(f"Variants violating HWE: {len(hwe_violations)}")

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

additive_gwas()
"""""""""""""""

**NEW in v0.1.1** - Perform standard additive GWAS for comparison with EDGE.

.. code-block:: python

   from edge_gwas.utils import additive_gwas
   
   results_df = additive_gwas(
       genotype_df,
       phenotype_df,
       outcome='disease',
       covariates=['age', 'sex', 'PC1', 'PC2'],
       outcome_type='binary'
   )

**Parameters:**

* ``genotype_df`` (pandas.DataFrame): Genotype data
* ``phenotype_df`` (pandas.DataFrame): Phenotype data
* ``outcome`` (str): Name of outcome column
* ``covariates`` (list): List of covariate column names
* ``outcome_type`` (str): ``'binary'`` for logistic regression, ``'continuous'`` for linear regression

**Returns:**

* ``pandas.DataFrame``: Results with columns:
  
  * ``variant_id``: SNP identifier
  * ``coef``: Effect coefficient
  * ``std_err``: Standard error
  * ``pval``: P-value

**Example:**

.. code-block:: python

   # Binary outcome (logistic regression)
   additive_results = additive_gwas(
       geno, pheno, 'disease', 
       ['age', 'sex'], 
       outcome_type='binary'
   )
   
   # Continuous outcome (linear regression)
   additive_results = additive_gwas(
       geno, pheno, 'BMI', 
       ['age', 'sex'], 
       outcome_type='continuous'
   )

cross_validated_edge_analysis()
""""""""""""""""""""""""""""""""

**NEW in v0.1.1** - Perform k-fold cross-validation for EDGE analysis.

.. code-block:: python

   from edge_gwas.utils import cross_validated_edge_analysis
   
   avg_alpha, meta_gwas, all_alpha, all_gwas = cross_validated_edge_analysis(
       genotype_df,
       phenotype_df,
       outcome='disease',
       covariates=['age', 'sex', 'PC1', 'PC2'],
       outcome_type='binary',
       n_folds=5,
       n_jobs=8,
       random_state=42
   )

**Parameters:**

* ``genotype_df`` (pandas.DataFrame): Genotype data
* ``phenotype_df`` (pandas.DataFrame): Phenotype data
* ``outcome`` (str): Name of outcome column
* ``covariates`` (list): List of covariate column names
* ``outcome_type`` (str): ``'binary'`` or ``'continuous'``
* ``n_folds`` (int): Number of cross-validation folds
* ``n_jobs`` (int): Number of parallel jobs for EDGE analysis
* ``random_state`` (int): Random seed for reproducibility

**Returns:**

* ``tuple``: (avg_alpha, meta_gwas_df, combined_alpha, combined_gwas)
  
  * ``avg_alpha``: Averaged alpha values across folds
  * ``meta_gwas_df``: Meta-analyzed GWAS results (Fisher's method)
  * ``combined_alpha``: All alpha values from all folds
  * ``combined_gwas``: All GWAS results from all folds

**Example:**

.. code-block:: python

   # Run 5-fold cross-validation
   avg_alpha, meta_gwas, all_alpha, all_gwas = cross_validated_edge_analysis(
       genotype_df=geno,
       phenotype_df=pheno,
       outcome='disease',
       covariates=['age', 'sex'] + [f'PC{i}' for i in range(1, 11)],
       outcome_type='binary',
       n_folds=5
   )
   
   # Check alpha stability
   print(f"Mean alpha std: {avg_alpha['alpha_std'].mean():.3f}")
   
   # Find unstable variants
   unstable = avg_alpha[avg_alpha['alpha_std'] > 0.3]
   print(f"Unstable variants: {len(unstable)}")

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

download_test_files()
"""""""""""""""""""""

**NEW in v0.1.1** - Download test data files from GitHub repository for tutorials and testing.

.. code-block:: python

   from edge_gwas import download_test_files
   
   results = download_test_files(
       output_dir='tests',
       version='v0.1.1',
       overwrite=False,
       verbose=True
   )

**Parameters:**

* ``output_dir`` (str): Directory to save test files (default: 'tests')
* ``version`` (str): GitHub release version tag to download from (default: 'v0.1.1')
* ``overwrite`` (bool): If True, overwrite existing files (default: False)
* ``verbose`` (bool): Print download progress (default: True)

**Returns:**

* ``dict``: Dictionary with download results containing:
  
  * ``downloaded`` (list): List of successfully downloaded files
  * ``skipped`` (list): List of skipped (already existing) files
  * ``failed`` (list): List of failed downloads

**Downloaded Files:**

* ``test.bed``: PLINK binary genotype file (3,925 samples × 1,000 variants)
* ``test.bim``: PLINK variant information file
* ``test.fam``: PLINK sample information file
* ``test.phen``: Phenotype file with disease status and age

**Examples:**

.. code-block:: python

   # Example 1: Basic usage - download to default 'tests' directory
   from edge_gwas import download_test_files
   
   download_test_files()
   # Files saved to: tests/test.bed, tests/test.bim, tests/test.fam, tests/test.phen

   # Example 2: Download to custom directory
   download_test_files(output_dir='data/examples')

   # Example 3: Force re-download (overwrite existing files)
   results = download_test_files(overwrite=True)
   print(f"Downloaded: {results['downloaded']}")
   print(f"Skipped: {results['skipped']}")

   # Example 4: Download specific version
   download_test_files(version='main')  # Download from main branch

   # Example 5: Verify download results
   results = download_test_files()
   if results['failed']:
       print(f"Warning: Failed to download {results['failed']}")
   else:
       print("✓ All test files ready!")

**Complete Tutorial Workflow:**

.. code-block:: python

   from edge_gwas import (
       download_test_files,
       load_plink_data,
       EDGEAnalysis
   )
   import pandas as pd
   
   # Step 1: Download test data
   download_test_files()
   
   # Step 2: Load test data
   geno, info = load_plink_data('tests/test.bed', 'tests/test.bim', 'tests/test.fam')
   pheno = pd.read_csv('tests/test.phen', sep=' ')
   
   # Step 3: Prepare data
   pheno = pheno.set_index('IID')
   
   # Step 4: Split and analyze
   from edge_gwas.utils import stratified_train_test_split
   
   train_g, test_g, train_p, test_p = stratified_train_test_split(
       geno, pheno, 'disease', test_size=0.5
   )
   
   # Step 5: Run EDGE analysis
   edge = EDGEAnalysis(outcome_type='binary')
   alpha_df, gwas_df = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       outcome_col='disease',
       covariate_cols=['age']
   )
   
   print(f"✓ Analysis complete! Found {len(gwas_df)} variants")

**File Verification:**

.. code-block:: python

   import os
   
   # Download files
   results = download_test_files()
   
   # Verify file sizes and existence
   test_files = {
       'tests/test.bed': 'Genotype data (binary)',
       'tests/test.bim': 'Variant information',
       'tests/test.fam': 'Sample information',
       'tests/test.phen': 'Phenotype data'
   }
   
   for filepath, description in test_files.items():
       if os.path.exists(filepath):
           size_kb = os.path.getsize(filepath) / 1024
           print(f"✓ {filepath}")
           print(f"  {description} - {size_kb:.1f} KB")
       else:
           print(f"✗ {filepath} - MISSING!")

**Re-downloading Files:**

.. code-block:: python

   # Method 1: Use overwrite parameter
   download_test_files(overwrite=True)
   
   # Method 2: Delete directory first
   import shutil
   if os.path.exists('tests'):
       shutil.rmtree('tests')
   download_test_files()
   
   # Method 3: Download to new location
   download_test_files(output_dir='tests_fresh')

**Troubleshooting:**

If download fails:

.. code-block:: python

   # Check internet connection and retry
   results = download_test_files(overwrite=True, verbose=True)
   
   # If still failing, download manually:
   # 1. Visit: https://github.com/nicenzhou/edge-gwas/tree/v0.1.1/tests
   # 2. Download files: test.bed, test.bim, test.fam, test.phen
   # 3. Place in 'tests' directory

**Note:**

* Requires internet connection to download from GitHub
* Files are downloaded from the specified GitHub release version
* Total download size: ~1-2 MB
* Test data contains simulated genotypes for demonstration purposes only

**See Also:**

* :ref:`quickstart` - Quick start tutorial using test data
* :ref:`tutorials` - Detailed tutorials with test data examples
* ``load_plink_data()`` - Load the downloaded PLINK files
* ``prepare_phenotype_data()`` - Load the downloaded phenotype file

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

* Index: Sample IDs (as strings)
* Columns: Variant IDs
* Values: Genotype encoding (0, 1, 2, or NaN for missing)

  * 0 = Reference homozygote (REF/REF)
  * 1 = Heterozygote (REF/ALT)
  * 2 = Alternate homozygote (ALT/ALT)
  * NaN = Missing
  * For dosages: continuous values 0-2

**Example:**

.. code-block:: python

   #              rs123  rs456  rs789
   # Sample1        0      1      2
   # Sample2        1      2      0
   # Sample3        2      0      1

Phenotype DataFrame Format
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Expected format:**

* Index: Sample IDs (as strings, matching genotype data)
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

PCA DataFrame Format (NEW in v0.1.1)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Output format from PCA functions:**

.. code-block:: python

   #             PC1    PC2    PC3    PC4    PC5
   # IID                                        
   # Sample1   0.12  -0.23   0.45  -0.12   0.34
   # Sample2  -0.34   0.56  -0.12   0.23  -0.45
   # Sample3   0.23  -0.12   0.34  -0.45   0.12

GRM Sample IDs Format (NEW in v0.1.1)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Output format from load_grm_gcta():**

.. code-block:: python

   #      FID      IID
   # 0  FAM001  Sample1
   # 1  FAM001  Sample2
   # 2  FAM002  Sample3

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
   DEFAULT_CALL_RATE = 0.95       # Minimum sample call rate
   
   # Analysis defaults
   DEFAULT_TEST_SIZE = 0.5        # Train/test split ratio
   DEFAULT_N_JOBS = -1            # Use all CPU cores
   DEFAULT_MAX_ITER = 1000        # Maximum iterations
   
   # GRM defaults (NEW in v0.1.1)
   DEFAULT_GRM_MAF = 0.01         # MAF for GRM calculation
   DEFAULT_KIN_THRESHOLD = 0.0884 # Kinship threshold (3rd degree)
   
   # PCA defaults (NEW in v0.1.1)
   DEFAULT_N_PCS = 10             # Number of PCs
   DEFAULT_LD_R2 = 0.2            # LD pruning threshold
   DEFAULT_LD_WINDOW = 50         # LD window size (kb)

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
   
   class GRMError(EdgeGWASError):
       """Raised when GRM-related operations fail."""
       pass
   
   class ToolNotFoundError(EdgeGWASError):
       """Raised when external tool (PLINK2, GCTA, R) is not found."""
       pass

Version Information
-------------------

.. code-block:: python

   import edge_gwas
   
   # Get version string
   print(edge_gwas.__version__)  # '0.1.1'
   
   # Get author information
   print(edge_gwas.__author__)   # 'Your Name'
   
   # Get license
   print(edge_gwas.__license__)  # 'GPL-3.0'

Command-Line Tools (NEW in v0.1.1)
----------------------------------

edge-gwas-install-tools
~~~~~~~~~~~~~~~~~~~~~~~

Interactive installer for external tools (PLINK2, GCTA, R packages).

.. code-block:: bash

   # Install all tools interactively
   edge-gwas-install-tools
   
   # The installer will:
   # 1. Detect your operating system and architecture
   # 2. Download and install PLINK2
   # 3. Download and install GCTA
   # 4. Install R packages (GENESIS, SNPRelate, gdsfmt)
   # 5. Configure PATH if needed

**Supported platforms:**

* Linux (x86_64)
* macOS (Intel x86_64 and Apple Silicon ARM64)

**Features:**

* Automatic architecture detection
* SSL-safe downloads with retry logic
* Verification of installed tools
* PATH configuration guidance

edge-gwas-check-tools
~~~~~~~~~~~~~~~~~~~~~

Verify that external tools are properly installed.

.. code-block:: bash

   # Check all tools
   edge-gwas-check-tools
   
   # Output example:
   # ======================================================================
   # EDGE-GWAS External Tools Check
   # ======================================================================
   # 
   # Python Packages:
   # ----------------------------------------------------------------------
   # ✓ numpy: Installed
   # ✓ pandas: Installed
   # ✓ scipy: Installed
   # ✓ statsmodels: Installed
   # ✓ sklearn: Installed
   # ✓ matplotlib: Installed
   # ✓ pandas_plink: Installed
   # 
   # External Tools:
   # ----------------------------------------------------------------------
   # ✓ PLINK2: PLINK v2.00a3.7 64-bit Intel (24 Jan 2024)
   # ✓ GCTA: GCTA 1.95.0
   # 
   # R and Packages:
   # ----------------------------------------------------------------------
   # ✓ R: R version 4.3.0
   # ✓ R package GENESIS: Installed
   # ✓ R package SNPRelate: Installed
   # ✓ R package gdsfmt: Installed
   # 
   # ======================================================================
   # ✓ All tools and packages are properly installed!
   # ======================================================================

Function Reference Summary
--------------------------

Core Module (edge_gwas.core)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Class: EDGEAnalysis**

.. list-table::
   :header-rows: 1
   :widths: 20 30 25 25

   * - Method
     - Purpose
     - Key Input
     - Output
   * - ``calculate_alpha()``
     - Calculate encoding parameters from training data
     - Genotype, phenotype, covariates, GRM (optional)
     - Alpha DataFrame
   * - ``apply_alpha()``
     - Apply alpha values to test data for GWAS
     - Genotype, phenotype, alpha values, GRM (optional)
     - GWAS results DataFrame
   * - ``run_full_analysis()``
     - Complete two-stage workflow
     - Train/test genotype & phenotype, GRM (optional)
     - (alpha_df, gwas_df)
   * - ``get_skipped_snps()``
     - Get list of SNPs that failed convergence
     - None
     - List of variant IDs

Utilities Module (edge_gwas.utils)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Population Structure Control (NEW in v0.1.1):**

.. list-table::
   :header-rows: 1
   :widths: 25 35 40

   * - Function
     - Purpose
     - Output
   * - ``calculate_grm_gcta()``
     - Calculate GRM using GCTA
     - GRM file prefix
   * - ``load_grm_gcta()``
     - Load GRM from GCTA files
     - (grm_matrix, sample_ids)
   * - ``identify_related_samples()``
     - Find related sample pairs
     - DataFrame with related pairs
   * - ``filter_related_samples()``
     - Remove related samples
     - Filtered phenotype DataFrame

**Principal Component Analysis (NEW in v0.1.1):**

.. list-table::
   :header-rows: 1
   :widths: 25 35 20 20

   * - Function
     - Purpose
     - Output
     - Requirements
   * - ``calculate_pca_plink()``
     - PCA using PLINK2
     - PCA DataFrame
     - PLINK2
   * - ``calculate_pca_pcair()``
     - PC-AiR for related samples
     - PCA DataFrame
     - R + GENESIS
   * - ``calculate_pca_sklearn()``
     - Basic PCA
     - PCA DataFrame
     - scikit-learn
   * - ``attach_pcs_to_phenotype()``
     - Merge PCs with phenotype
     - Phenotype DataFrame with PCs
     - \-
   * - ``get_pc_covariate_list()``
     - Generate PC covariate names
     - List of PC names
     - \-

**Data Loading Functions:**

.. list-table::
   :header-rows: 1
   :widths: 20 25 20 20 15

   * - Function
     - Purpose
     - Key Input
     - Output
     - Requirements
   * - ``load_plink_data()``
     - Load PLINK binary files
     - .bed, .bim, .fam paths
     - (genotype_df, variant_info_df)
     - pandas-plink
   * - ``load_pgen_data()``
     - Load PLINK 2 binary files
     - .pgen, .pvar, .psam paths
     - (genotype_df, variant_info_df)
     - pgenlib
   * - ``load_bgen_data()``
     - Load BGEN files
     - .bgen, .sample paths
     - (genotype_df, variant_info_df)
     - bgen-reader
   * - ``load_vcf_data()``
     - Load VCF/VCF.GZ files
     - .vcf/.vcf.gz path
     - (genotype_df, variant_info_df)
     - cyvcf2
   * - ``prepare_phenotype_data()``
     - Load and prepare phenotype data
     - File path, column names
     - Phenotype DataFrame
     - \-

**Data Processing Functions:**

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Function
     - Purpose
     - Output
   * - ``stratified_train_test_split()``
     - Split data into train/test sets
     - (train_g, test_g, train_p, test_p)
   * - ``filter_variants_by_maf()``
     - Filter by minor allele frequency
     - Filtered DataFrame
   * - ``filter_variants_by_missing()``
     - Filter by missingness rate
     - Filtered DataFrame
   * - ``filter_variants_by_hwe()``
     - Filter by HWE p-value
     - Filtered DataFrame
   * - ``filter_samples_by_call_rate()``
     - Filter samples by call rate
     - (filtered_geno, filtered_pheno)
   * - ``merge_alpha_with_gwas()``
     - Merge GWAS and alpha results
     - Merged DataFrame

**Statistical Functions:**

.. list-table::
   :header-rows: 1
   :widths: 30 40 30

   * - Function
     - Purpose
     - Output
   * - ``calculate_genomic_inflation()``
     - Calculate lambda (λ)
     - Lambda (float)
   * - ``calculate_hwe_pvalues()``
     - Calculate HWE p-values
     - HWE p-values Series
   * - ``check_case_control_balance()``
     - Check case/control balance
     - Dictionary with counts and ratio

**Analysis Functions (NEW in v0.1.1):**

.. list-table::
   :header-rows: 1
   :widths: 35 35 30

   * - Function
     - Purpose
     - Output
   * - ``additive_gwas()``
     - Standard additive GWAS for comparison
     - GWAS results DataFrame
   * - ``cross_validated_edge_analysis()``
     - K-fold cross-validation for EDGE
     - (avg_alpha, meta_gwas, all_alpha, all_gwas)

Complete Workflow Example
--------------------------

Here's a complete example using new v0.1.1 features:

.. code-block:: python

   from edge_gwas import EDGEAnalysis
   from edge_gwas.utils import (
       load_plink_data,
       prepare_phenotype_data,
       stratified_train_test_split,
       filter_variants_by_maf,
       calculate_grm_gcta,
       load_grm_gcta,
       calculate_pca_plink,
       attach_pcs_to_phenotype,
       get_pc_covariate_list
   )
   from edge_gwas.visualize import manhattan_plot, qq_plot
   
   # 1. Load data
   geno_df, var_info = load_plink_data('data.bed', 'data.bim', 'data.fam')
   pheno_df = prepare_phenotype_data(
       'phenotype.txt',
       outcome_col='disease',
       covariate_cols=['age', 'sex'],
       sample_id_col='IID'
   )
   
   # 2. QC filtering
   geno_df = filter_variants_by_maf(geno_df, min_maf=0.01)
   
   # 3. Calculate GRM for population structure control
   grm_prefix = calculate_grm_gcta('data', output_prefix='grm/output')
   grm_matrix, grm_ids = load_grm_gcta('grm/output')
   
   # 4. Calculate PCA
   pca_df = calculate_pca_plink('data', n_pcs=10)
   pheno_df = attach_pcs_to_phenotype(pheno_df, pca_df, n_pcs=10)
   
   # 5. Split data
   train_g, test_g, train_p, test_p = stratified_train_test_split(
       geno_df, pheno_df, outcome_col='disease', test_size=0.5
   )
   
   # 6. Run EDGE analysis with outcome transformation and GRM
   edge = EDGEAnalysis(
       outcome_type='continuous',
       outcome_transform='rank_inverse_normal'
   )
   
   covariates = ['age', 'sex'] + get_pc_covariate_list(10)
   
   alpha_df, gwas_df = edge.run_full_analysis(
       train_g, train_p,
       test_g, test_p,
       outcome='trait',
       covariates=covariates,
       grm_matrix=grm_matrix,
       grm_sample_ids=grm_ids,
       output_prefix='results/edge'
   )
   
   # 7. Visualize results
   manhattan_plot(gwas_df, output='results/manhattan.png')
   lambda_gc = qq_plot(gwas_df, output='results/qq.png')
   
   print(f"Analysis complete. Lambda GC: {lambda_gc:.3f}")

Migration Guide (v0.1.0 → v0.1.1)
----------------------------------

**Breaking Changes:**

1. **Koalas removed**: Replace ``.to_koalas()`` with pandas

   .. code-block:: python
   
      # Old (v0.1.0)
      import databricks.koalas as ks
      df = data.to_koalas()
      
      # New (v0.1.1)
      import pandas as pd
      df = data  # Already pandas

2. **PCA functions return indexed DataFrames**: All PCA functions now return DataFrames with IID as index

**New Features to Adopt:**

1. **Use GRM for population structure**:

   .. code-block:: python
   
      grm_matrix, grm_ids = load_grm_gcta('grm_prefix')
      alpha_df = edge.calculate_alpha(..., grm_matrix=grm_matrix, grm_sample_ids=grm_ids)

2. **Use outcome transformations**:

   .. code-block:: python
   
      edge = EDGEAnalysis(outcome_type='continuous', outcome_transform='rank_inverse_normal')

3. **Add PCs as covariates**:

   .. code-block:: python
   
      pca_df = calculate_pca_plink('genotypes', n_pcs=10)
      pheno_df = attach_pcs_to_phenotype(pheno_df, pca_df, n_pcs=10)
      covariates = ['age', 'sex'] + get_pc_covariate_list(10)

Index of All Functions
----------------------

**Alphabetical listing:**

* ``additive_gwas()`` - Standard additive GWAS
* ``apply_alpha()`` - Apply EDGE encoding to test data
* ``attach_pcs_to_phenotype()`` - Merge PCs with phenotype
* ``calculate_alpha()`` - Calculate EDGE encoding parameters
* ``calculate_genomic_inflation()`` - Calculate lambda GC
* ``calculate_grm_gcta()`` - Calculate GRM with GCTA
* ``calculate_hwe_pvalues()`` - Calculate HWE p-values
* ``calculate_pca_pcair()`` - PC-AiR for related samples
* ``calculate_pca_plink()`` - PCA using PLINK2
* ``calculate_pca_sklearn()`` - Basic PCA with scikit-learn
* ``check_case_control_balance()`` - Check case/control ratio
* ``create_summary_report()`` - Generate text summary
* ``cross_validated_edge_analysis()`` - K-fold cross-validation
* ``filter_related_samples()`` - Remove related samples
* ``filter_samples_by_call_rate()`` - Filter by sample call rate
* ``filter_variants_by_hwe()`` - Filter by HWE p-value
* ``filter_variants_by_maf()`` - Filter by minor allele frequency
* ``filter_variants_by_missing()`` - Filter by missingness
* ``get_pc_covariate_list()`` - Generate PC covariate names
* ``get_skipped_snps()`` - Get list of failed variants
* ``identify_related_samples()`` - Find related pairs
* ``load_alpha_values()`` - Load pre-calculated alphas
* ``load_bgen_data()`` - Load BGEN format
* ``load_grm_gcta()`` - Load GRM from GCTA files
* ``load_pgen_data()`` - Load PLINK2 PGEN format
* ``load_plink_data()`` - Load PLINK binary format
* ``load_vcf_data()`` - Load VCF/VCF.GZ format
* ``manhattan_plot()`` - Create Manhattan plot
* ``merge_alpha_with_gwas()`` - Merge alpha and GWAS results
* ``plot_alpha_distribution()`` - Plot alpha distribution
* ``prepare_phenotype_data()`` - Load and prepare phenotypes
* ``qq_plot()`` - Create QQ plot
* ``run_full_analysis()`` - Complete two-stage EDGE workflow
* ``save_results()`` - Save GWAS and alpha results
* ``stratified_train_test_split()`` - Train/test split with stratification

**By Category:**

*Core Analysis:*

* ``EDGEAnalysis`` (class)

  * ``calculate_alpha()``
  * ``apply_alpha()``
  * ``run_full_analysis()``
  * ``get_skipped_snps()``

*Data Loading:*

* ``load_plink_data()``
* ``load_pgen_data()`` (NEW in v0.1.1)
* ``load_bgen_data()`` (NEW in v0.1.1)
* ``load_vcf_data()`` (NEW in v0.1.1)
* ``prepare_phenotype_data()``

*Population Structure (NEW in v0.1.1):*

* ``calculate_grm_gcta()``
* ``load_grm_gcta()``
* ``identify_related_samples()``
* ``filter_related_samples()``
* ``calculate_pca_plink()``
* ``calculate_pca_pcair()``
* ``calculate_pca_sklearn()``
* ``attach_pcs_to_phenotype()``
* ``get_pc_covariate_list()``

*Quality Control:*

* ``filter_variants_by_maf()``
* ``filter_variants_by_missing()``
* ``filter_variants_by_hwe()`` (NEW in v0.1.1)
* ``filter_samples_by_call_rate()`` (NEW in v0.1.1)
* ``check_case_control_balance()`` (NEW in v0.1.1)
* ``calculate_hwe_pvalues()`` (NEW in v0.1.1)

*Statistical Functions:*

* ``calculate_genomic_inflation()``
* ``merge_alpha_with_gwas()``
* ``additive_gwas()`` (NEW in v0.1.1)
* ``cross_validated_edge_analysis()`` (NEW in v0.1.1)

*Data Processing:*

* ``stratified_train_test_split()``

*Visualization:*

* ``manhattan_plot()``
* ``qq_plot()``
* ``plot_alpha_distribution()``

*I/O:*

* ``save_results()``
* ``load_alpha_values()``
* ``create_summary_report()``

**Command-Line Tools (NEW in v0.1.1):**

* ``edge-gwas-install-tools`` - Install external tools
* ``edge-gwas-check-tools`` - Verify tool installation

Quick Function Finder
----------------------

**"I want to..."**

* **Load genetic data**:

  * PLINK binary (.bed/.bim/.fam): ``load_plink_data()``
  * PLINK2 (.pgen/.pvar/.psam): ``load_pgen_data()``
  * BGEN: ``load_bgen_data()``
  * VCF: ``load_vcf_data()``

* **Control for population structure**:

  * Calculate GRM: ``calculate_grm_gcta()``
  * Load GRM: ``load_grm_gcta()``
  * Calculate PCs (unrelated): ``calculate_pca_plink()``
  * Calculate PCs (related): ``calculate_pca_pcair()``
  * Add PCs to phenotype: ``attach_pcs_to_phenotype()``

* **Quality control**:

  * Filter by MAF: ``filter_variants_by_maf()``
  * Filter by missingness: ``filter_variants_by_missing()``
  * Filter by HWE: ``filter_variants_by_hwe()``
  * Filter samples: ``filter_samples_by_call_rate()``
  * Check case/control balance: ``check_case_control_balance()``
  * Find related samples: ``identify_related_samples()``
  * Remove related samples: ``filter_related_samples()``

* **Run EDGE analysis**:

  * Calculate alpha: ``calculate_alpha()``
  * Apply alpha: ``apply_alpha()``
  * Full workflow: ``run_full_analysis()``
  * Cross-validation: ``cross_validated_edge_analysis()``

* **Compare with standard GWAS**:

  * Run additive model: ``additive_gwas()``

* **Visualize results**:

  * Manhattan plot: ``manhattan_plot()``
  * QQ plot: ``qq_plot()``
  * Alpha distribution: ``plot_alpha_distribution()``

* **Save/load results**:

  * Save results: ``save_results()``
  * Load alpha values: ``load_alpha_values()``
  * Create summary: ``create_summary_report()``

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
