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
       outcome_type='binary',
       n_jobs=-1,
       max_iter=1000,
       verbose=True
   )

**Parameters:**

* ``outcome_type`` (str): Type of outcome variable - ``'binary'`` or ``'continuous'``
* ``n_jobs`` (int): Number of CPU cores (``-1`` uses all cores)
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

* ``genotype_data`` (pandas.DataFrame): Genotype matrix (samples by variants)
* ``phenotype_df`` (pandas.DataFrame): Phenotype data with outcome and covariates
* ``outcome`` (str): Name of outcome column
* ``covariates`` (list): List of covariate column names
* ``variant_info`` (pandas.DataFrame, optional): Variant information

**Returns:**

* ``pandas.DataFrame``: Alpha values with variant_id, alpha_value, eaf, coef_het, coef_hom

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

* ``pandas.DataFrame``: GWAS results with variant_id, chr, pos, pval, coef, std_err, stat, alpha_value

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

load_pgen_data()
""""""""""""""""

Load PLINK 2 binary format data.

.. code-block:: python

   from edge_gwas.utils import load_pgen_data
   
   genotype_df, variant_info = load_pgen_data(
       pgen_file='data.pgen',
       pvar_file='data.pvar',
       psam_file='data.psam',
       verbose=True
   )

**Note:** Requires ``pgenlib`` package

load_bgen_data()
""""""""""""""""

Load BGEN format data with dosages.

.. code-block:: python

   from edge_gwas.utils import load_bgen_data
   
   genotype_df, variant_info = load_bgen_data(
       bgen_file='data.bgen',
       sample_file='data.sample',
       verbose=True
   )

**Note:** Requires ``bgen-reader`` package

load_vcf_data()
"""""""""""""""

Load VCF format data.

.. code-block:: python

   from edge_gwas.utils import load_vcf_data
   
   genotype_df, variant_info = load_vcf_data(
       vcf_file='data.vcf.gz',
       dosage=True,
       verbose=True
   )

**Note:** Requires ``cyvcf2`` package

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

Data Processing Functions
~~~~~~~~~~~~~~~~~~~~~~~~~

stratified_train_test_split()
"""""""""""""""""""""""""""""

Split data into training and test sets.

filter_variants_by_maf()
""""""""""""""""""""""""

Filter variants by minor allele frequency.

filter_variants_by_missing()
""""""""""""""""""""""""""""

Filter variants by missingness rate.

filter_samples_by_call_rate()
"""""""""""""""""""""""""""""

Filter samples by genotype call rate.

filter_variants_by_hwe()
""""""""""""""""""""""""

Filter variants by Hardy-Weinberg Equilibrium p-value.

check_case_control_balance()
""""""""""""""""""""""""""""

Check case/control balance in binary outcome.

calculate_hwe_pvalues()
"""""""""""""""""""""""

Calculate Hardy-Weinberg Equilibrium p-values for each variant.

Statistical Functions
~~~~~~~~~~~~~~~~~~~~~

calculate_genomic_inflation()
"""""""""""""""""""""""""""""

Calculate genomic inflation factor.

merge_alpha_with_gwas()
"""""""""""""""""""""""

Merge GWAS results with alpha values.

additive_gwas()
"""""""""""""""

Perform standard additive GWAS for comparison with EDGE.

cross_validated_edge_analysis()
""""""""""""""""""""""""""""""""

Perform k-fold cross-validation for EDGE analysis.

Visualization Module
--------------------

manhattan_plot()
~~~~~~~~~~~~~~~~

Create Manhattan plot of GWAS results.

qq_plot()
~~~~~~~~~

Create QQ plot and calculate genomic inflation factor.

plot_alpha_distribution()
~~~~~~~~~~~~~~~~~~~~~~~~~

Plot distribution of alpha values.

I/O Handlers Module
-------------------

save_results()
~~~~~~~~~~~~~~

Save GWAS and alpha results to files.

load_alpha_values()
~~~~~~~~~~~~~~~~~~~

Load pre-calculated alpha values from file.

create_summary_report()
~~~~~~~~~~~~~~~~~~~~~~~

Generate text summary of GWAS results.

See Also
--------

* :ref:`quickstart` - Quick start guide
* :ref:`examples` - Example workflows
* :ref:`statistical_model` - Statistical methodology
* :ref:`user_guide` - Detailed user guide
