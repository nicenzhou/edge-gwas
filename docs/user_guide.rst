.. _user_guide:

User Guide
==========

This comprehensive guide covers all aspects of using edge-gwas.

Overview
--------

edge-gwas is designed for researchers performing large-scale GWAS analysis
with modern computational infrastructure.

Core Concepts
-------------

GWAS Analysis
~~~~~~~~~~~~~

Genome-Wide Association Studies (GWAS) identify genetic variants associated
with traits or diseases.

Edge Computing
~~~~~~~~~~~~~~

edge-gwas leverages edge computing for distributed processing of genomic data.

Data Formats
------------

Supported Input Formats
~~~~~~~~~~~~~~~~~~~~~~~

* **PLINK** - Binary PED/MAP files
* **TXT/CSV/TSV** - Phenotype data

Output Formats
~~~~~~~~~~~~~~

* **CSV** - Association results
* **PNG/PDF** - Visualization outputs

Workflow Components
-------------------

1. Data Loading
~~~~~~~~~~~~~~~

Load and prepare genomic and phenotype data.

2. Quality Control
~~~~~~~~~~~~~~~~~~

Filter variants and samples based on:

* Minor Allele Frequency (MAF)
* Genotyping call rate

3. Association Testing
~~~~~~~~~~~~~~~~~~~~~~

Perform statistical tests:

* Linear regression (quantitative traits)
* Logistic regression (binary traits)

4. Results Interpretation
~~~~~~~~~~~~~~~~~~~~~~~~~

* P-value thresholds
* Multiple testing correction
* Effect size interpretation

Best Practices
--------------

Quality Control
~~~~~~~~~~~~~~~

1. Filter variants with MAF < 0.01

Statistical Analysis
~~~~~~~~~~~~~~~~~~~~

1. Always include covariates (age, sex, PCs)
2. Use an appropriate multiple testing correction
3. Validate top hits in independent cohorts
4. Check for confounding factors

Performance Optimization
~~~~~~~~~~~~~~~~~~~~~~~~

1. Use binary file formats when possible
2. Filter data before analysis
3. Utilize parallel processing
4. Monitor memory usage

Advanced Topics
---------------

Meta-Analysis
~~~~~~~~~~~~~

Combine results from multiple GWAS studies.

Troubleshooting
---------------

Common Issues
~~~~~~~~~~~~~

**Memory errors**: Reduce batch size or filter more aggressively

**Slow performance**: Enable parallel processing

**Missing values**: Check data preprocessing

Getting Help
------------

* Check the :ref:`api` documentation
* Review :ref:`examples`
* Open an issue on GitHub

References
----------

For more information on GWAS methodology:

* Bush WS, Moore JH. Chapter 11: Genome-wide association studies. PLoS Comput Biol. 2012.
* Visscher PM, et al. 10 Years of GWAS Discovery. Am J Hum Genet. 2017.
