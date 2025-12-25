.. _api:

API Reference
=============

This page contains the complete API reference for edge-gwas.

.. note::
   Full API implementation coming in v0.1.1.
   This documentation shows the planned API structure.

Main Modules
------------

edge_gwas.gwas
~~~~~~~~~~~~~~

Core GWAS analysis functions.

.. code-block:: python

   from edge_gwas import gwas
   
   # Main association testing function
   results = gwas.run_association(genotypes, phenotypes, covariates)

**Functions:**

* ``run_association()`` - Perform GWAS analysis
* ``linear_regression()`` - Linear regression for quantitative traits
* ``logistic_regression()`` - Logistic regression for binary traits
* ``mixed_model()`` - Mixed model analysis

edge_gwas.data_loader
~~~~~~~~~~~~~~~~~~~~~

Data loading and preprocessing utilities.

.. code-block:: python

   from edge_gwas import data_loader
   
   # Load genotype data
   genotypes = data_loader.load_genotypes('data.vcf')
   
   # Load phenotype data
   phenotypes = data_loader.load_phenotypes('phenotypes.csv')

**Functions:**

* ``load_genotypes()`` - Load genotype data (VCF, PLINK)
* ``load_phenotypes()`` - Load phenotype data
