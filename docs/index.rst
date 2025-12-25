Welcome to edge-gwas Documentation
===================================

**EDGE-GWAS** (Elastic Data-Driven Encoding GWAS) identifies nonadditive SNP effects 
using flexible genetic encoding, rather than assuming additive inheritance.

.. warning::
   **Version 0.1.0** is currently under active development and public testing.

.. note::
   The original EDGE implementation (v0.0.0) is available at 
   `https://github.com/nicenzhou/EDGE <https://github.com/nicenzhou/EDGE>`_.
   Version 0.0.0 is **no longer maintained**; users are encouraged to migrate to v0.1.0.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   quickstart
   user_guide
   api_reference
   examples
   statistical_model
   visualization
   changelog
   contributing
   citation

Key Features
------------

* **Two-stage analysis**: Calculate alpha on training data, apply to test data
* **Flexible encoding**: Detects under-recessive, recessive, additive, dominant, and over-dominant effects
* **Multiple outcomes**: Handles binary and quantitative outcomes
* **PLINK support**: Native support for PLINK format data
* **Built-in visualization**: Manhattan plots, QQ plots, and alpha distribution plots

Quick Example
-------------

.. code-block:: python

   from edge_gwas import EDGEAnalysis
   from edge_gwas.utils import load_plink_data, prepare_phenotype_data
   
   # Load data
   geno, info = load_plink_data('data.bed', 'data.bim', 'data.fam')
   pheno = prepare_phenotype_data('pheno.txt', 'disease', ['age', 'sex'])
   
   # Run EDGE analysis
   edge = EDGEAnalysis(outcome_type='binary')
   alpha_df, gwas_df = edge.run_full_analysis(
       train_geno, train_pheno, test_geno, test_pheno,
       outcome='disease', covariates=['age', 'sex']
   )

Getting Help
------------

* `GitHub Issues <https://github.com/YOUR-USERNAME/edge-gwas/issues>`_
* `GitHub Discussions <https://github.com/YOUR-USERNAME/edge-gwas/discussions>`_
* Email: jyzhou@stanford.edu (Code Questions)
* Email: molly.hall@pennmedicine.upenn.edu (Research Questions)

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
