
Welcome to edge-gwas Documentation
===================================

**EDGE-GWAS** (Elastic Data-Driven Encoding GWAS) identifies nonadditive SNP effects 
using flexible genetic encoding, rather than assuming additive inheritance.

.. warning::
   ⚠️ **Current Version 0.1.1** - Under Public Testing
   
   **Recommended to use v0.1.1 - more stable and more functions.**

.. note::
   The original EDGE implementation (v0.0.0) is available at 
   `https://github.com/nicenzhou/EDGE <https://github.com/nicenzhou/EDGE>`_.
   Version 0.0.0 is **no longer maintained**; users are encouraged to migrate to v0.1.0+.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   quickstart
   user_guide
   api_reference
   statistical_model
   examples
   visualization
   changelog
   citation
   troubleshooting
   faq
   futureupdates

Key Features
------------

* **Two-stage analysis**: Calculate alpha on training data, apply to test data
* **Flexible encoding**: Detects recessive, additive, dominant, and over-dominant effects
* **Multiple outcomes**: Binary and quantitative traits
* **PLINK support**: Native .bed/.bim/.fam format
* **Visualization**: Manhattan, QQ, and alpha distribution plots

Quick Start
-----------

.. code-block:: python

   from edge_gwas import EDGEAnalysis
   from edge_gwas.utils import load_plink_data, prepare_phenotype_data
   
   # Load data
   geno, info = load_plink_data('data.bed', 'data.bim', 'data.fam')
   pheno = prepare_phenotype_data('pheno.txt', 'disease', ['age', 'sex'])
   
   # Run analysis
   edge = EDGEAnalysis(outcome_type='binary')
   alpha_df, gwas_df = edge.run_full_analysis(
       train_geno, train_pheno, test_geno, test_pheno,
       outcome='disease', covariates=['age', 'sex']
   )

Support
-------

* **Issues**: `GitHub Issues <https://github.com/YOUR-USERNAME/edge-gwas/issues>`_
* **Discussions**: `GitHub Discussions <https://github.com/YOUR-USERNAME/edge-gwas/discussions>`_
* **Code Questions**: jyzhou@stanford.edu
* **Research Questions**: molly.hall@pennmedicine.upenn.edu

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
