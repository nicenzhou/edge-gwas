.. _changelog:

Changelog
=========

All notable changes to edge-gwas will be documented in this file.

The format is based on `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_,
and this project adheres to `Semantic Versioning <https://semver.org/spec/v2.0.0.html>`_.

Version 0.1.1 (2025-12-25) - Current Release
--------------------------------------------

**Status:** ⚠️ Public Testing Phase

Changed
^^^^^^^

**Data Handling:**

* **Breaking Change:** Replaced Koalas DataFrames with pandas DataFrames throughout codebase
* Updated ``load_participant_data()`` to use pandas instead of ``participant.retrieve_fields().to_koalas()``
* Simplified DataFrame operations by removing Koalas dependency
* Improved compatibility with standard pandas ecosystem
* Better performance for typical dataset sizes

**Dependencies:**

* Removed dependency on ``databricks.koalas``
* All operations now use native pandas (>= 1.2.0)

Migration from v0.1.0
^^^^^^^^^^^^^^^^^^^^^

**Breaking Changes:**

* Code using ``.to_koalas()`` will need to be updated to use pandas DataFrames
* Remove any ``import databricks.koalas`` statements
* Replace ``koalas`` DataFrame operations with standard pandas operations

**Migration Example:**

.. code-block:: python

   # Old code (v0.1.0)
   import databricks.koalas as ks
   participant_df = participant.retrieve_fields(fields=fields, engine=dxdata.connect()).to_koalas()
   
   # New code (v0.1.1)
   import pandas as pd
   participant_df = participant.retrieve_fields(fields=fields, engine=dxdata.connect())

Bug Fixes
^^^^^^^^^

* Resolved compatibility issues with recent pandas versions
* Improved memory efficiency in data loading operations

Version 0.1.0 (2025-12-24)
--------------------------

**Status:** Public Testing Phase

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
* ``merge_alpha_with_gwas()`` - Merge alpha and GWAS results

**Statistical Functions:**

* ``calculate_genomic_inflation()`` - Calculate genomic inflation factor (λ)
* Two-stage analysis framework (calculate alpha → apply alpha)
* Parallel processing support with configurable CPU cores
* Convergence monitoring and skipped SNP tracking

**Visualization:**

* ``manhattan_plot()`` - Create Manhattan plots with customization
* ``qq_plot()`` - Create QQ plots with lambda calculation
* ``plot_alpha_distribution()`` - Visualize alpha value distribution
* Multi-chromosome support for all plots
* Customizable colors, thresholds, and figure sizes

**I/O Handlers:**

* ``save_results()`` - Save GWAS and alpha results
* ``load_alpha_values()`` - Load pre-calculated alpha values
* ``create_summary_report()`` - Generate text summary reports
* Standardized output formats (tab-delimited text)

**Documentation:**

* Complete Sphinx documentation with Read the Docs integration
* Installation guide with multiple methods
* Quick start guide with minimal example
* Comprehensive user guide
* Statistical model documentation with mathematical equations
* API reference with all functions and parameters
* Example workflows covering common use cases
* Visualization guide with custom plotting examples

**Quality Control:**

* MAF filtering (default: 0.01)
* Missingness filtering (default: 0.05)
* Hardy-Weinberg equilibrium support
* Sample call rate filtering
* Automatic quality checks during analysis

**Performance:**

* Multi-core parallel processing (``n_jobs`` parameter)
* Memory-efficient chromosome-by-chromosome processing
* Batch processing support for large datasets
* Progress tracking with ``verbose`` option

Features
^^^^^^^^

* **Flexible Encoding:** Detects recessive, additive, dominant, and over-dominant effects
* **Two-Stage Design:** Reduces overfitting through train/test split
* **PLINK Integration:** Native support for standard genomic file formats
* **Automated QC:** Built-in quality control filters
* **Rich Visualization:** Publication-quality plots
* **Comprehensive Documentation:** Full API reference and examples

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

Installation Methods
^^^^^^^^^^^^^^^^^^^^

* GitHub installation: ``pip install git+https://github.com/YOUR-USERNAME/edge-gwas.git``
* Development installation: ``pip install -e .``
* Virtual environment support (venv, conda)

Known Issues
^^^^^^^^^^^^

* Package is under active development and testing
* Some edge cases in convergence may require manual inspection
* Large datasets may require batch processing for memory efficiency
* Documentation is being continuously expanded
* Koalas dependency will be removed in v0.1.1

Migration from v0.0.0
^^^^^^^^^^^^^^^^^^^^^

Users of the original EDGE implementation (v0.0.0, standalone scripts) should:

1. Uninstall or archive old scripts
2. Install edge-gwas v0.1.0 via pip
3. Update code to use new API (see :ref:`quickstart`)
4. Review new documentation for best practices

**Breaking Changes from v0.0.0:**

* New package structure (not compatible with old scripts)
* Different function names and signatures
* New data format requirements
* Updated output file formats

Notes
^^^^^

* This is a **testing release** - please report issues on GitHub
* Version 0.0.0 (standalone scripts) is no longer maintained
* See https://github.com/nicenzhou/EDGE for legacy version
* Recommended to use v0.1.1 or later for new projects

---

Version 0.0.0 (2024-04-02) - Legacy
------------------------------------

**Status:** 🔒 Closed - No Further Maintenance

**Repository:** https://github.com/nicenzhou/EDGE

Legacy Release
~~~~~~~~~~~~~~

Original EDGE implementation as standalone Python scripts.

Features
^^^^^^^^

* Basic EDGE analysis implementation
* Individual Python scripts for each analysis step
* Manual workflow execution
* Support for binary and quantitative outcomes
* Basic visualization functions

Limitations
^^^^^^^^^^^

* No package structure
* Manual dependency management
* No automated installation
* Limited documentation
* No quality control automation
* Requires manual file path configuration

Deprecation Notice
^^^^^^^^^^^^^^^^^^

.. warning::
   Version 0.0.0 is **deprecated** and no longer maintained.
   
   * No bug fixes will be provided
   * No new features will be added
   * Users should migrate to v0.1.1 or later
   * Legacy repository remains available for reference only

Migration Guide
^^^^^^^^^^^^^^^

To migrate from v0.0.0 to v0.1.1:

**Old code (v0.0.0):**

.. code-block:: python

   # Manual script execution
   python calculate_alpha.py --input data.txt --output alpha.txt
   python apply_alpha.py --alpha alpha.txt --input test.txt

**New code (v0.1.1):**

.. code-block:: python

   from edge_gwas import EDGEAnalysis
   
   edge = EDGEAnalysis(outcome_type='binary')
   alpha_df, gwas_df = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       outcome='disease', covariates=['age', 'sex']
   )

---

Upcoming Features (Roadmap)
----------------------------

Version 0.2.0 (Planned)
~~~~~~~~~~~~~~~~~~~~~~~

Enhancements
^^^^^^^^^^^^

* Additional statistical tests and corrections
* Support for more genetic data formats
* Enhanced parallel processing capabilities
* Advanced filtering options
* Interactive visualization tools

Bug Fixes
^^^^^^^^^

* Improved convergence handling
* Better error messages and validation
* Memory optimization for large datasets

---

Version History Summary
-----------------------

.. list-table::
   :header-rows: 1
   :widths: 15 15 20 50

   * - 0.1.1
     - 2025-12-25
     - Current
     - Pandas migration, removed Koalas dependency
   * - 0.1.0
     - 2025-12-24
     - Superseded
     - First packaged release, full API, documentation
   * - 0.0.0
     - 2024-04-02
     - Deprecated
     - Original standalone scripts (no longer maintained)

---

See Also
--------

* :ref:`installation` - Installation instructions
* :ref:`quickstart` - Getting started guide
* GitHub Releases - https://github.com/YOUR-USERNAME/edge-gwas/releases
* Issue Tracker - https://github.com/YOUR-USERNAME/edge-gwas/issues
