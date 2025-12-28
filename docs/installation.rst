.. _installation:

Installation Guide
==================

System Requirements
-------------------

* Python 3.7+
* pip package manager
* 4GB RAM minimum (8GB+ recommended)
* 16GB+ RAM for biobank-scale data with GRM

Quick Installation
------------------

**Install edge-gwas:**

.. code-block:: bash

   pip install git+https://github.com/nicenzhou/edge-gwas.git

**Install external tools:**

.. code-block:: bash

   edge-gwas-install-tools

**Verify installation:**

.. code-block:: bash

   edge-gwas-check-tools

Virtual Environment (Recommended)
----------------------------------

.. code-block:: bash

   # Create environment
   python3 -m venv edge-gwas-env
   source edge-gwas-env/bin/activate
   
   # Install
   pip install git+https://github.com/nicenzhou/edge-gwas.git
   edge-gwas-install-tools

Optional Dependencies
---------------------

**File format support:**

.. code-block:: bash

   pip install pgenlib bgen-reader cyvcf2

Upgrading from v0.1.0
---------------------

.. code-block:: bash

   pip install --upgrade --force-reinstall git+https://github.com/nicenzhou/edge-gwas.git

**Breaking changes:**

* Koalas removed (now uses pandas)
* PCA functions return DataFrames with IID as index
* All sample IDs converted to strings

Troubleshooting
---------------

**PLINK2/GCTA not found:**

.. code-block:: bash

   edge-gwas-install-tools
   export PATH="$HOME/.local/bin:$PATH"

**Memory errors with GRM:**

.. code-block:: python

   from edge_gwas.utils import calculate_grm_gcta
   
   grm_prefix = calculate_grm_gcta('genotypes', method='grm-sparse')

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
