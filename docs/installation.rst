.. _installation:

Installation Guide
==================

System Requirements
-------------------

* Python 3.7 or higher
* pip package manager
* 4GB RAM minimum (8GB+ recommended for large datasets)

Dependencies
------------

edge-gwas requires the following packages:

* numpy >= 1.19.0
* pandas >= 1.2.0
* scipy >= 1.6.0
* statsmodels >= 0.12.0
* scikit-learn >= 0.24.0
* matplotlib >= 3.3.0
* pandas-plink >= 2.0.0

Installation from GitHub
------------------------

**Recommended method for v0.1.0:**

.. code-block:: bash

   pip install git+https://github.com/YOUR-USERNAME/edge-gwas.git

Installation on Mac/Linux
-------------------------

.. code-block:: bash

   # Clone the repository
   git clone https://github.com/YOUR-USERNAME/edge-gwas.git
   cd edge-gwas
   
   # Install dependencies
   pip install -r requirements.txt
   
   # Install the package in development mode
   pip install -e .

Virtual Environment Setup
-------------------------

**Using venv:**

.. code-block:: bash

   python3 -m venv edge-gwas-env
   source edge-gwas-env/bin/activate  # On Windows: edge-gwas-env\Scripts\activate
   pip install git+https://github.com/YOUR-USERNAME/edge-gwas.git

**Using conda:**

.. code-block:: bash

   conda create -n edge-gwas python=3.9
   conda activate edge-gwas
   pip install git+https://github.com/YOUR-USERNAME/edge-gwas.git

Verify Installation
-------------------

**In Python/Jupyter Notebook:**

.. code-block:: python

   from edge_gwas import EDGEAnalysis, manhattan_plot, qq_plot
   print("✓ edge-gwas installed successfully")

**Using command line (Mac/Linux):**

.. code-block:: bash

   python3 -c "from edge_gwas import EDGEAnalysis, manhattan_plot, qq_plot; print('✓ Installed successfully')"

Troubleshooting
---------------

**ImportError: No module named 'edge_gwas'**

Ensure virtual environment is activated and package is installed.

**Dependency conflicts**

Create a fresh environment:

.. code-block:: bash

   conda create -n edge-gwas-fresh python=3.9
   conda activate edge-gwas-fresh
   pip install git+https://github.com/YOUR-USERNAME/edge-gwas.git

**Permission errors**

Use ``--user`` flag:

.. code-block:: bash

   pip install --user git+https://github.com/YOUR-USERNAME/edge-gwas.git

Platform-Specific Notes
-----------------------

**Windows:**
- Use Command Prompt or PowerShell
- May require Visual C++ Build Tools for scientific packages

**macOS:**
- May need Xcode Command Line Tools: ``xcode-select --install``

**Linux:**
- Ubuntu/Debian: ``sudo apt-get install python3-dev``

Next Steps
----------

After installation:

* :ref:`quickstart` - Quick start guide
* :ref:`examples` - Example workflows
* :ref:`api_reference` - API documentation
