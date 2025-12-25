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

Installation Methods
--------------------

Install from GitHub (Recommended for v0.1.0)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pip install git+https://github.com/YOUR-USERNAME/edge-gwas.git

Install from PyPI (Coming Soon)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pip install edge-gwas

Development Installation
~~~~~~~~~~~~~~~~~~~~~~~~

For contributors and developers:

.. code-block:: bash

   # Clone the repository
   git clone https://github.com/YOUR-USERNAME/edge-gwas.git
   cd edge-gwas
   
   # Install in editable mode with development dependencies
   pip install -e ".[dev]"

Virtual Environment Setup
-------------------------

We recommend using a virtual environment:

**Using venv:**

.. code-block:: bash

   python -m venv edge-gwas-env
   source edge-gwas-env/bin/activate  # On Windows: edge-gwas-env\Scripts\activate
   pip install edge-gwas

**Using conda:**

.. code-block:: bash

   conda create -n edge-gwas python=3.9
   conda activate edge-gwas
   pip install edge-gwas

Verify Installation
-------------------

To verify your installation:

.. code-block:: python

   import edge_gwas
   print(edge_gwas.__version__)
   # Should print: 0.1.0

Run a simple test:

.. code-block:: python

   from edge_gwas import gwas
   # Test basic functionality
   print("edge-gwas is ready!")

Troubleshooting
---------------

Common Issues
~~~~~~~~~~~~~

**ImportError: No module named 'edge_gwas'**

Make sure you've activated your virtual environment and installed the package.

**Dependency conflicts**

Try creating a fresh virtual environment:

.. code-block:: bash

   python -m venv fresh-env
   source fresh-env/bin/activate
   pip install edge-gwas

**Permission errors**

Use the ``--user`` flag:

.. code-block:: bash

   pip install --user edge-gwas

Platform-Specific Notes
-----------------------

**Windows:**
- Use Command Prompt or PowerShell
- Some scientific packages may require Visual C++ Build Tools

**macOS:**
- May need to install Xcode Command Line Tools: ``xcode-select --install``

**Linux:**
- Usually works out of the box
- On Ubuntu/Debian: ``sudo apt-get install python3-dev``

Next Steps
----------

After installation, check out:

* :ref:`quickstart` - Quick start guide
* :ref:`examples` - Example workflows
* :ref:`api` - API reference
