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

edge-gwas v0.1.1 includes all dependencies for complete functionality:

**Core Analysis:**

* numpy >= 1.19.0
* pandas >= 1.1.0
* scipy >= 1.5.0
* statsmodels >= 0.12.0
* scikit-learn >= 0.24.0
* pandas-plink >= 2.0.0
* joblib >= 1.0.0
* matplotlib >= 3.3.0

**File Format Support:**

* pgenlib >= 0.81.0 (PLINK 2 format: .pgen/.pvar/.psam)
* bgen-reader >= 4.0.8 (BGEN format with dosages)
* cyvcf2 >= 0.30.0 (VCF/VCF.GZ format)

**Development Tools:**

* pytest >= 6.2.0 (testing)
* black >= 21.0 (code formatting)
* sphinx >= 4.0.0 (documentation)

All dependencies are installed automatically!

Installation from GitHub
------------------------

**Recommended method for v0.1.1:**

.. code-block:: bash

   pip install git+https://github.com/nicenzhou/edge-gwas.git

This single command installs edge-gwas with **all file format support** included.

Installation from Source
------------------------

.. code-block:: bash

   # Clone the repository
   git clone https://github.com/nicenzhou/edge-gwas.git
   cd edge-gwas
   
   # Install all dependencies
   pip install -r requirements.txt
   
   # Install the package in development mode
   pip install -e .

**Or install directly:**

.. code-block:: bash

   cd edge-gwas
   pip install .

Virtual Environment Setup
-------------------------

We recommend using a virtual environment to avoid dependency conflicts.

**Using venv:**

.. code-block:: bash

   # Create virtual environment
   python3 -m venv edge-gwas-env
   
   # Activate (Linux/Mac)
   source edge-gwas-env/bin/activate
   
   # Activate (Windows)
   edge-gwas-env\Scripts\activate
   
   # Install edge-gwas
   pip install git+https://github.com/nicenzhou/edge-gwas.git

**Using conda:**

.. code-block:: bash

   # Create conda environment
   conda create -n edge-gwas python=3.9
   
   # Activate environment
   conda activate edge-gwas
   
   # Install edge-gwas
   pip install git+https://github.com/nicenzhou/edge-gwas.git

Verify Installation
-------------------

**In Python/Jupyter Notebook:**

.. code-block:: python

   from edge_gwas import EDGEAnalysis
   from edge_gwas.utils import (
       load_plink_data,
       load_pgen_data,
       load_bgen_data,
       load_vcf_data
   )
   from edge_gwas.visualize import manhattan_plot, qq_plot
   
   print("✓ edge-gwas v0.1.1 installed successfully!")
   print("✓ All file formats supported!")

**Using command line (Mac/Linux):**

.. code-block:: bash

   python3 -c "from edge_gwas import EDGEAnalysis; print('✓ Installed successfully')"

**Check version:**

.. code-block:: python

   import edge_gwas
   print(f"edge-gwas version: {edge_gwas.__version__}")

Supported File Formats
----------------------

edge-gwas v0.1.1 supports all major genotype file formats out of the box:

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - Format
     - Function
     - Description
   * - **PLINK**
     - ``load_plink_data()``
     - Standard PLINK binary (.bed/.bim/.fam)
   * - **PLINK 2**
     - ``load_pgen_data()``
     - New PLINK 2 format (.pgen/.pvar/.psam)
   * - **BGEN**
     - ``load_bgen_data()``
     - UK Biobank format with dosages (.bgen)
   * - **VCF**
     - ``load_vcf_data()``
     - Variant Call Format (.vcf/.vcf.gz)

All file formats work immediately after installation!

Troubleshooting
---------------

Common Issues
~~~~~~~~~~~~~

**ImportError: No module named 'edge_gwas'**

Ensure virtual environment is activated and package is installed:

.. code-block:: bash

   # Check if package is installed
   pip list | grep edge-gwas
   
   # Reinstall if needed
   pip install --force-reinstall git+https://github.com/nicenzhou/edge-gwas.git

**Dependency conflicts**

Create a fresh environment:

.. code-block:: bash

   # Using conda
   conda create -n edge-gwas-fresh python=3.9
   conda activate edge-gwas-fresh
   pip install git+https://github.com/nicenzhou/edge-gwas.git
   
   # Using venv
   python3 -m venv fresh-env
   source fresh-env/bin/activate
   pip install git+https://github.com/nicenzhou/edge-gwas.git

**Permission errors**

Use ``--user`` flag:

.. code-block:: bash

   pip install --user git+https://github.com/nicenzhou/edge-gwas.git

**File format specific errors**

If you encounter issues with specific file formats:

.. code-block:: bash

   # PLINK 2 format issues
   pip install --upgrade pgenlib
   
   # BGEN format issues
   pip install --upgrade bgen-reader
   
   # VCF format issues
   pip install --upgrade cyvcf2

Platform-Specific Notes
-----------------------

**Windows:**

- Use Command Prompt or PowerShell
- May require `Microsoft C++ Build Tools <https://visualstudio.microsoft.com/visual-cpp-build-tools/>`_
- For cyvcf2, you may need to use WSL (Windows Subsystem for Linux)

**macOS:**

- May need Xcode Command Line Tools:

.. code-block:: bash

   xcode-select --install

- For Apple Silicon (M1/M2), use:

.. code-block:: bash

   # Create conda environment with x86_64 architecture
   CONDA_SUBDIR=osx-64 conda create -n edge-gwas python=3.9
   conda activate edge-gwas
   conda config --env --set subdir osx-64
   pip install git+https://github.com/nicenzhou/edge-gwas.git

**Linux:**

- Ubuntu/Debian:

.. code-block:: bash

   sudo apt-get update
   sudo apt-get install python3-dev python3-pip build-essential

- CentOS/RHEL:

.. code-block:: bash

   sudo yum install python3-devel gcc gcc-c++

Upgrading from v0.1.0
---------------------

If you have v0.1.0 installed, upgrade to v0.1.1:

.. code-block:: bash

   pip install --upgrade git+https://github.com/nicenzhou/edge-gwas.git

**What's New in v0.1.1:**

* PLINK 2 format support (.pgen/.pvar/.psam)
* BGEN format support with dosages
* VCF/VCF.GZ format support
* Hardy-Weinberg Equilibrium filtering
* Sample call rate filtering
* Case/control balance checking
* Standard additive GWAS for comparison
* Cross-validation support
* Enhanced quality control tools

Development Installation
------------------------

For contributors and developers:

.. code-block:: bash

   # Clone repository
   git clone https://github.com/nicenzhou/edge-gwas.git
   cd edge-gwas
   
   # Create development environment
   conda create -n edge-gwas-dev python=3.9
   conda activate edge-gwas-dev
   
   # Install all dependencies (includes dev tools)
   pip install -r requirements.txt
   
   # Install in editable mode
   pip install -e .
   
   # Run tests to verify
   pytest tests/

See :ref:`contributing` for detailed development setup.

Docker Installation (Optional)
-------------------------------

For isolated environments:

.. code-block:: bash

   # Create Dockerfile
   FROM python:3.9-slim
   
   RUN apt-get update && apt-get install -y git build-essential
   RUN pip install git+https://github.com/nicenzhou/edge-gwas.git
   
   WORKDIR /workspace
   CMD ["python"]

.. code-block:: bash

   # Build and run
   docker build -t edge-gwas .
   docker run -it -v $(pwd):/workspace edge-gwas

Getting Help
------------

If you encounter installation issues:

1. Check the `GitHub Issues <https://github.com/nicenzhou/edge-gwas/issues>`_
2. Search `GitHub Discussions <https://github.com/nicenzhou/edge-gwas/discussions>`_
3. Email: jyzhou@stanford.edu (code questions)
4. Email: molly.hall@pennmedicine.upenn.edu (research questions)

Next Steps
----------

After successful installation:

* :ref:`quickstart` - Quick start guide with minimal example
* :ref:`examples` - Complete example workflows
* :ref:`user_guide` - Detailed usage guide
* :ref:`api_reference` - Full API documentation

**Ready to start?** Try the :ref:`quickstart` tutorial!

---

**Version:** 0.1.1

**Last Updated:** January 2025
