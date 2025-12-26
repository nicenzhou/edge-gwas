.. _installation:

Installation Guide
==================

System Requirements
-------------------

* Python 3.7 or higher
* pip package manager
* 4GB RAM minimum (8GB+ recommended for large datasets)
* 16GB+ RAM recommended for biobank-scale data with GRM

Dependencies
------------

edge-gwas v0.1.1 requires the following dependencies:

**Core Analysis:**

* numpy >= 1.19.0
* pandas >= 1.2.0 (v0.1.1 uses pandas, not Koalas)
* scipy >= 1.6.0
* statsmodels >= 0.12.0
* scikit-learn >= 0.24.0
* pandas-plink >= 2.0.0
* joblib >= 1.0.0
* matplotlib >= 3.3.0

**Optional File Format Support:**

* pgenlib >= 0.81.0 (for PLINK 2 format: .pgen/.pvar/.psam)
* bgen-reader >= 4.0.8 (for BGEN format with dosages)
* cyvcf2 >= 0.30.0 (for VCF/VCF.GZ format)

**Optional External Tools (NEW in v0.1.1):**

* PLINK2 (for PCA calculation)
* GCTA (for GRM calculation)
* R with Bioconductor packages (for PC-AiR):
  
  * GENESIS
  * SNPRelate
  * gdsfmt

Core dependencies are installed automatically. Optional dependencies can be installed as needed.

Quick Installation
------------------

**Install edge-gwas with core dependencies:**

.. code-block:: bash

   pip install git+https://github.com/nicenzhou/edge-gwas.git

**Install with all optional file format support:**

.. code-block:: bash

   pip install git+https://github.com/nicenzhou/edge-gwas.git
   pip install pgenlib bgen-reader cyvcf2

**Install external tools for PCA and GRM (NEW in v0.1.1):**

After installing edge-gwas, run the interactive installer:

.. code-block:: bash

   edge-gwas-install-tools

This will install:

* PLINK2 (for PCA)
* GCTA (for GRM calculation)
* R packages (GENESIS, SNPRelate, gdsfmt for PC-AiR)

Installation from GitHub
------------------------

**Standard installation:**

.. code-block:: bash

   pip install git+https://github.com/nicenzhou/edge-gwas.git

**With specific version:**

.. code-block:: bash

   pip install git+https://github.com/nicenzhou/edge-gwas.git@v0.1.1

Installation from Source
------------------------

.. code-block:: bash

   # Clone the repository
   git clone https://github.com/nicenzhou/edge-gwas.git
   cd edge-gwas
   
   # Install core dependencies
   pip install -r requirements.txt
   
   # Install the package in development mode
   pip install -e .
   
   # (Optional) Install file format support
   pip install pgenlib bgen-reader cyvcf2
   
   # (Optional) Install external tools
   edge-gwas-install-tools

Virtual Environment Setup
-------------------------

We **strongly recommend** using a virtual environment to avoid dependency conflicts.

**Using venv (Python built-in):**

.. code-block:: bash

   # Create virtual environment
   python3 -m venv edge-gwas-env
   
   # Activate (Linux/Mac)
   source edge-gwas-env/bin/activate
   
   # Activate (Windows)
   edge-gwas-env\Scripts\activate
   
   # Install edge-gwas
   pip install git+https://github.com/nicenzhou/edge-gwas.git
   
   # Install optional dependencies
   pip install pgenlib bgen-reader cyvcf2
   
   # Install external tools
   edge-gwas-install-tools

**Using conda (recommended for complex dependencies):**

.. code-block:: bash

   # Create conda environment
   conda create -n edge-gwas python=3.9
   
   # Activate environment
   conda activate edge-gwas
   
   # Install edge-gwas
   pip install git+https://github.com/nicenzhou/edge-gwas.git
   
   # Install optional dependencies
   pip install pgenlib bgen-reader cyvcf2
   
   # Install external tools
   edge-gwas-install-tools

External Tools Installation (NEW in v0.1.1)
--------------------------------------------

edge-gwas v0.1.1 includes an interactive installer for external tools.

**Automatic Installation:**

.. code-block:: bash

   edge-gwas-install-tools

This interactive installer will:

1. Detect your operating system and architecture
2. Download and install PLINK2 for your system
3. Download and install GCTA for your system
4. Install R packages (GENESIS, SNPRelate, gdsfmt) if R is available
5. Configure PATH if needed

**Supported Platforms:**

* Linux (x86_64)
* macOS (Intel x86_64)
* macOS (Apple Silicon ARM64)

**What Gets Installed:**

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - Tool
     - Purpose
     - Installation Location
   * - PLINK2
     - PCA calculation
     - ``~/.local/bin/plink2``
   * - GCTA
     - GRM calculation
     - ``~/.local/bin/gcta64``
   * - R packages
     - PC-AiR analysis
     - R library path

**Manual Installation (if automatic fails):**

*PLINK2:*

.. code-block:: bash

   # Download from https://www.cog-genomics.org/plink/2.0/
   # Linux
   wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20251205.zip
   unzip plink2_linux_x86_64_20251205.zip
   chmod +x plink2
   mv plink2 ~/.local/bin/
   
   # macOS (see website for latest)
   # https://www.cog-genomics.org/plink/2.0/

*GCTA:*

.. code-block:: bash

   # Linux
   wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.95.0-linux-kernel-3-x86_64.zip
   unzip gcta-1.95.0-linux-kernel-3-x86_64.zip
   chmod +x gcta-1.95.0-linux-kernel-3-x86_64/gcta
   mv gcta-1.95.0-linux-kernel-3-x86_64/gcta ~/.local/bin/gcta64
   
   # macOS (see website for latest)
   # https://yanglab.westlake.edu.cn/software/gcta/

*R packages:*

.. code-block:: r

   # In R console
   if (!requireNamespace("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   
   BiocManager::install(c("GENESIS", "SNPRelate", "gdsfmt"))

**Add to PATH:**

If tools aren't found after installation, add to your PATH:

.. code-block:: bash

   # Add to ~/.bashrc (Linux) or ~/.zshrc (macOS)
   export PATH="$$HOME/.local/bin:$$PATH"
   
   # Reload configuration
   source ~/.bashrc  # or source ~/.zshrc

Verify Installation
-------------------

**Check core installation:**

.. code-block:: python

   import edge_gwas
   from edge_gwas import EDGEAnalysis
   from edge_gwas.utils import (
       load_plink_data,
       calculate_pca_plink,
       calculate_grm_gcta,
       load_grm_gcta
   )
   
   print(f"edge-gwas version: {edge_gwas.__version__}")
   print("✓ Core package installed successfully!")

**Check file format support:**

.. code-block:: python

   from edge_gwas.utils import (
       load_pgen_data,  # Requires pgenlib
       load_bgen_data,  # Requires bgen-reader
       load_vcf_data    # Requires cyvcf2
   )
   print("✓ All file formats supported!")

**Check external tools (NEW in v0.1.1):**

.. code-block:: bash

   # Verify all tools
   edge-gwas-check-tools
   
   # This will show status of:
   # - Python packages
   # - PLINK2
   # - GCTA
   # - R and R packages

**Expected output:**

.. code-block:: text

   ======================================================================
   EDGE-GWAS External Tools Check
   ======================================================================
   
   Python Packages:
   ----------------------------------------------------------------------
   ✓ numpy: Installed
   ✓ pandas: Installed
   ✓ scipy: Installed
   ✓ statsmodels: Installed
   ✓ sklearn: Installed
   ✓ matplotlib: Installed
   ✓ pandas_plink: Installed
   
   External Tools:
   ----------------------------------------------------------------------
   ✓ PLINK2: PLINK v2.00a3.7 64-bit Intel (24 Jan 2024)
   ✓ GCTA: GCTA 1.95.0
   
   R and Packages:
   ----------------------------------------------------------------------
   ✓ R: R version 4.3.0
   ✓ R package GENESIS: Installed
   ✓ R package SNPRelate: Installed
   ✓ R package gdsfmt: Installed
   
   ======================================================================
   ✓ All tools and packages are properly installed!
   ======================================================================

Supported File Formats
----------------------

edge-gwas v0.1.1 supports all major genotype file formats:

.. list-table::
   :header-rows: 1
   :widths: 15 25 30 30

   * - Format
     - Function
     - Dependency
     - Description
   * - **PLINK**
     - ``load_plink_data()``
     - pandas-plink
     - Standard PLINK binary (.bed/.bim/.fam)
   * - **PLINK 2**
     - ``load_pgen_data()``
     - pgenlib
     - New PLINK 2 format (.pgen/.pvar/.psam)
   * - **BGEN**
     - ``load_bgen_data()``
     - bgen-reader
     - UK Biobank format with dosages
   * - **VCF**
     - ``load_vcf_data()``
     - cyvcf2
     - Variant Call Format (.vcf/.vcf.gz)

Upgrading from v0.1.0
---------------------

If you have v0.1.0 installed, upgrade to v0.1.1:

.. code-block:: bash

   pip install --upgrade --force-reinstall git+https://github.com/nicenzhou/edge-gwas.git

**Breaking Changes in v0.1.1:**

* **Koalas removed**: All code now uses pandas instead of Koalas
* **PCA functions**: Now return DataFrames with IID as index
* **Sample ID handling**: All IDs are converted to strings for consistency

**Migration steps:**

1. Remove any Koalas imports:

   .. code-block:: python
   
      # OLD (v0.1.0)
      import databricks.koalas as ks
      df = data.to_koalas()
      
      # NEW (v0.1.1)
      import pandas as pd
      df = data  # Already pandas

2. Update PCA usage:

   .. code-block:: python
   
      # v0.1.1: PCA DataFrames have IID as index
      pca_df = calculate_pca_plink('genotypes', n_pcs=10)
      # pca_df.index contains sample IDs

3. Install external tools:

   .. code-block:: bash
   
      edge-gwas-install-tools

**New Features in v0.1.1:**

* GRM support for population structure control
* Multiple PCA methods (PLINK2, PC-AiR, sklearn)
* Outcome transformations for continuous traits
* Additional file format support (PGEN, BGEN, VCF)
* Enhanced quality control functions
* Cross-validation support
* Additive GWAS for comparison

Troubleshooting
---------------

Common Issues
~~~~~~~~~~~~~

**1. ImportError: No module named 'edge_gwas'**

.. code-block:: bash

   # Check if installed
   pip list | grep edge-gwas
   
   # Reinstall
   pip install --force-reinstall git+https://github.com/nicenzhou/edge-gwas.git

**2. ImportError: No module named 'databricks.koalas' (from v0.1.0)**

.. code-block:: bash

   # v0.1.1 no longer uses Koalas
   # Upgrade to v0.1.1
   pip install --upgrade --force-reinstall git+https://github.com/nicenzhou/edge-gwas.git

**3. PLINK2/GCTA not found**

.. code-block:: bash

   # Run installer
   edge-gwas-install-tools
   
   # Or check PATH
   echo $PATH
   
   # Add to PATH if needed
   export PATH="$$HOME/.local/bin:$$PATH"

**4. R package installation fails**

.. code-block:: r

   # Install manually in R
   if (!requireNamespace("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   
   BiocManager::install(c("GENESIS", "SNPRelate", "gdsfmt"))

**5. File format specific errors**

.. code-block:: bash

   # PLINK 2 format
   pip install --upgrade pgenlib
   
   # BGEN format
   pip install --upgrade bgen-reader
   
   # VCF format
   pip install --upgrade cyvcf2

**6. Memory errors with GRM**

.. code-block:: bash

   # Use sparse GRM for large datasets
   from edge_gwas.utils import calculate_grm_gcta
   
   grm_prefix = calculate_grm_gcta(
       'genotypes',
       method='grm-sparse'  # Use sparse instead of full
   )

**7. Permission errors**

.. code-block:: bash

   # Install to user directory
   pip install --user git+https://github.com/nicenzhou/edge-gwas.git

Platform-Specific Notes
-----------------------

**Windows:**

- External tool installer (``edge-gwas-install-tools``) not supported on Windows
- Install PLINK2 and GCTA manually from their websites
- For cyvcf2, consider using WSL (Windows Subsystem for Linux)
- Core Python package works on Windows

.. code-block:: bash

   # Windows installation (core only)
   pip install git+https://github.com/nicenzhou/edge-gwas.git

**macOS:**

- May need Xcode Command Line Tools:

.. code-block:: bash

   xcode-select --install

- Apple Silicon (M1/M2/M3):

.. code-block:: bash

   # Tool installer automatically detects ARM64
   edge-gwas-install-tools
   
   # If issues, try Rosetta 2:
   arch -x86_64 pip install git+https://github.com/nicenzhou/edge-gwas.git

**Linux:**

- Ubuntu/Debian:

.. code-block:: bash

   sudo apt-get update
   sudo apt-get install python3-dev python3-pip build-essential

- CentOS/RHEL:

.. code-block:: bash

   sudo yum install python3-devel gcc gcc-c++

- For HPC clusters, use conda:

.. code-block:: bash

   module load anaconda3  # If using module system
   conda create -n edge-gwas python=3.9
   conda activate edge-gwas
   pip install git+https://github.com/nicenzhou/edge-gwas.git

**UK Biobank Research Analysis Platform (RAP):**

.. code-block:: bash

   # On RAP, install in user space
   pip install --user git+https://github.com/nicenzhou/edge-gwas.git
   
   # External tools may already be available
   which plink2
   which gcta64
   
   # If not, install to user bin
   mkdir -p ~/.local/bin
   edge-gwas-install-tools

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
   pip install -r requirements-dev.txt
   
   # Install in editable mode
   pip install -e .
   
   # Install external tools
   edge-gwas-install-tools
   
   # Run tests to verify
   pytest tests/
   
   # Check code style
   black --check edge_gwas/
   
   # Build documentation
   cd docs
   make html

See :ref:`contributing` for detailed development guidelines.

Docker Installation (Optional)
-------------------------------

For containerized deployments:

**Dockerfile:**

.. code-block:: dockerfile

   FROM python:3.9-slim
   
   # Install system dependencies
   RUN apt-get update && apt-get install -y \
       git \
       build-essential \
       curl \
       wget \
       && rm -rf /var/lib/apt/lists/*
   
   # Install edge-gwas
   RUN pip install --no-cache-dir \
       git+https://github.com/nicenzhou/edge-gwas.git
   
   # Install optional file format support
   RUN pip install --no-cache-dir \
       pgenlib \
       bgen-reader \
       cyvcf2
   
   # Install external tools
   RUN edge-gwas-install-tools
   
   # Set working directory
   WORKDIR /workspace
   
   # Default command
   CMD ["python"]

**Build and run:**

.. code-block:: bash

   # Build image
   docker build -t edge-gwas:0.1.1 .
   
   # Run interactively
   docker run -it -v $(pwd):/workspace edge-gwas:0.1.1
   
   # Run specific script
   docker run -v $(pwd):/workspace edge-gwas:0.1.1 python my_analysis.py

**Docker Compose (for complex workflows):**

.. code-block:: yaml

   # docker-compose.yml
   version: '3.8'
   
   services:
     edge-gwas:
       build: .
       volumes:
         - ./data:/workspace/data
         - ./output:/workspace/output
       environment:
         - PYTHONUNBUFFERED=1
       command: python analysis.py

.. code-block:: bash

   docker-compose up

Singularity/Apptainer (for HPC)
--------------------------------

For HPC environments that use Singularity/Apptainer:

**Definition file (edge-gwas.def):**

.. code-block:: singularity

   Bootstrap: docker
   From: python:3.9-slim
   
   %post
       apt-get update
       apt-get install -y git build-essential curl wget
       pip install git+https://github.com/nicenzhou/edge-gwas.git
       pip install pgenlib bgen-reader cyvcf2
       edge-gwas-install-tools
       apt-get clean
   
   %environment
       export PATH=/root/.local/bin:$PATH
   
   %runscript
       exec python "$@"

**Build and run:**

.. code-block:: bash

   # Build container
   singularity build edge-gwas.sif edge-gwas.def
   
   # Run Python script
   singularity exec edge-gwas.sif python my_analysis.py
   
   # Interactive shell
   singularity shell edge-gwas.sif

Testing Your Installation
--------------------------

**Quick test:**

.. code-block:: python

   #!/usr/bin/env python3
   """Quick installation test for edge-gwas v0.1.1"""
   
   import sys
   
   def test_installation():
       """Test edge-gwas installation."""
       
       print("Testing edge-gwas v0.1.1 installation...")
       print("-" * 70)
       
       # Test core import
       try:
           import edge_gwas
           print(f"✓ edge-gwas version: {edge_gwas.__version__}")
       except ImportError as e:
           print(f"✗ Failed to import edge-gwas: {e}")
           return False
       
       # Test core functionality
       try:
           from edge_gwas import EDGEAnalysis
           edge = EDGEAnalysis(outcome_type='binary')
           print("✓ EDGEAnalysis class available")
       except Exception as e:
           print(f"✗ Failed to create EDGEAnalysis: {e}")
           return False
       
       # Test utilities
       try:
           from edge_gwas.utils import (
               load_plink_data,
               calculate_pca_plink,
               calculate_grm_gcta,
               load_grm_gcta,
               filter_variants_by_maf
           )
           print("✓ Core utility functions available")
       except ImportError as e:
           print(f"✗ Failed to import utilities: {e}")
           return False
       
       # Test optional file formats
       optional_formats = {
           'PGEN': 'load_pgen_data',
           'BGEN': 'load_bgen_data',
           'VCF': 'load_vcf_data'
       }
       
       for fmt, func in optional_formats.items():
           try:
               exec(f"from edge_gwas.utils import {func}")
               print(f"✓ {fmt} format support available")
           except ImportError:
               print(f"⚠ {fmt} format support not available (install optional dependency)")
       
       # Test visualization
       try:
           from edge_gwas.visualize import manhattan_plot, qq_plot
           print("✓ Visualization functions available")
       except ImportError as e:
           print(f"✗ Failed to import visualization: {e}")
           return False
       
       print("-" * 70)
       print("✓ Core installation test passed!")
       print("\nRun 'edge-gwas-check-tools' to verify external tools.")
       return True
   
   if __name__ == '__main__':
       success = test_installation()
       sys.exit(0 if success else 1)

Save as ``test_install.py`` and run:

.. code-block:: bash

   python test_install.py

**Full system check:**

.. code-block:: bash

   # Check all components
   edge-gwas-check-tools
   
   # This verifies:
   # - Python packages
   # - External tools (PLINK2, GCTA)
   # - R and R packages

Installation Checklist
----------------------

Use this checklist to ensure complete installation:

.. list-table::
   :header-rows: 1
   :widths: 50 25 25

   * - Component
     - Required
     - Installed?
   * - Python 3.7+
     - Yes
     - ☐
   * - edge-gwas package
     - Yes
     - ☐
   * - Core dependencies (numpy, pandas, etc.)
     - Yes
     - ☐
   * - PLINK binary format support
     - Yes
     - ☐
   * - PLINK2 format support (pgenlib)
     - Optional
     - ☐
   * - BGEN format support (bgen-reader)
     - Optional
     - ☐
   * - VCF format support (cyvcf2)
     - Optional
     - ☐
   * - PLINK2 (for PCA)
     - Optional
     - ☐
   * - GCTA (for GRM)
     - Optional
     - ☐
   * - R + GENESIS (for PC-AiR)
     - Optional
     - ☐

Minimum Installation (Core Only)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For basic EDGE analysis without external tools:

.. code-block:: bash

   pip install git+https://github.com/nicenzhou/edge-gwas.git

**Features available:**

- Basic EDGE analysis (calculate_alpha, apply_alpha)
- PLINK binary format support
- Standard quality control
- Visualization (Manhattan, QQ plots)

**Features NOT available:**

- PCA calculation (use pre-calculated PCs)
- GRM calculation (use pre-calculated GRM)
- PC-AiR (use standard PCA)
- PGEN/BGEN/VCF formats

Recommended Installation (Full Features)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For complete functionality including population structure control:

.. code-block:: bash

   # Install package with optional file formats
   pip install git+https://github.com/nicenzhou/edge-gwas.git
   pip install pgenlib bgen-reader cyvcf2
   
   # Install external tools
   edge-gwas-install-tools

**All features available:**

- Complete EDGE analysis pipeline
- All file format support
- PCA calculation (PLINK2, PC-AiR, sklearn)
- GRM calculation and usage
- Advanced quality control
- Cross-validation
- Comprehensive visualization

Configuration
-------------

**Environment variables (optional):**

.. code-block:: bash

   # Add to ~/.bashrc or ~/.zshrc
   
   # Path to external tools (if not using default location)
   export PLINK2_PATH="$HOME/.local/bin/plink2"
   export GCTA_PATH="$HOME/.local/bin/gcta64"
   
   # Temporary directory for large datasets
   export TMPDIR="/scratch/$USER/tmp"
   
   # Number of threads for parallel operations
   export OMP_NUM_THREADS=8

**Python configuration:**

.. code-block:: python

   # In your analysis script
   import os
   
   # Set number of threads for numpy/scipy
   os.environ['OMP_NUM_THREADS'] = '8'
   os.environ['OPENBLAS_NUM_THREADS'] = '8'
   os.environ['MKL_NUM_THREADS'] = '8'
   
   import edge_gwas

Uninstallation
--------------

**Remove edge-gwas:**

.. code-block:: bash

   pip uninstall edge-gwas

**Remove external tools:**

.. code-block:: bash

   # Remove binaries
   rm ~/.local/bin/plink2
   rm ~/.local/bin/gcta64
   
   # Remove R packages (in R console)
   remove.packages(c("GENESIS", "SNPRelate", "gdsfmt"))

**Remove conda environment:**

.. code-block:: bash

   conda deactivate
   conda env remove -n edge-gwas

Getting Help
------------

If you encounter installation issues:

**1. Check the documentation:**

* :ref:`troubleshooting` - Common issues and solutions
* :ref:`faq` - Frequently asked questions

**2. Search existing issues:**

* `GitHub Issues <https://github.com/nicenzhou/edge-gwas/issues>`_
* `GitHub Discussions <https://github.com/nicenzhou/edge-gwas/discussions>`_

**3. Get support:**

* Report bugs: https://github.com/nicenzhou/edge-gwas/issues/new
* Ask questions: https://github.com/nicenzhou/edge-gwas/discussions/new
* Email support:
  
  * Code questions: jyzhou@stanford.edu
  * Research questions: molly.hall@pennmedicine.upenn.edu

**4. Provide details when reporting issues:**

.. code-block:: bash

   # Include this information in bug reports:
   python --version
   pip list | grep edge-gwas
   edge-gwas-check-tools  # Output from this command
   
   # Operating system
   uname -a  # Linux/Mac
   
   # Error message and full traceback

Next Steps
----------

After successful installation:

1. **Verify installation**: Run ``edge-gwas-check-tools``
2. **Quick start**: Try the :ref:`quickstart` tutorial
3. **Learn the basics**: Read the :ref:`user_guide`
4. **Explore examples**: See :ref:`examples` for complete workflows
5. **API reference**: Check :ref:`api_reference` for detailed documentation

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
