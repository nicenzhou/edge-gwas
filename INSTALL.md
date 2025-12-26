# Installation Guide for EDGE-GWAS

This guide provides detailed installation instructions for EDGE-GWAS and all its dependencies.

## Quick Install (Recommended)

The easiest way to install EDGE-GWAS with all dependencies:

```bash
# Clone the repository
git clone https://github.com/nicenzhou/edge-gwas.git
cd edge-gwas

# Run the installation script
chmod +x install.sh
./install.sh
```

This will:
1. Create a virtual environment
2. Install Python dependencies
3. Install PLINK2, GCTA, and R packages automatically
4. Verify the installation

## Manual Installation

### Step 1: Python Package

```bash
# Create virtual environment
python3 -m venv edge-gwas-env
source edge-gwas-env/bin/activate

# Install from PyPI (when available)
pip install edge-gwas

# Or install from source
git clone https://github.com/nicenzhou/edge-gwas.git
cd edge-gwas
pip install -e .
```

### Step 2: External Tools

The external tools will be installed automatically during package installation.

If automatic installation fails, you can:

1. **Install manually** following the instructions below
2. **Run the interactive installer**:
   ```bash
   edge-gwas-install-tools
   ```

#### Manual Installation of External Tools

**PLINK2:**
```bash
# Linux
wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20251205.zip
unzip plink2_linux_x86_64_20251205.zip
sudo mv plink2 /usr/local/bin/
```

**GCTA:**
```bash
# Linux
wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.95.0-linux-kernel-3-x86_64.zip
unzip gcta-1.95.0-linux-kernel-3-x86_64.zip
sudo mv gcta-1.95.0-linux-kernel-3-x86_64/gcta64 /usr/local/bin/
```

**R Packages:**
```R
# In R console
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("GENESIS", "SNPRelate", "gdsfmt"))
```

### Step 3: Verify Installation

```bash
edge-gwas-check-tools
```

## Installation Options

### Minimal Installation (Core Features Only)

```bash
pip install edge-gwas
```

This installs only Python dependencies. You can still use:
- Basic EDGE analysis
- sklearn-based PCA
- Data loading and processing

### Full Installation (All Features)

```bash
# Install with all optional dependencies
pip install edge-gwas[all]

# Or install optional packages individually
pip install pgenlib  # PGEN format support
pip install bgen-reader  # BGEN format support
pip install cyvcf2  # VCF format support
pip install seaborn plotly  # Enhanced visualizations
```

### Development Installation

For contributing to EDGE-GWAS:

```bash
git clone https://github.com/nicenzhou/edge-gwas.git
cd edge-gwas

# Create development environment
python3 -m venv dev-env
source dev-env/bin/activate

# Install in development mode with dev dependencies
pip install -e ".[dev]"

# Install pre-commit hooks
pre-commit install
```

## Platform-Specific Instructions

### Linux (Ubuntu/Debian)

```bash
# Install system dependencies
sudo apt-get update
sudo apt-get install -y python3-dev python3-pip python3-venv
sudo apt-get install -y build-essential libcurl4-openssl-dev libssl-dev libxml2-dev
sudo apt-get install -y r-base r-base-dev

# Install EDGE-GWAS
pip install edge-gwas

# Install external tools
edge-gwas-install-tools
```

### macOS

```bash
# Install Homebrew if not already installed
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install dependencies
brew install python3 r

# Install EDGE-GWAS
pip3 install edge-gwas

# Install external tools
edge-gwas-install-tools
```

### Windows (via WSL)

EDGE-GWAS is best used on Linux/macOS. For Windows users, we recommend Windows Subsystem for Linux (WSL):

```bash
# In WSL Ubuntu
sudo apt-get update
sudo apt-get install -y python3 python3-pip python3-venv r-base

# Install EDGE-GWAS
pip3 install edge-gwas

# Install external tools
edge-gwas-install-tools
```

## Docker Installation

For a containerized environment:

```dockerfile
# Dockerfile
FROM python:3.9-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    r-base \
    r-base-dev \
    wget \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# Install EDGE-GWAS
RUN pip install edge-gwas

# Install external tools
RUN edge-gwas-install-tools

WORKDIR /data
CMD ["/bin/bash"]
```

Build and run:

```bash
# Build the Docker image
docker build -t edge-gwas .

# Run container
docker run -it --rm -v $(pwd):/data edge-gwas
```

Or use pre-built image (when available):

```bash
# Pull the Docker image
docker pull nicenzhou/edge-gwas:latest

# Run container
docker run -it --rm -v $(pwd):/data nicenzhou/edge-gwas:latest
```

## Conda Installation

```bash
# Create conda environment
conda create -n edge-gwas python=3.9
conda activate edge-gwas

# Install from conda-forge (when available)
conda install -c conda-forge edge-gwas

# Or install with pip in conda environment
pip install edge-gwas

# Install external tools
edge-gwas-install-tools
```

## Installation on HPC/Cluster

For high-performance computing environments:

```bash
# Load required modules (example for SLURM-based cluster)
module load python/3.9
module load r/4.2

# Install in user directory
pip install --user edge-gwas

# Install external tools to user directory
export PATH="$$HOME/.local/bin:$$PATH"
edge-gwas-install-tools

# Add to PATH in ~/.bashrc
echo 'export PATH="$$HOME/.local/bin:$$PATH"' >> ~/.bashrc
source ~/.bashrc
```

### HPC-Specific Configuration

```bash
# Create a job script (example: install_edge.sh)
#!/bin/bash
#SBATCH --job-name=edge_install
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

module load python/3.9
module load r/4.2

pip install --user edge-gwas
edge-gwas-install-tools

# Submit the job
sbatch install_edge.sh
```

## Troubleshooting

### Common Issues

**1. Permission denied errors:**
```bash
# Install to user directory instead of system-wide
pip install --user edge-gwas

# Or use virtual environment
python3 -m venv edge-env
source edge-env/bin/activate
pip install edge-gwas
```

**2. PLINK2/GCTA not found after installation:**
```bash
# Add to PATH manually
export PATH="$$HOME/.local/bin:$$PATH"

# Make permanent by adding to ~/.bashrc or ~/.zshrc
echo 'export PATH="$$HOME/.local/bin:$$PATH"' >> ~/.bashrc
source ~/.bashrc

# Or create symlinks
ln -s ~/.local/bin/plink2 /usr/local/bin/plink2
ln -s ~/.local/bin/gcta64 /usr/local/bin/gcta64
```

**3. R package installation fails:**
```bash
# Install system dependencies first
# Ubuntu/Debian:
sudo apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev libbz2-dev liblzma-dev

# macOS:
brew install curl openssl libxml2 xz

# Then retry R package installation
Rscript -e "BiocManager::install(c('GENESIS', 'SNPRelate', 'gdsfmt'))"
```

**4. Compilation errors for Python packages:**
```bash
# Install build tools
# Ubuntu/Debian:
sudo apt-get install -y build-essential python3-dev

# macOS:
xcode-select --install

# CentOS/RHEL:
sudo yum install gcc gcc-c++ python3-devel
```

**5. cyvcf2 installation fails:**
```bash
# Install dependencies
# Ubuntu/Debian:
sudo apt-get install -y zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev

# macOS:
brew install zlib bzip2 xz curl

# Then retry
pip install cyvcf2
```

**6. numpy/pandas version conflicts:**
```bash
# Upgrade pip and setuptools
pip install --upgrade pip setuptools wheel

# Install with specific versions
pip install "numpy>=1.19.0,<2.0.0" "pandas>=1.2.0,<3.0.0"

# Then install edge-gwas
pip install edge-gwas
```

**7. SSL certificate errors:**
```bash
# Use trusted host
pip install --trusted-host pypi.org --trusted-host files.pythonhosted.org edge-gwas

# Or update certificates
pip install --upgrade certifi
```

### Verification Steps

After installation, verify everything is working:

```bash
# Check Python package
python -c "import edge_gwas; edge_gwas.print_info()"

# Check external tools
edge-gwas-check-tools

# Run a simple test
python -c "
from edge_gwas import EDGEAnalysis
import numpy as np
import pandas as pd

# Create dummy data
np.random.seed(42)
geno = pd.DataFrame(np.random.randint(0, 3, (100, 10)))
pheno = pd.DataFrame({
    'outcome': np.random.randint(0, 2, 100), 
    'age': np.random.randn(100)
})

# Initialize EDGE
edge = EDGEAnalysis(outcome_type='binary')
print('✓ EDGE-GWAS is working correctly!')
"
```

### Getting Help

If you encounter issues:

1. **Check the documentation**: https://edge-gwas.readthedocs.io/
2. **Search existing issues**: https://github.com/nicenzhou/edge-gwas/issues
3. **Create a new issue**: https://github.com/nicenzhou/edge-gwas/issues/new

Include in your issue:
- Operating system and version (`uname -a` or `cat /etc/os-release`)
- Python version (`python --version`)
- Installation method used
- Complete error message
- Output of `edge-gwas-check-tools`
- Output of `pip list | grep edge`

## Uninstallation

To completely remove EDGE-GWAS:

```bash
# Uninstall Python package
pip uninstall edge-gwas

# Remove external tools (if installed to ~/.local/bin)
rm -f ~/.local/bin/plink2
rm -f ~/.local/bin/gcta64

# Remove R packages (in R console)
remove.packages(c("GENESIS", "SNPRelate", "gdsfmt"))

# Remove virtual environment (if used)
rm -rf edge-gwas-env

# Remove from PATH in ~/.bashrc
# Manually edit ~/.bashrc and remove the line:
# export PATH="$$HOME/.local/bin:$$PATH"
```

## Updating

To update to the latest version:

```bash
# Update from PyPI
pip install --upgrade edge-gwas

# Update from source
cd edge-gwas
git pull
pip install --upgrade -e .

# Update external tools
edge-gwas-install-tools
```

## Version Information

To check your installed version:

```bash
# Command line
python -c "import edge_gwas; print(edge_gwas.__version__)"

# Or check package info
pip show edge-gwas

# In Python
import edge_gwas
print(edge_gwas.__version__)
edge_gwas.print_info()
```

## System Requirements

### Minimum Requirements
- Python 3.7 or higher
- 4 GB RAM
- 1 GB disk space
- Linux, macOS, or Windows (via WSL)

### Recommended Requirements
- Python 3.9 or higher
- 16 GB RAM (for large datasets)
- 10 GB disk space
- Multi-core CPU for parallel processing
- SSD for faster I/O

### For Large-Scale Analyses
- 32+ GB RAM
- 100+ GB disk space
- High-performance computing cluster
- PLINK2 approximate PCA for >5000 samples
- Consider chromosome-by-chromosome analysis

## Next Steps

After successful installation:

1. **Read the Quick Start guide**: See `docs/quickstart.rst` or https://edge-gwas.readthedocs.io/en/latest/quickstart.html
2. **Try the examples**: See `examples/` directory or https://edge-gwas.readthedocs.io/en/latest/examples.html
3. **Review the API documentation**: See `docs/api.rst` or https://edge-gwas.readthedocs.io/en/latest/api.html
4. **Run the test suite**: `pytest tests/` (if you installed from source)

## Support

- **Documentation**: https://edge-gwas.readthedocs.io/
- **GitHub**: https://github.com/nicenzhou/edge-gwas
- **Issues**: https://github.com/nicenzhou/edge-gwas/issues
- **Discussions**: https://github.com/nicenzhou/edge-gwas/discussions
- **Email**: jyzhou@stanford.edu (code questions), molly.hall@pennmedicine.upenn.edu (method questions)

## Citation

If you use EDGE-GWAS in your research, please cite:

```bibtex
@article{zhou2023edgegwas,
  title={Flexibly encoded genome-wide association study identifies novel nonadditive 
         genetic risk variants for cardiometabolic traits},
  author={Zhou, Jiayan and Rico, Andre Luis Garao and Guare, Lindsay and others},
  journal={medRxiv},
  year={2023},
  doi={10.1101/2023.06.01.23290857}
}
```

## Additional Resources

- **PLINK2 documentation**: https://www.cog-genomics.org/plink/2.0/
- **GCTA documentation**: https://yanglab.westlake.edu.cn/software/gcta/
- **GENESIS documentation**: https://bioconductor.org/packages/release/bioc/html/GENESIS.html
- **pandas-plink documentation**: https://pandas-plink.readthedocs.io/

## Installation Tips by Use Case

### For Small Studies (<1000 samples, <100K variants)

```bash
# Minimal installation is sufficient
pip install edge-gwas

# Use sklearn-based PCA
python -c "
from edge_gwas.utils import calculate_pca_sklearn
# Use in your analysis
"
```

### For Medium Studies (1000-10000 samples, 100K-1M variants)

```bash
# Full installation recommended
pip install edge-gwas[all]

# Install external tools
edge-gwas-install-tools

# Use PLINK2 exact PCA
python -c "
from edge_gwas.utils import calculate_pca_plink
pca_df = calculate_pca_plink('data', n_pcs=10, approx=False)
"
```

### For Large Studies (>10000 samples, >1M variants)

```bash
# Full installation required
pip install edge-gwas[all]
edge-gwas-install-tools

# Use approximate PCA and parallel processing
python -c "
from edge_gwas.utils import calculate_pca_plink
from edge_gwas import EDGEAnalysis

# Approximate PCA for speed
pca_df = calculate_pca_plink('data', n_pcs=10, approx=True, approx_samples=5000)

# Use parallel processing
edge = EDGEAnalysis(outcome_type='binary', n_jobs=16)
"
```

### For Studies with Related Individuals

```bash
# Install with R packages
pip install edge-gwas[all]
edge-gwas-install-tools

# Ensure R packages are installed
Rscript -e "library(GENESIS); library(SNPRelate); library(gdsfmt)"

# Use PC-AiR
python -c "
from edge_gwas.utils import calculate_pca_pcair, calculate_grm_gcta

# Calculate GRM
grm_prefix = calculate_grm_gcta('data')

# Run PC-AiR
pca_df = calculate_pca_pcair('data', n_pcs=10, kinship_matrix=grm_prefix)
"
```

### For UK Biobank or Large Biobanks

```bash
# Full installation
pip install edge-gwas[all]

# Install BGEN support
pip install bgen-reader

# Install external tools
edge-gwas-install-tools

# Use approximate PCA
# Process by chromosome
# Use batch processing
```

## Advanced Configuration

### Setting Number of Threads

```bash
# Set environment variables before running
export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8

# Then run your analysis
python your_analysis.py
```

### Custom Installation Directory

```bash
# Install to custom location
pip install --prefix=/custom/path edge-gwas

# Update PYTHONPATH
export PYTHONPATH="/custom/path/lib/python3.9/site-packages:$PYTHONPATH"

# Install external tools to custom location
mkdir -p /custom/path/bin
export PATH="/custom/path/bin:$PATH"

# Run installer with custom prefix
python -c "
from setup import ExternalToolsInstaller
from pathlib import Path
installer = ExternalToolsInstaller()
installer.bin_dir = Path('/custom/path/bin')
installer.install_all()
"
```

### Offline Installation

If you need to install on a system without internet access:

```bash
# On a system with internet:
# Download wheel and dependencies
pip download edge-gwas -d edge-gwas-offline

# Download external tools manually
wget https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20251205.zip
wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.95.0-linux-kernel-3-x86_64.zip

# Transfer the directory to offline system

# On offline system:
pip install --no-index --find-links edge-gwas-offline edge-gwas

# Install external tools manually
unzip plink2_linux_x86_64_latest.zip
mv plink2 ~/.local/bin/

unzip gcta-1.94.1-linux-kernel-3-x86_64.zip
mv gcta-1.94.1-linux-kernel-3-x86_64/gcta64 ~/.local/bin/
```

## Testing Installation

### Quick Test

```bash
# Run basic import test
python -c "import edge_gwas; print('✓ Import successful')"

# Check version
python -c "import edge_gwas; print(f'Version: {edge_gwas.__version__}')"

# Check external tools
edge-gwas-check-tools
```

### Comprehensive Test

```bash
# If installed from source, run test suite
cd edge-gwas
pytest tests/ -v

# Run specific test
pytest tests/test_core.py -v

# Run with coverage
pytest tests/ --cov=edge_gwas --cov-report=html
```

### Integration Test

```python
# test_installation.py
import numpy as np
import pandas as pd
from edge_gwas import EDGEAnalysis
from edge_gwas.utils import (
    calculate_pca_sklearn,
    attach_pcs_to_phenotype,
    stratified_train_test_split
)

def test_basic_workflow():
    """Test basic EDGE workflow."""
    # Create synthetic data
    np.random.seed(42)
    n_samples = 200
    n_variants = 50
    
    # Genotype data
    geno = pd.DataFrame(
        np.random.randint(0, 3, (n_samples, n_variants)),
        columns=[f'SNP{i}' for i in range(n_variants)]
    )
    
    # Phenotype data
    pheno = pd.DataFrame({
        'outcome': np.random.randint(0, 2, n_samples),
        'age': np.random.randn(n_samples),
        'sex': np.random.randint(0, 2, n_samples)
    })
    
    # Calculate PCs
    pca_df = calculate_pca_sklearn(geno, n_pcs=5)
    
    # Attach PCs
    pheno = attach_pcs_to_phenotype(pheno, pca_df, n_pcs=5)
    
    # Split data
    train_g, test_g, train_p, test_p = stratified_train_test_split(
        geno, pheno, 'outcome', test_size=0.5
    )
    
    # Run EDGE
    edge = EDGEAnalysis(outcome_type='binary', n_jobs=2)
    alpha_df, gwas_df = edge.run_full_analysis(
        train_g, train_p, test_g, test_p,
        outcome='outcome',
        covariates=['age', 'sex', 'PC1', 'PC2']
    )
    
    # Check results
    assert len(alpha_df) > 0, "No alpha values calculated"
    assert len(gwas_df) > 0, "No GWAS results"
    assert 'pval' in gwas_df.columns, "Missing p-values"
    
    print("✓ Basic workflow test passed")
    return True

if __name__ == '__main__':
    test_basic_workflow()
    print("\n✓ All installation tests passed!")
```

Run the test:

```bash
python test_installation.py
```

## Environment Management

### Using virtualenv

```bash
# Create environment
python3 -m venv edge-env

# Activate
source edge-env/bin/activate  # Linux/macOS
# or
edge-env\Scripts\activate  # Windows

# Install
pip install edge-gwas

# Deactivate when done
deactivate
```

### Using conda

```bash
# Create environment with specific Python version
conda create -n edge-gwas python=3.9

# Activate
conda activate edge-gwas

# Install
pip install edge-gwas

# Or create from environment file
conda env create -f environment.yml

# Deactivate
conda deactivate
```

### Environment File

Create `environment.yml`:

```yaml
name: edge-gwas
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.9
  - pip
  - numpy>=1.19.0
  - pandas>=1.2.0
  - scipy>=1.6.0
  - scikit-learn>=0.24.0
  - matplotlib>=3.3.0
  - r-base>=4.0
  - pip:
    - edge-gwas
    - pandas-plink
    - statsmodels
```

## Performance Optimization

### Install with Intel MKL

For better performance on Intel CPUs:

```bash
# Using conda (recommended)
conda install -c intel mkl mkl-service

# Or with pip
pip install mkl
```

### Install with OpenBLAS

For AMD CPUs or general use:

```bash
# Using conda
conda install -c conda-forge openblas

# System installation (Ubuntu)
sudo apt-get install libopenblas-dev
```

## Security Considerations

### Verify Package Integrity

```bash
# Check package hash
pip hash edge-gwas

# Verify signatures (when available)
pip install edge-gwas --require-hashes
```

### Install from Trusted Sources

```bash
# Only install from official PyPI
pip install edge-gwas --index-url https://pypi.org/simple

# Or from official GitHub
pip install git+https://github.com/nicenzhou/edge-gwas.git@v0.1.1
```

## Frequently Asked Questions

### Q: Do I need to install external tools?

A: It depends on your use case:
- **Basic EDGE analysis**: No external tools required
- **PCA with PLINK2**: Install PLINK2
- **GRM calculation**: Install GCTA  
- **PC-AiR**: Install R with GENESIS packages

### Q: Can I use EDGE-GWAS without root/sudo access?

A: Yes! Use the `--user` flag or virtual environments:
```bash
pip install --user edge-gwas
```

### Q: What if automatic tool installation fails?

A: Install tools manually following the guide, or use pre-installed versions on your HPC cluster.

### Q: Is EDGE-GWAS compatible with Python 2?

A: No, Python 3.7+ is required.

### Q: Can I install on Windows?

A: Use Windows Subsystem for Linux (WSL) for best compatibility.

### Q: How much disk space do I need?

A: Minimum 1 GB, recommended 10 GB for large analyses.

### Q: Can I run EDGE-GWAS in Jupyter notebooks?

A: Yes! Install Jupyter in the same environment:
```bash
pip install jupyter
jupyter notebook
```

## Troubleshooting FAQ

### Problem: ImportError: No module named 'edge_gwas'

**Solution:**
```bash
# Check if installed
pip list | grep edge-gwas

# Reinstall if necessary
pip install --force-reinstall edge-gwas

# Check Python path
python -c "import sys; print(sys.path)"
```

### Problem: Command 'edge-gwas-check-tools' not found

**Solution:**
```bash
# Add to PATH
export PATH="$HOME/.local/bin:$PATH"

# Or run directly
python -m edge_gwas.check_tools
```

### Problem: R packages won't install

**Solution:**
```bash
# Install system dependencies
sudo apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev

# Install R packages with sudo (if needed)
sudo Rscript -e "BiocManager::install(c('GENESIS', 'SNPRelate', 'gdsfmt'))"
```

### Problem: Out of memory errors

**Solution:**
- Use approximate PCA for large datasets
- Process chromosomes separately
- Increase swap space
- Use a machine with more RAM

### Problem: Slow installation

**Solution:**
```bash
# Use a faster mirror
pip install --index-url https://pypi.tuna.tsinghua.edu.cn/simple edge-gwas

# Or build without testing
pip install --no-deps edge-gwas
pip install -r requirements.txt
```

## Contact for Installation Help

If you're still having issues after following this guide:

- **Email**: jyzhou@stanford.edu
- **GitHub Issues**: https://github.com/nicenzhou/edge-gwas/issues
- **Include**:
  - OS and version
  - Python version
  - Complete error log
  - Output of `pip list`
  - Output of `edge-gwas-check-tools`
