#!/bin/bash
# Installation script for EDGE-GWAS with all dependencies

set -e

echo "================================================================================"
echo "EDGE-GWAS Installation Script"
echo "================================================================================"
echo ""

# Check Python version
echo "Checking Python version..."
python_version=$(python3 --version 2>&1 | awk '{print $2}')
required_version="3.7"

if [ "$(printf '%s\n' "$required_version" "$python_version" | sort -V | head -n1)" != "$required_version" ]; then 
    echo "Error: Python 3.7 or higher is required. Found: $python_version"
    exit 1
fi
echo "✓ Python $python_version found"
echo ""

# Create virtual environment
echo "Creating virtual environment..."
python3 -m venv edge-gwas-env
source edge-gwas-env/bin/activate
echo "✓ Virtual environment created and activated"
echo ""

# Upgrade pip
echo "Upgrading pip..."
pip install --upgrade pip setuptools wheel
echo "✓ pip upgraded"
echo ""

# Install package
echo "Installing EDGE-GWAS..."
pip install -e .
echo "✓ EDGE-GWAS installed"
echo ""

# Install optional dependencies
echo "Would you like to install optional dependencies? (y/N)"
read -r response
if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]; then
    echo "Installing optional dependencies..."
    
    # Install PGEN support
    pip install pgenlib || echo "Warning: pgenlib installation failed"
    
    # Install BGEN support
    pip install bgen-reader || echo "Warning: bgen-reader installation failed"
    
    # Install VCF support
    pip install cyvcf2 || echo "Warning: cyvcf2 installation failed"
    
    # Install visualization enhancements
    pip install seaborn plotly || echo "Warning: visualization packages installation failed"
    
    echo "✓ Optional dependencies installed"
fi
echo ""

# Install external tools
echo "Installing external tools (PLINK2, GCTA, R packages)..."
python -c "from setup import ExternalToolsInstaller; ExternalToolsInstaller().install_all()"
echo ""

# Check installation
echo "Verifying installation..."
edge-gwas-check-tools
echo ""

echo "================================================================================"
echo "Installation complete!"
echo "================================================================================"
echo ""
echo "To activate the environment in the future, run:"
echo "  source edge-gwas-env/bin/activate"
echo ""
echo "To verify installation, run:"
echo "  edge-gwas-check-tools"
echo ""
echo "To get started, see the documentation at:"
echo "  https://edge-gwas.readthedocs.io/"
echo ""
echo "================================================================================"
