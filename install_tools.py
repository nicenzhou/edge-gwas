"""
Interactive installer for external tools.
"""

import sys
import subprocess
from pathlib import Path


def main():
    """Interactive installation of external tools."""
    print("\n" + "="*70)
    print("EDGE-GWAS External Tools Installer")
    print("="*70 + "\n")
    
    print("This will install the following tools:")
    print("  1. PLINK2 - For PCA calculation")
    print("  2. GCTA - For GRM calculation")
    print("  3. R packages (GENESIS, SNPRelate, gdsfmt) - For PC-AiR")
    print()
    
    response = input("Do you want to proceed? [y/N]: ")
    
    if response.lower() not in ['y', 'yes']:
        print("Installation cancelled.")
        return 1
    
    print("\nStarting installation...\n")
    
    # Import and use the installer
    from setup import ExternalToolsInstaller
    
    installer = ExternalToolsInstaller()
    installer.install_all()
    
    print("\nInstallation complete!")
    print("Run 'edge-gwas-check-tools' to verify installation.")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
