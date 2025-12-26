from setuptools import setup, find_packages
from setuptools.command.install import install
from setuptools.command.develop import develop
import subprocess
import sys
import os
import platform
import urllib.request
import zipfile
import tarfile
import shutil
from pathlib import Path

VERSION = "0.1.1"

class PostInstallCommand(install):
    """Post-installation for installation mode."""
    def run(self):
        install.run(self)
        self.install_external_tools()
    
    def install_external_tools(self):
        """Install external tools after package installation."""
        print("\n" + "="*70)
        print("Installing external tools for EDGE-GWAS...")
        print("="*70 + "\n")
        
        installer = ExternalToolsInstaller()
        installer.install_all()


class PostDevelopCommand(develop):
    """Post-installation for development mode."""
    def run(self):
        develop.run(self)
        installer = ExternalToolsInstaller()
        installer.install_all()


class ExternalToolsInstaller:
    """Handle installation of external tools."""
    
    def __init__(self):
        self.system = platform.system()
        self.home = Path.home()
        self.bin_dir = self.home / '.local' / 'bin'
        self.bin_dir.mkdir(parents=True, exist_ok=True)
        
        # Add to PATH if not already there
        self.update_path()
    
    def update_path(self):
        """Add local bin directory to PATH."""
        path_addition = f'\n# Added by edge-gwas\nexport PATH="$HOME/.local/bin:$PATH"\n'
        
        if self.system in ['Linux', 'Darwin']:
            shell_configs = [
                self.home / '.bashrc',
                self.home / '.zshrc',
                self.home / '.bash_profile'
            ]
            
            for config in shell_configs:
                if config.exists():
                    with open(config, 'r') as f:
                        content = f.read()
                    
                    if 'edge-gwas' not in content:
                        with open(config, 'a') as f:
                            f.write(path_addition)
                        print(f"Updated PATH in {config}")
    
    def install_all(self):
        """Install all external tools."""
        print("\nInstalling external tools...")
        print(f"System: {self.system}")
        print(f"Installation directory: {self.bin_dir}\n")
        
        # Install PLINK2
        try:
            self.install_plink2()
        except Exception as e:
            print(f"Warning: PLINK2 installation failed: {e}")
        
        # Install GCTA
        try:
            self.install_gcta()
        except Exception as e:
            print(f"Warning: GCTA installation failed: {e}")
        
        # Install R packages
        try:
            self.install_r_packages()
        except Exception as e:
            print(f"Warning: R packages installation failed: {e}")
        
        print("\n" + "="*70)
        print("External tools installation complete!")
        print("="*70)
        print("\nIMPORTANT: Please restart your terminal or run:")
        print(f"  source ~/.bashrc  # or source ~/.zshrc")
        print("\nTo verify installation, run:")
        print("  edge-gwas-check-tools")
        print("="*70 + "\n")
    
    def install_plink2(self):
        """Install PLINK2."""
        print("Installing PLINK2...")
        
        if self.system == 'Linux':
            url = 'https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_x86_64_latest.zip'
            archive_name = 'plink2_linux.zip'
            binary_name = 'plink2'
        elif self.system == 'Darwin':  # macOS
            url = 'https://s3.amazonaws.com/plink2-assets/alpha5/plink2_mac_latest.zip'
            archive_name = 'plink2_mac.zip'
            binary_name = 'plink2'
        else:
            print(f"  PLINK2 auto-installation not supported on {self.system}")
            print(f"  Please install manually from: https://www.cog-genomics.org/plink/2.0/")
            return
        
        # Download
        download_path = self.bin_dir / archive_name
        print(f"  Downloading from {url}...")
        urllib.request.urlretrieve(url, download_path)
        
        # Extract
        print(f"  Extracting...")
        with zipfile.ZipFile(download_path, 'r') as zip_ref:
            zip_ref.extractall(self.bin_dir)
        
        # Make executable
        binary_path = self.bin_dir / binary_name
        binary_path.chmod(0o755)
        
        # Clean up
        download_path.unlink()
        
        # Verify
        result = subprocess.run([str(binary_path), '--version'], 
                              capture_output=True, text=True)
        if result.returncode == 0:
            print(f"  ✓ PLINK2 installed successfully at {binary_path}")
        else:
            print(f"  ✗ PLINK2 installation verification failed")
    
    def install_gcta(self):
        """Install GCTA."""
        print("Installing GCTA...")
        
        if self.system == 'Linux':
            url = 'https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-linux-kernel-3-x86_64.zip'
            archive_name = 'gcta_linux.zip'
            binary_name = 'gcta64'
            extract_dir = 'gcta-1.94.1-linux-kernel-3-x86_64'
        elif self.system == 'Darwin':  # macOS
            url = 'https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-MacOS-x86_64.zip'
            archive_name = 'gcta_mac.zip'
            binary_name = 'gcta64'
            extract_dir = 'gcta-1.94.1-MacOS-x86_64'
        else:
            print(f"  GCTA auto-installation not supported on {self.system}")
            print(f"  Please install manually from: https://yanglab.westlake.edu.cn/software/gcta/")
            return
        
        # Download
        download_path = self.bin_dir / archive_name
        print(f"  Downloading from {url}...")
        urllib.request.urlretrieve(url, download_path)
        
        # Extract
        print(f"  Extracting...")
        with zipfile.ZipFile(download_path, 'r') as zip_ref:
            zip_ref.extractall(self.bin_dir)
        
        # Move binary
        source_binary = self.bin_dir / extract_dir / binary_name
        dest_binary = self.bin_dir / binary_name
        
        if source_binary.exists():
            shutil.move(str(source_binary), str(dest_binary))
            dest_binary.chmod(0o755)
            
            # Clean up
            shutil.rmtree(self.bin_dir / extract_dir)
            download_path.unlink()
            
            # Verify
            result = subprocess.run([str(dest_binary), '--version'], 
                                  capture_output=True, text=True)
            if result.returncode == 0:
                print(f"  ✓ GCTA installed successfully at {dest_binary}")
            else:
                print(f"  ✗ GCTA installation verification failed")
        else:
            print(f"  ✗ GCTA binary not found after extraction")
    
    def install_r_packages(self):
        """Install R packages for PC-AiR."""
        print("Installing R packages for PC-AiR...")
        
        # Check if R is installed
        try:
            result = subprocess.run(['R', '--version'], 
                                  capture_output=True, text=True)
            if result.returncode != 0:
                print("  R is not installed. Skipping R package installation.")
                print("  Install R from: https://www.r-project.org/")
                return
        except FileNotFoundError:
            print("  R is not installed. Skipping R package installation.")
            print("  Install R from: https://www.r-project.org/")
            return
        
        # Create R installation script
        r_script = """
# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos="https://cloud.r-project.org")
}

# Install packages
packages <- c("GENESIS", "SNPRelate", "gdsfmt")

for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat(paste("Installing", pkg, "...\\n"))
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
        cat(paste(pkg, "already installed\\n"))
    }
}

cat("\\nR packages installation complete!\\n")
"""
        
        script_path = self.bin_dir / 'install_r_packages.R'
        with open(script_path, 'w') as f:
            f.write(r_script)
        
        # Run R script
        print("  Installing GENESIS, SNPRelate, and gdsfmt...")
        print("  (This may take several minutes...)")
        
        result = subprocess.run(['Rscript', str(script_path)],
                              capture_output=True, text=True)
        
        if result.returncode == 0:
            print("  ✓ R packages installed successfully")
        else:
            print("  ✗ R package installation failed")
            print(f"  Error: {result.stderr}")
        
        # Clean up
        script_path.unlink()


def read_requirements():
    """Read requirements from requirements.txt."""
    with open('requirements.txt') as f:
        return [line.strip() for line in f if line.strip() and not line.startswith('#')]


setup(
    name='edge-gwas',
    version=VERSION,
    description='EDGE: Encoding for Detecting Genetic Effects in GWAS',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Your Name',
    author_email='your.email@example.com',
    url='https://github.com/yourusername/edge-gwas',
    packages=find_packages(),
    install_requires=read_requirements(),
    python_requires='>=3.7',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    keywords='gwas genetics genomics bioinformatics edge encoding',
    entry_points={
        'console_scripts': [
            'edge-gwas-check-tools=edge_gwas.check_tools:main',
            'edge-gwas-install-tools=edge_gwas.install_tools:main',
        ],
    },
    cmdclass={
        'install': PostInstallCommand,
        'develop': PostDevelopCommand,
    },
    include_package_data=True,
)
