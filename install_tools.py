"""
Interactive installer for external tools.
"""

import sys
import subprocess
import platform
import urllib.request
import zipfile
import shutil
import os
from pathlib import Path


def main():
    """Interactive installation of external tools."""
    
    installer = ExternalToolsInstaller()
    
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
    
    installer.install_all()
    
    print("\nInstallation complete!")
    print("Run 'edge-gwas-check-tools' to verify installation.")
    
    return 0


class ExternalToolsInstaller:
    """Handle installation of external tools."""
    
    def __init__(self):
        self.system = platform.system()
        self.home = Path.home()
        self.bin_dir = self.home / '.local' / 'bin'
        self.bin_dir.mkdir(parents=True, exist_ok=True)
    
    def check_path(self):
        """Check if bin_dir is in PATH and provide instructions if not."""
        current_path = os.environ.get('PATH', '')
        bin_dir_str = str(self.bin_dir)
        
        if bin_dir_str not in current_path:
            print("\n" + "!"*70)
            print("IMPORTANT: Please add the following to your shell configuration:")
            print("!"*70)
            
            if self.system == 'Darwin':  # macOS
                shell_config = '~/.zshrc'
            else:  # Linux
                shell_config = '~/.bashrc'
            
            print(f'\nAdd this line to {shell_config}:')
            print(f'  export PATH="$HOME/.local/bin:$PATH"')
            print('\nThen run:')
            print(f'  source {shell_config}')
            print("!"*70 + "\n")
    
    def install_all(self):
        """Install all external tools."""
        print("\n" + "="*70)
        print("Installing external tools for EDGE-GWAS...")
        print("="*70)
        print(f"\nSystem: {self.system}")
        print(f"Installation directory: {self.bin_dir}\n")
        
        success_count = 0
        total_count = 3
        
        # Install PLINK2
        try:
            self.install_plink2()
            success_count += 1
        except Exception as e:
            print(f"\nWarning: PLINK2 installation failed: {e}")
            print("You can install manually from: https://www.cog-genomics.org/plink/2.0/\n")
        
        # Install GCTA
        try:
            self.install_gcta()
            success_count += 1
        except Exception as e:
            print(f"\nWarning: GCTA installation failed: {e}")
            print("You can install manually from: https://yanglab.westlake.edu.cn/software/gcta/\n")
        
        # Install R packages
        try:
            self.install_r_packages()
            success_count += 1
        except Exception as e:
            print(f"\nWarning: R packages installation failed: {e}")
            print("You can install manually in R:")
            print('  BiocManager::install(c("GENESIS", "SNPRelate", "gdsfmt"))\n')
        
        print("\n" + "="*70)
        print(f"External tools installation complete! ({success_count}/{total_count} successful)")
        print("="*70)
        
        # Check PATH
        self.check_path()
        
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
        
        try:
            urllib.request.urlretrieve(url, download_path)
        except Exception as e:
            raise Exception(f"Download failed: {e}")
        
        # Extract
        print(f"  Extracting...")
        try:
            with zipfile.ZipFile(download_path, 'r') as zip_ref:
                zip_ref.extractall(self.bin_dir)
        except Exception as e:
            download_path.unlink(missing_ok=True)
            raise Exception(f"Extraction failed: {e}")
        
        # Make executable
        binary_path = self.bin_dir / binary_name
        try:
            binary_path.chmod(0o755)
        except Exception as e:
            print(f"  Warning: Could not set executable permissions: {e}")
        
        # Clean up
        download_path.unlink(missing_ok=True)
        
        # Verify
        try:
            result = subprocess.run([str(binary_path), '--version'], 
                                  capture_output=True, text=True, timeout=5)
            if result.returncode == 0:
                print(f"  ✓ PLINK2 installed successfully at {binary_path}")
            else:
                print(f"  ✗ PLINK2 installation verification failed")
        except Exception as e:
            print(f"  Warning: Could not verify installation: {e}")
    
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
        
        try:
            urllib.request.urlretrieve(url, download_path)
        except Exception as e:
            raise Exception(f"Download failed: {e}")
        
        # Extract
        print(f"  Extracting...")
        try:
            with zipfile.ZipFile(download_path, 'r') as zip_ref:
                zip_ref.extractall(self.bin_dir)
        except Exception as e:
            download_path.unlink(missing_ok=True)
            raise Exception(f"Extraction failed: {e}")
        
        # Move binary
        source_binary = self.bin_dir / extract_dir / binary_name
        dest_binary = self.bin_dir / binary_name
        
        if source_binary.exists():
            try:
                shutil.move(str(source_binary), str(dest_binary))
                dest_binary.chmod(0o755)
                
                # Clean up
                shutil.rmtree(self.bin_dir / extract_dir, ignore_errors=True)
                download_path.unlink(missing_ok=True)
                
                # Verify
                result = subprocess.run([str(dest_binary), '--version'], 
                                      capture_output=True, text=True, timeout=5)
                if result.returncode == 0:
                    print(f"  ✓ GCTA installed successfully at {dest_binary}")
                else:
                    print(f"  ✗ GCTA installation verification failed")
            except Exception as e:
                print(f"  Warning: Installation completed but verification failed: {e}")
        else:
            print(f"  ✗ GCTA binary not found after extraction")
    
    def install_r_packages(self):
        """Install R packages for PC-AiR."""
        print("Installing R packages for PC-AiR...")
        
        # Check if R is installed
        try:
            result = subprocess.run(['R', '--version'], 
                                  capture_output=True, text=True, timeout=5)
            if result.returncode != 0:
                print("  R is not installed. Skipping R package installation.")
                print("  Install R from: https://www.r-project.org/")
                return
        except FileNotFoundError:
            print("  R is not installed. Skipping R package installation.")
            print("  Install R from: https://www.r-project.org/")
            return
        except Exception as e:
            print(f"  Could not check R installation: {e}")
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
        tryCatch({
            BiocManager::install(pkg, update = FALSE, ask = FALSE)
            cat(paste(pkg, "installed successfully\\n"))
        }, error = function(e) {
            cat(paste("Error installing", pkg, ":", e$message, "\\n"))
        })
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
        
        try:
            result = subprocess.run(['Rscript', str(script_path)],
                                  capture_output=True, text=True, timeout=600)  # 10 min timeout
            
            if result.returncode == 0:
                print("  ✓ R packages installed successfully")
            else:
                print("  ✗ R package installation failed")
                if result.stderr:
                    print(f"  Error: {result.stderr[:500]}")  # Print first 500 chars
        except subprocess.TimeoutExpired:
            print("  ✗ R package installation timed out (>10 minutes)")
        except Exception as e:
            print(f"  ✗ R package installation failed: {e}")
        
        # Clean up
        script_path.unlink(missing_ok=True)


if __name__ == '__main__':
    sys.exit(main())
