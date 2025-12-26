"""Interactive installer for external tools."""
import sys
import subprocess
import platform
import urllib.request
import zipfile
import shutil
import os
import ssl
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
        
        # Create SSL context that doesn't verify certificates (for downloads only)
        try:
            self.ssl_context = ssl._create_unverified_context()
        except:
            self.ssl_context = None
    
    def download_file(self, url, destination):
        """Download file handling SSL issues."""
        try:
            # Try with SSL context first
            if self.ssl_context:
                with urllib.request.urlopen(url, context=self.ssl_context) as response:
                    with open(destination, 'wb') as out_file:
                        out_file.write(response.read())
            else:
                urllib.request.urlretrieve(url, destination)
        except Exception as e:
            # Fallback: try using curl command
            print(f"  urllib failed, trying curl... ({e})")
            result = subprocess.run(['curl', '-L', '-k', '-o', str(destination), url], 
                                  capture_output=True, timeout=300)
            if result.returncode != 0:
                raise Exception(f"Download failed with curl: {result.stderr.decode()}")
    
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
            url = 'https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20251205.zip'
        elif self.system == 'Darwin':  # macOS
            url = 'https://s3.amazonaws.com/plink2-assets/plink2_mac_20251205.zip'
        else:
            print(f"  PLINK2 auto-installation not supported on {self.system}")
            raise Exception(f"Unsupported OS: {self.system}")
        
        zip_path = self.bin_dir / 'plink2.zip'
        print(f"  Downloading from {url}...")
        
        self.download_file(url, zip_path)
        
        print(f"  Extracting...")
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(self.bin_dir)
        
        binary_path = self.bin_dir / 'plink2'
        binary_path.chmod(0o755)
        zip_path.unlink()
        
        # Verify installation
        try:
            result = subprocess.run([str(binary_path), '--version'], 
                                  capture_output=True, text=True, timeout=5)
            if result.returncode == 0:
                print(f"  ✓ PLINK2 installed successfully at {binary_path}")
            else:
                print(f"  ⚠ PLINK2 installed but verification failed")
        except Exception as e:
            print(f"  ⚠ PLINK2 installed at {binary_path} (verification skipped: {e})")
    
    def install_gcta(self):
        """Install GCTA."""
        print("Installing GCTA...")
        
        if self.system == 'Linux':
            url = 'https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.95.0-linux-kernel-3-x86_64.zip'
            extract_dir = 'gcta-1.95.0-linux-kernel-3-x86_64'
        elif self.system == 'Darwin':  # macOS
            url = 'https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.95.0-macOS-arm64.zip'
            extract_dir = 'gcta-1.95.0-macOS-arm64'
        else:
            print(f"  GCTA auto-installation not supported on {self.system}")
            raise Exception(f"Unsupported OS: {self.system}")
        
        zip_path = self.bin_dir / 'gcta.zip'
        print(f"  Downloading from {url}...")
        
        self.download_file(url, zip_path)
        
        print(f"  Extracting...")
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(self.bin_dir)
        
        # Try both gcta64 and gcta (newer versions might use different names)
        source_binary = self.bin_dir / extract_dir / 'gcta64'
        if not source_binary.exists():
            source_binary = self.bin_dir / extract_dir / 'gcta'
        
        dest_binary = self.bin_dir / 'gcta64'
        
        if not source_binary.exists():
            raise Exception(f"GCTA binary not found after extraction in {extract_dir}")
        
        shutil.move(str(source_binary), str(dest_binary))
        dest_binary.chmod(0o755)
        shutil.rmtree(self.bin_dir / extract_dir, ignore_errors=True)
        zip_path.unlink()
        
        # Verify installation
        try:
            result = subprocess.run([str(dest_binary), '--version'], 
                                  capture_output=True, text=True, timeout=5)
            if result.returncode == 0:
                print(f"  ✓ GCTA installed successfully at {dest_binary}")
            else:
                print(f"  ⚠ GCTA installed but verification failed")
        except Exception as e:
            print(f"  ⚠ GCTA installed at {dest_binary} (verification skipped: {e})")
    
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
                raise Exception("R not installed")
        except FileNotFoundError:
            print("  R is not installed. Skipping R package installation.")
            print("  Install R from: https://www.r-project.org/")
            raise Exception("R not found")
        
        r_script = '''
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos="https://cloud.r-project.org")
}

packages <- c("GENESIS", "SNPRelate", "gdsfmt")

for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        cat(paste("Installing", pkg, "...\\n"))
        tryCatch({
            BiocManager::install(pkg, update = FALSE, ask = FALSE)
            cat(paste("✓", pkg, "installed\\n"))
        }, error = function(e) {
            cat(paste("✗", pkg, "installation failed:", e$message, "\\n"))
        })
    } else {
        cat(paste("✓", pkg, "already installed\\n"))
    }
}
cat("\\nR packages installation complete!\\n")
'''
        
        script_file = self.bin_dir / 'install_r_packages.R'
        script_file.write_text(r_script)
        
        print("  Installing GENESIS, SNPRelate, and gdsfmt...")
        print("  (This may take several minutes...)")
        
        try:
            result = subprocess.run(['Rscript', str(script_file)], 
                                  capture_output=True, text=True, timeout=600)
            print(result.stdout)
            if result.returncode == 0:
                print("  ✓ R packages installed successfully")
            else:
                print(f"  ⚠ R packages installation completed with warnings")
                if result.stderr:
                    print(f"  Errors: {result.stderr[:500]}")
        finally:
            script_file.unlink(missing_ok=True)


if __name__ == '__main__':
    sys.exit(main())
