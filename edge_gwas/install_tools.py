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
        self.machine = platform.machine()  # arm64, x86_64, etc.
        self.home = Path.home()
        self.bin_dir = self.home / '.local' / 'bin'
        self.bin_dir.mkdir(parents=True, exist_ok=True)
        
        # Create SSL context that doesn't verify certificates (for downloads only)
        try:
            self.ssl_context = ssl._create_unverified_context()
        except:
            self.ssl_context = None
    
    def download_file(self, url, destination):
        """Download file handling SSL issues and verifying the download."""
        max_retries = 3
        
        for attempt in range(max_retries):
            try:
                if attempt > 0:
                    print(f"    Retry {attempt}/{max_retries}...")
                
                # Remove existing file if any
                if destination.exists():
                    destination.unlink()
                
                # Try with curl (most reliable for these URLs)
                result = subprocess.run(
                    ['curl', '-L', '-k', '-f', '--retry', '3', '--retry-delay', '2',
                     '-o', str(destination), url], 
                    capture_output=True,
                    text=True,
                    timeout=300
                )
                
                if result.returncode != 0:
                    print(f"    curl error: {result.stderr}")
                    if attempt < max_retries - 1:
                        continue
                    raise Exception(f"curl failed: {result.stderr}")
                
                # Verify file exists and has content
                if not destination.exists():
                    raise Exception("Downloaded file does not exist")
                
                file_size = destination.stat().st_size
                if file_size == 0:
                    raise Exception("Downloaded file is empty")
                
                # Verify it's a zip file
                try:
                    with zipfile.ZipFile(destination, 'r') as test_zip:
                        # Just test if we can read the file list
                        _ = test_zip.namelist()
                    print(f"    ✓ Downloaded and verified ({file_size:,} bytes)")
                    return
                except zipfile.BadZipFile:
                    print(f"    Downloaded file is not a valid zip ({file_size} bytes)")
                    if attempt < max_retries - 1:
                        continue
                    raise Exception("Downloaded file is not a valid zip file")
                
            except subprocess.TimeoutExpired:
                print(f"    Download timed out")
                if attempt < max_retries - 1:
                    continue
                raise Exception("Download timed out after multiple retries")
            
            except Exception as e:
                if attempt < max_retries - 1:
                    print(f"    Error: {e}")
                    continue
                raise
        
        raise Exception(f"Failed to download after {max_retries} attempts")
    
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
            print(f'  export PATH="$$HOME/.local/bin:$$PATH"')
            print('\nThen run:')
            print(f'  source {shell_config}')
            print("!"*70 + "\n")
    
    def choose_mac_architecture(self, tool_name):
        """Let user choose architecture for macOS."""
        print(f"\nDetected system: macOS ({self.machine})")
        
        if tool_name == "PLINK2":
            print("\nAvailable PLINK2 versions for macOS:")
            print("  1. ARM64 (M1/M2/M3 Macs - native, recommended for Apple Silicon)")
            print("  2. AVX2 (Intel Macs with AVX2 support - fastest for Intel)")
            print("  3. Standard (Intel Macs - compatible with all)")
            
            if self.machine == 'arm64':
                print("\nRecommended: Option 1 (ARM64)")
                default = '1'
            else:
                print("\nRecommended: Option 2 (AVX2) if your Intel Mac supports it, otherwise 3")
                default = '2'
            
            choice = input(f"\nSelect version [1/2/3, default={default}]: ").strip()
            
            if not choice:
                choice = default
            
            if choice == '1':
                return 'arm64', 'https://s3.amazonaws.com/plink2-assets/plink2_mac_arm64_20251205.zip'
            elif choice == '2':
                return 'avx2', 'https://s3.amazonaws.com/plink2-assets/plink2_mac_avx2_20251205.zip'
            else:
                return 'standard', 'https://s3.amazonaws.com/plink2-assets/plink2_mac_20251205.zip'
        
        elif tool_name == "GCTA":
            print("\nAvailable GCTA versions for macOS:")
            print("  1. ARM64 (M1/M2/M3 Macs - native, v1.95.0)")
            print("  2. x86_64 (Intel Macs or Rosetta 2 - stable, v1.94.1)")
            print("  3. Linux version (fallback - works via compatibility)")
            
            if self.machine == 'arm64':
                print("\nRecommended: Option 1 (ARM64) for native performance")
                print("Note: If option 1 fails, will automatically try Linux version")
                default = '1'
            else:
                print("\nRecommended: Option 2 (x86_64)")
                print("Note: If option 2 fails, will automatically try Linux version")
                default = '2'
            
            choice = input(f"\nSelect version [1/2/3, default={default}]: ").strip()
            
            if not choice:
                choice = default
            
            if choice == '1':
                return 'arm64', 'https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.95.0-macOS-arm64.zip', 'gcta-1.95.0-macOS-arm64'
            elif choice == '2':
                return 'x86_64', 'https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-macOS-x86_64.zip', 'gcta-1.94.1-macOS-x86_64'
            else:
                return 'linux', 'https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.95.0-linux-kernel-3-x86_64.zip', 'gcta-1.95.0-linux-kernel-3-x86_64'
    
    def install_all(self):
        """Install all external tools."""
        print("\n" + "="*70)
        print("Installing external tools for EDGE-GWAS...")
        print("="*70)
        print(f"\nSystem: {self.system}")
        print(f"Architecture: {self.machine}")
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
            arch_type = 'linux'
        elif self.system == 'Darwin':  # macOS
            arch_type, url = self.choose_mac_architecture("PLINK2")
            print(f"  Selected: PLINK2 macOS {arch_type}")
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
                version_info = result.stdout.strip().split('\n')[0]
                print(f"  ✓ PLINK2 installed successfully: {version_info}")
            else:
                print(f"  ⚠ PLINK2 installed but verification failed")
        except Exception as e:
            print(f"  ⚠ PLINK2 installed at {binary_path} (verification skipped: {e})")
    
    def install_gcta(self):
        """Install GCTA with fallback to Linux version on macOS."""
        print("\nInstalling GCTA...")
        
        if self.system == 'Linux':
            url = 'https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-linux-kernel-3-x86_64.zip'
            extract_dir = 'gcta-1.94.1-linux-kernel-3-x86_64'
            arch_type = 'linux'
            use_fallback = False
        elif self.system == 'Darwin':  # macOS
            arch_type, url, extract_dir = self.choose_mac_architecture("GCTA")
            print(f"  Selected: GCTA macOS {arch_type}")
            use_fallback = (arch_type != 'linux')  # Only use fallback if not already Linux version
        else:
            print(f"  GCTA auto-installation not supported on {self.system}")
            raise Exception(f"Unsupported OS: {self.system}")
        
        # Try the selected version
        try:
            self._install_gcta_version(url, extract_dir, arch_type)
        except Exception as e:
            if use_fallback and self.system == 'Darwin':
                print(f"\n  macOS version failed: {e}")
                print(f"  Trying Linux version as fallback...")
                # Fallback to Linux version
                fallback_url = 'https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-linux-kernel-3-x86_64.zip'
                fallback_dir = 'gcta-1.94.1-linux-kernel-3-x86_64'
                self._install_gcta_version(fallback_url, fallback_dir, 'linux-fallback')
            else:
                raise
    
    def _install_gcta_version(self, url, extract_dir, arch_type):
    """Install a specific GCTA version."""
    zip_path = self.bin_dir / 'gcta.zip'
    print(f"  Downloading from {url}...")
    
    try:
        self.download_file(url, zip_path)
    except Exception as e:
        print(f"\n  Download failed: {e}")
        raise
    
    print(f"  Extracting...")
    try:
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            # List contents for debugging
            contents = zip_ref.namelist()
            print(f"    Archive contains {len(contents)} files")
            zip_ref.extractall(self.bin_dir)
    except zipfile.BadZipFile as e:
        print(f"  Error: Invalid zip file")
        print(f"  File size: {zip_path.stat().st_size} bytes")
        zip_path.unlink(missing_ok=True)
        raise Exception(f"Invalid zip file: {e}")
    
    # Look for the binary in the extracted directory
    # Try multiple possible binary locations and names
    possible_binaries = [
        # ARM64 version has bin/gcta64
        self.bin_dir / extract_dir / 'bin' / 'gcta64',
        # Root directory locations
        self.bin_dir / extract_dir / 'gcta64',
        self.bin_dir / extract_dir / 'gcta-1.94.1',  # macOS x86_64
        self.bin_dir / extract_dir / 'gcta-1.95.0',
        self.bin_dir / extract_dir / 'gcta',
        self.bin_dir / extract_dir / 'gcta_1.94.1',
        self.bin_dir / extract_dir / 'gcta_1.95.0',
    ]
    
    source_binary = None
    for candidate in possible_binaries:
        if candidate.exists() and candidate.is_file():
            # Check if it's actually an executable binary (not a text file)
            try:
                # Try to read first few bytes
                with open(candidate, 'rb') as f:
                    header = f.read(4)
                    # Check for ELF magic number (Linux) or Mach-O (macOS)
                    # Mach-O: \xfe\xed\xfa\xce (32-bit), \xfe\xed\xfa\xcf (64-bit)
                    # Mach-O (reverse): \xce\xfa\xed\xfe, \xcf\xfa\xed\xfe
                    # ELF: \x7fELF
                    if header[:4] in [b'\x7fELF', b'\xcf\xfa\xed\xfe', b'\xce\xfa\xed\xfe', 
                                     b'\xfe\xed\xfa\xce', b'\xfe\xed\xfa\xcf']:
                        source_binary = candidate
                        print(f"    Found binary: {candidate.relative_to(self.bin_dir)}")
                        break
            except Exception as e:
                continue
    
    if not source_binary:
        # Debug: show what was extracted
        print(f"  Debug: Looking for binary in {extract_dir}")
        if (self.bin_dir / extract_dir).exists():
            print(f"  Contents of {extract_dir}:")
            for item in (self.bin_dir / extract_dir).iterdir():
                if item.is_file():
                    size = item.stat().st_size
                    print(f"    - {item.name} (file, {size:,} bytes)")
                elif item.is_dir():
                    print(f"    - {item.name}/ (dir)")
                    # If it's a bin directory, show its contents
                    if item.name == 'bin':
                        print(f"      Contents of bin/:")
                        for subitem in item.iterdir():
                            if subitem.is_file():
                                subsize = subitem.stat().st_size
                                print(f"        - {subitem.name} (file, {subsize:,} bytes)")
        else:
            print(f"  Directory not found: {extract_dir}")
            print(f"  Available directories in {self.bin_dir}:")
            for item in self.bin_dir.iterdir():
                if item.is_dir() and 'gcta' in item.name.lower():
                    print(f"    - {item.name}")
        
        # Clean up before raising exception
        zip_path.unlink(missing_ok=True)
        shutil.rmtree(self.bin_dir / extract_dir, ignore_errors=True)
        raise Exception(f"GCTA binary not found in extracted directory")
    
    dest_binary = self.bin_dir / 'gcta64'
    
    print(f"  Moving {source_binary.name} to {dest_binary.name}")
    
    # Remove existing binary if present
    if dest_binary.exists():
        dest_binary.unlink()
    
    shutil.move(str(source_binary), str(dest_binary))
    dest_binary.chmod(0o755)
    
    # Clean up
    shutil.rmtree(self.bin_dir / extract_dir, ignore_errors=True)
    zip_path.unlink(missing_ok=True)
    
    # Verify installation
    try:
        result = subprocess.run([str(dest_binary), '--version'], 
                              capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            version_info = result.stdout.strip().split('\n')[0]
            print(f"  ✓ GCTA installed successfully ({arch_type}): {version_info}")
        else:
            print(f"  ⚠ GCTA installed but verification failed")
            print(f"     Return code: {result.returncode}")
            if result.stderr:
                print(f"     Error: {result.stderr[:200]}")
    except OSError as e:
        if 'Exec format error' in str(e):
            print(f"  ✗ Error: This binary is not compatible with your system")
            print(f"     Binary type mismatch (wrong architecture)")
            dest_binary.unlink()  # Remove the incompatible binary
            raise Exception("Incompatible binary format")
        else:
            print(f"  ⚠ GCTA installed at {dest_binary} (verification skipped: {e})")
    
    def install_r_packages(self):
        """Install R packages for PC-AiR."""
        print("\nInstalling R packages for PC-AiR...")
        
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
