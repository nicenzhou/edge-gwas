"""Check if external tools are properly installed."""
import subprocess
import sys


def check_tool(name, command, version_flag='--version'):
    """Check if a tool is installed and accessible."""
    try:
        result = subprocess.run(
            [command, version_flag],
            capture_output=True,
            text=True,
            timeout=10
        )
        if result.returncode == 0:
            version = result.stdout.strip().split('\n')[0]
            print(f"✓ {name}: {version}")
            return True
        else:
            print(f"✗ {name}: Installed but returned error")
            return False
    except FileNotFoundError:
        print(f"✗ {name}: Not found")
        return False
    except subprocess.TimeoutExpired:
        print(f"✗ {name}: Timeout")
        return False
    except Exception as e:
        print(f"✗ {name}: Error - {e}")
        return False


def check_r_package(package_name):
    """Check if an R package is installed."""
    r_code = f'if (requireNamespace("{package_name}", quietly = TRUE)) {{ cat("installed") }} else {{ cat("missing") }}'
    
    try:
        result = subprocess.run(
            ['Rscript', '-e', r_code],
            capture_output=True,
            text=True,
            timeout=10
        )
        
        if 'installed' in result.stdout:
            print(f"✓ R package {package_name}: Installed")
            return True
        else:
            print(f"✗ R package {package_name}: Not installed")
            return False
    except FileNotFoundError:
        print(f"✗ R package {package_name}: R not found")
        return False
    except subprocess.TimeoutExpired:
        print(f"✗ R package {package_name}: Timeout")
        return False
    except Exception as e:
        print(f"✗ R package {package_name}: Error - {e}")
        return False


def main():
    """Main function to check all tools."""
    print("\n" + "="*70)
    print("EDGE-GWAS External Tools Check")
    print("="*70 + "\n")
    
    all_ok = True
    
    # Check Python packages
    print("Python Packages:")
    print("-" * 70)
    
    python_packages = [
        'numpy', 'pandas', 'scipy', 'statsmodels', 
        'sklearn', 'matplotlib', 'pandas_plink'
    ]
    
    for package in python_packages:
        try:
            if package == 'sklearn':
                __import__('sklearn')
            else:
                __import__(package)
            print(f"✓ {package}: Installed")
        except ImportError:
            print(f"✗ {package}: Not installed")
            all_ok = False
    
    print()
    
    # Check external tools
    print("External Tools:")
    print("-" * 70)
    
    plink2_ok = check_tool('PLINK2', 'plink2')
    gcta_ok = check_tool('GCTA', 'gcta64') or check_tool('GCTA', 'gcta')
    
    print()
    
    # Check R
    print("R and Packages:")
    print("-" * 70)
    
    r_ok = check_tool('R', 'R', '--version')
    
    if r_ok:
        genesis_ok = check_r_package('GENESIS')
        snprelate_ok = check_r_package('SNPRelate')
        gdsfmt_ok = check_r_package('gdsfmt')
    else:
        print("  Skipping R package checks (R not installed)")
        genesis_ok = snprelate_ok = gdsfmt_ok = False
    
    print()
    print("="*70)
    
    # Summary
    core_ok = all([plink2_ok, gcta_ok, r_ok, genesis_ok, snprelate_ok, gdsfmt_ok])
    
    if core_ok and all_ok:
        print("✓ All tools and packages are properly installed!")
        print("\nYou can now use all EDGE-GWAS features including:")
        print("  - Basic EDGE analysis")
        print("  - PCA calculation with PLINK2")
        print("  - GRM calculation with GCTA")
        print("  - PC-AiR analysis with R/GENESIS")
    else:
        print("⚠ Some tools or packages are missing")
        print("\nMissing components will limit functionality:")
        
        if not plink2_ok:
            print("  - PLINK2: Cannot use calculate_pca_plink()")
        if not gcta_ok:
            print("  - GCTA: Cannot use calculate_grm_gcta()")
        if not (r_ok and genesis_ok and snprelate_ok and gdsfmt_ok):
            print("  - R/GENESIS: Cannot use calculate_pca_pcair()")
        
        print("\nTo install missing tools, run:")
        print("  edge-gwas-install-tools")
        print("\nOr install manually following the guide at:")
        print("  https://github.com/nicenzhou/edge-gwas/blob/main/INSTALL.md")
    
    print("="*70 + "\n")
    
    return 0 if (core_ok and all_ok) else 1


if __name__ == '__main__':
    sys.exit(main())
