"""
EDGE-GWAS: Encoding for Detecting Genetic Effects in GWAS
"""

__version__ = "0.1.2"
__author__ = "Jiayan Zhou"
__email__ = "jyzhou@stanford.edu"

# Import main classes
from .core import EDGEAnalysis

# Import utility functions
from .utils import (
    load_plink_data,
    load_pgen_data,
    load_bgen_data,
    load_vcf_data,
    prepare_phenotype_data,
    stratified_train_test_split,
    filter_variants_by_maf,
    filter_variants_by_missing,
    filter_samples_by_call_rate,
    calculate_hwe_pvalues,
    filter_variants_by_hwe,
    # Population structure functions
    calculate_pca_sklearn,
    calculate_pca_plink,
    calculate_pca_pcair,
    calculate_grm_gcta,
    load_grm_gcta,
    attach_pcs_to_phenotype,
    get_pc_covariate_list,
    identify_related_samples,
    filter_related_samples,
    calculate_genomic_inflation,
)

# Import IO handler functions
from .io_handlers import (
    save_results,
    load_alpha_values,
    format_gwas_output,
    create_summary_report,
    download_test_files,
)

# Import visualize functions
from .visualize import (
    manhattan_plot,
    qq_plot,
    plot_alpha_distribution,
)


__all__ = [
    # Core
    'EDGEAnalysis',
    
    # Data loading
    'load_plink_data',
    'load_pgen_data',
    'load_bgen_data',
    'load_vcf_data',
    'prepare_phenotype_data',
    
    # Data processing
    'stratified_train_test_split',
    'filter_variants_by_maf',
    'filter_variants_by_missing',
    'filter_samples_by_call_rate',
    'calculate_hwe_pvalues',
    'filter_variants_by_hwe',
    
    # IO handlers
    'save_results',
    'load_alpha_values',
    'format_gwas_output',
    'create_summary_report',
    'download_test_files',
    
    # Population structure
    'calculate_pca_sklearn',
    'calculate_pca_plink',
    'calculate_pca_pcair',
    'calculate_grm_gcta',
    'load_grm_gcta',
    'attach_pcs_to_phenotype',
    'get_pc_covariate_list',
    'identify_related_samples',
    'filter_related_samples',
    'calculate_genomic_inflation',
    
    # Visualization
    'manhattan_plot',
    'qq_plot',
    'plot_alpha_distribution',
]


def check_external_tools():
    """
    Check if external tools are installed.
    
    Returns:
        dict: Status of each tool (True/False)
    """
    from .check_tools import check_tool, check_r_package
    
    status = {
        'plink2': check_tool('PLINK2', 'plink2'),
        'gcta': check_tool('GCTA', 'gcta64') or check_tool('GCTA', 'gcta'),
        'r': check_tool('R', 'R', '--version'),
        'genesis': check_r_package('GENESIS'),
        'snprelate': check_r_package('SNPRelate'),
        'gdsfmt': check_r_package('gdsfmt'),
    }
    
    return status


def get_info():
    """
    Get package information.
    
    Returns:
        dict: Package information including version, tools status, etc.
    """
    import sys
    import platform
    
    info = {
        'version': __version__,
        'python_version': sys.version,
        'platform': platform.platform(),
        'external_tools': check_external_tools(),
    }
    
    return info


def print_info():
    """Print package information."""
    info = get_info()
    
    print("\n" + "="*70)
    print(f"EDGE-GWAS version {info['version']}")
    print("="*70)
    print(f"\nPython: {info['python_version']}")
    print(f"Platform: {info['platform']}")
    print("\nExternal Tools:")
    print("-" * 70)
    
    for tool, installed in info['external_tools'].items():
        status = "✓ Installed" if installed else "✗ Not installed"
        print(f"  {tool}: {status}")
    
    print("="*70 + "\n")
