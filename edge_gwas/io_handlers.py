"""
Input/Output handlers for EDGE GWAS analysis.
"""

import os
import urllib.request
import pandas as pd
import numpy as np
from typing import Optional, Dict, List
import logging

logger = logging.getLogger(__name__)


def save_results(
    gwas_df: pd.DataFrame,
    alpha_df: Optional[pd.DataFrame] = None,
    output_prefix: str = 'edge_gwas',
    save_alpha: bool = True
) -> Dict[str, str]:
    """
    Save EDGE GWAS results to files.
    
    Args:
        gwas_df: GWAS results DataFrame
        alpha_df: Alpha values DataFrame
        output_prefix: Prefix for output files
        save_alpha: Whether to save alpha values
        
    Returns:
        Dictionary with output file paths
    """
    output_files = {}
    
    # Save GWAS results
    gwas_file = f"{output_prefix}_gwas_results.tsv"
    gwas_df.to_csv(gwas_file, sep='\t', index=False, na_rep='NA')
    output_files['gwas'] = gwas_file
    logger.info(f"GWAS results saved to {gwas_file}")
    
    # Save alpha values
    if save_alpha and alpha_df is not None:
        alpha_file = f"{output_prefix}_alpha_values.tsv"
        alpha_df.to_csv(alpha_file, sep='\t', index=False, na_rep='NA')
        output_files['alpha'] = alpha_file
        logger.info(f"Alpha values saved to {alpha_file}")
    
    return output_files


def load_alpha_values(alpha_file: str) -> pd.DataFrame:
    """
    Load pre-calculated alpha values.
    
    Args:
        alpha_file: Path to alpha values file
        
    Returns:
        DataFrame with alpha values
    """
    logger.info(f"Loading alpha values from {alpha_file}")
    
    # Try different separators
    for sep in ['\t', ',', ' ']:
        try:
            alpha_df = pd.read_csv(alpha_file, sep=sep)
            if 'variant_id' in alpha_df.columns and 'alpha_value' in alpha_df.columns:
                logger.info(f"Successfully loaded {len(alpha_df)} alpha values")
                return alpha_df
        except:
            continue
    
    raise ValueError(f"Could not load alpha values from {alpha_file}")


def format_gwas_output(
    gwas_df: pd.DataFrame,
    include_alpha: bool = True,
    sort_by: str = 'pval'
) -> pd.DataFrame:
    """
    Format GWAS output for publication/reporting.
    
    Args:
        gwas_df: GWAS results DataFrame
        include_alpha: Include alpha-related columns
        sort_by: Column to sort by
        
    Returns:
        Formatted DataFrame
    """
    # Define column order
    base_cols = ['variant_id', 'coef', 'std_err', 'stat', 'pval', 
                 'conf_int_low', 'conf_int_high', 'n_samples']
    
    if include_alpha:
        alpha_cols = ['alpha_value', 'ref_allele', 'alt_allele', 'eaf']
        cols = base_cols[:1] + alpha_cols + base_cols[1:]
    else:
        cols = base_cols
    
    # Select and order columns
    available_cols = [col for col in cols if col in gwas_df.columns]
    formatted_df = gwas_df[available_cols].copy()
    
    # Sort results
    if sort_by in formatted_df.columns:
        formatted_df = formatted_df.sort_values(sort_by)
    
    return formatted_df


def create_summary_report(
    gwas_df: pd.DataFrame,
    alpha_df: Optional[pd.DataFrame] = None,
    significance_threshold: float = 5e-8,
    output_file: Optional[str] = None
) -> str:
    """
    Create a summary report of EDGE GWAS analysis.
    
    Args:
        gwas_df: GWAS results DataFrame
        alpha_df: Alpha values DataFrame
        significance_threshold: P-value threshold for significance
        output_file: Optional file to save report
        
    Returns:
        Summary report as string
    """
    report_lines = []
    report_lines.append("=" * 80)
    report_lines.append("EDGE GWAS Analysis Summary Report")
    report_lines.append("=" * 80)
    report_lines.append("")
    
    # GWAS summary
    report_lines.append("GWAS Results Summary:")
    report_lines.append(f"  Total variants tested: {len(gwas_df)}")
    
    if 'pval' in gwas_df.columns:
        significant = gwas_df[gwas_df['pval'] < significance_threshold]
        report_lines.append(f"  Genome-wide significant hits (p < {significance_threshold}): {len(significant)}")
        
        if len(significant) > 0:
            report_lines.append("")
            report_lines.append("Top significant variants:")
            top_variants = significant.nsmallest(10, 'pval')
            for idx, row in top_variants.iterrows():
                report_lines.append(f"    {row['variant_id']}: p={row['pval']:.2e}, "
                                  f"coef={row['coef']:.4f}")
        
        # Genomic inflation
        from .utils import calculate_genomic_inflation
        lambda_gc = calculate_genomic_inflation(gwas_df['pval'])
        report_lines.append("")
        report_lines.append(f"  Genomic inflation factor (λ): {lambda_gc:.4f}")
    
    # Alpha values summary
    if alpha_df is not None:
        report_lines.append("")
        report_lines.append("Alpha Values Summary:")
        report_lines.append(f"  Total variants with alpha values: {len(alpha_df)}")
        
        if 'alpha_value' in alpha_df.columns:
            valid_alpha = alpha_df['alpha_value'].dropna()
            report_lines.append(f"  Valid alpha values: {len(valid_alpha)}")
            report_lines.append(f"  Mean alpha: {valid_alpha.mean():.4f}")
            report_lines.append(f"  Median alpha: {valid_alpha.median():.4f}")
            report_lines.append(f"  Alpha range: [{valid_alpha.min():.4f}, {valid_alpha.max():.4f}]")
            
            # Distribution of alpha values
            alpha_ranges = {
                'Near additive (0.4-0.6)': ((valid_alpha >= 0.4) & (valid_alpha <= 0.6)).sum(),
                'Dominant-like (> 0.6)': (valid_alpha > 0.6).sum(),
                'Recessive-like (< 0.4)': (valid_alpha < 0.4).sum()
            }
            
            report_lines.append("")
            report_lines.append("  Alpha value distribution:")
            for range_name, count in alpha_ranges.items():
                pct = count / len(valid_alpha) * 100
                report_lines.append(f"    {range_name}: {count} ({pct:.1f}%)")
    
    report_lines.append("")
    report_lines.append("=" * 80)
    
    report = "\n".join(report_lines)
    
    # Save to file if requested
    if output_file:
        with open(output_file, 'w') as f:
            f.write(report)
        logger.info(f"Summary report saved to {output_file}")
    
    return report


def download_test_files(
    output_dir: str = 'tests',
    version: str = 'v0.1.1',
    overwrite: bool = False,
    verbose: bool = True
) -> dict:
    """
    Download test files from GitHub repository.
    
    Args:
        output_dir: Directory to save test files (default: 'tests')
        version: GitHub release version tag (default: 'v0.1.1')
        overwrite: If True, overwrite existing files (default: False)
        verbose: Print download progress (default: True)
    
    Returns:
        Dictionary with download results:
        {
            'downloaded': List of successfully downloaded files,
            'skipped': List of skipped (already existing) files,
            'failed': List of failed downloads
        }
    
    Examples:
        >>> from edge_gwas.io_handler import download_test_files
        >>> results = download_test_files()
        >>> print(f"Downloaded {len(results['downloaded'])} files")
        
        >>> # Force re-download
        >>> download_test_files(overwrite=True)
        
        >>> # Download to custom directory
        >>> download_test_files(output_dir='my_data/tests')
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    if verbose:
        logger.info(f"Downloading test files to {output_dir}")
        print(f"Downloading test files to: {output_dir}")
    
    # Base URL for raw files
    base_url = f"https://raw.githubusercontent.com/nicenzhou/edge-gwas/{version}/tests/"
    files = ['test.bed', 'test.bim', 'test.fam', 'test.phen', 'test.vcf']
    
    results = {
        'downloaded': [],
        'skipped': [],
        'failed': []
    }
    
    for filename in files:
        filepath = os.path.join(output_dir, filename)
        
        # Check if file exists
        if os.path.exists(filepath) and not overwrite:
            if verbose:
                print(f"✓ Already exists: {filename}")
            logger.debug(f"Skipped existing file: {filepath}")
            results['skipped'].append(filename)
            continue
        
        # Download file
        try:
            url = base_url + filename
            
            if verbose:
                print(f"⬇ Downloading: {filename}...", end=' ')
            
            urllib.request.urlretrieve(url, filepath)
            
            # Verify download
            if os.path.exists(filepath) and os.path.getsize(filepath) > 0:
                if verbose:
                    size_mb = os.path.getsize(filepath) / (1024 * 1024)
                    print(f"✓ ({size_mb:.2f} MB)")
                logger.info(f"Downloaded {filename} from {url}")
                results['downloaded'].append(filename)
            else:
                raise ValueError("Downloaded file is empty or doesn't exist")
                
        except Exception as e:
            if verbose:
                print(f"✗ Error: {e}")
            logger.error(f"Failed to download {filename}: {e}")
            results['failed'].append(filename)
    
    # Print summary
    if verbose:
        print(f"\n{'='*50}")
        print(f"Summary:")
        print(f"  ✓ Downloaded: {len(results['downloaded'])}")
        print(f"  ⊙ Skipped: {len(results['skipped'])}")
        print(f"  ✗ Failed: {len(results['failed'])}")
        print(f"{'='*50}")
    
    if results['failed']:
        logger.warning(f"Failed to download: {', '.join(results['failed'])}")
    
    return results
