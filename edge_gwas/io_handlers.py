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
        except Exception:
            continue
    
    raise ValueError(f"Could not load alpha values from {alpha_file}")


def format_gwas_output_for_locuszoom(
    gwas_df: pd.DataFrame,
    include_alpha: bool = True,
    sort_by: str = 'pval',
    format_for_locuszoom: bool = False
) -> pd.DataFrame:
    """
    Format GWAS output for publication/reporting or LocusZoom upload.
    
    Args:
        gwas_df: GWAS results DataFrame
        include_alpha: Include alpha-related columns (default: True)
        sort_by: Column to sort by (default: 'pval')
        format_for_locuszoom: If True, format for LocusZoom upload with required columns
                             (default: False)
        
    Returns:
        Formatted DataFrame
        
    Note:
        For LocusZoom format, the output will be tab-delimited with columns:
        - chrom (chromosome: 1-22, X, Y, M, MT)
        - pos (position)
        - ref (reference allele)
        - alt (alternate/effect allele)
        - pval (p-value)
        - beta (effect size, assuming alt is effect allele)
        - se (standard error)
        - eaf (effect allele frequency, oriented to alt allele)
        - alpha_value (optional, if include_alpha=True)
        
        The file should be sorted by chrom and pos, compressed with bgzip,
        and indexed with tabix for optimal LocusZoom performance.
    """
    if format_for_locuszoom:
        # LocusZoom required/recommended columns
        locuszoom_cols = {
            'chrom': 'chrom',      # Chromosome
            'pos': 'pos',          # Position
            'ref': 'ref_allele',   # Reference allele (source column name)
            'alt': 'alt_allele',   # Alternate allele (source column name)
            'pval': 'pval',        # P-value
            'beta': 'coef',        # Effect size (assuming alt is effect allele)
            'se': 'std_err',       # Standard error
            'eaf': 'eaf'           # Effect allele frequency (alt allele)
        }
        
        # Check for required columns
        required = ['chrom', 'pos', 'pval']
        missing = [col for col in required if locuszoom_cols.get(col) not in gwas_df.columns]
        if missing:
            raise ValueError(f"Missing required columns for LocusZoom format: {missing}")
        
        # Build column mapping
        col_mapping = {}
        for lz_col, source_col in locuszoom_cols.items():
            if source_col in gwas_df.columns:
                col_mapping[source_col] = lz_col
        
        # Select and rename columns
        available_source_cols = [source_col for source_col in col_mapping.keys() 
                                if source_col in gwas_df.columns]
        formatted_df = gwas_df[available_source_cols].copy()
        formatted_df = formatted_df.rename(columns=col_mapping)
        
        # Add alpha_value if requested and available
        if include_alpha and 'alpha_value' in gwas_df.columns:
            formatted_df['alpha_value'] = gwas_df['alpha_value']
        
        # Ensure proper column order for LocusZoom
        ordered_cols = ['chrom', 'pos', 'ref', 'alt', 'pval', 'beta', 'se', 'eaf']
        if include_alpha and 'alpha_value' in formatted_df.columns:
            ordered_cols.append('alpha_value')
        
        # Only include columns that exist
        final_cols = [col for col in ordered_cols if col in formatted_df.columns]
        formatted_df = formatted_df[final_cols]
        
        # Sort by chromosome and position (required for LocusZoom)
        # Convert chromosome to numeric for proper sorting
        formatted_df['chrom_sort'] = formatted_df['chrom'].apply(
            lambda x: 23 if str(x).upper() == 'X' else 
                     24 if str(x).upper() == 'Y' else 
                     25 if str(x).upper() in ['M', 'MT'] else 
                     int(x) if str(x).isdigit() else 99
        )
        formatted_df = formatted_df.sort_values(['chrom_sort', 'pos'])
        formatted_df = formatted_df.drop('chrom_sort', axis=1)
        
        return formatted_df
    
    else:
        # Standard format for publication/reporting
        base_cols = ['chrom', 'pos', 'variant_id', 'coef', 'std_err', 'stat', 'pval', 
                     'conf_int_low', 'conf_int_high', 'n_samples']
        
        if include_alpha:
            alpha_cols = ['alpha_value', 'ref_allele', 'alt_allele', 'eaf', 'maf']
            # Insert alpha columns after variant_id
            cols = base_cols[:3] + alpha_cols + base_cols[3:]
        else:
            cols = base_cols
        
        # Add n_cases and n_controls if available (for binary outcomes)
        if 'n_cases' in gwas_df.columns:
            cols = cols + ['n_cases', 'n_controls']
        
        # Select and order columns
        available_cols = [col for col in cols if col in gwas_df.columns]
        formatted_df = gwas_df[available_cols].copy()
        
        # Sort results
        if sort_by in formatted_df.columns:
            ascending = True if sort_by == 'pval' else False
            formatted_df = formatted_df.sort_values(sort_by, ascending=ascending)
        
        return formatted_df


# Public API alias (exported as format_gwas_output in __init__.py)
format_gwas_output = format_gwas_output_for_locuszoom


def save_for_locuszoom(
    gwas_df: pd.DataFrame,
    output_file: str,
    include_alpha: bool = True,
    compress: bool = True
) -> None:
    """
    Save GWAS results in LocusZoom-compatible format.
    
    Args:
        gwas_df: GWAS results DataFrame
        output_file: Output file path (will add .gz if compress=True)
        include_alpha: Include alpha_value column (default: True)
        compress: If True, compress with gzip (default: True)
        
    Returns:
        None
        
    Note:
        The output file will be tab-delimited and sorted by chromosome and position.
        For best performance with LocusZoom:
        1. Compress with bgzip: bgzip output_file.tsv
        2. Index with tabix: tabix -s 1 -b 2 -e 2 output_file.tsv.gz
        
    Example:
        >>> save_for_locuszoom(gwas_df, 'myresults.tsv', compress=True)
        >>> # Then run: bgzip -f myresults.tsv.gz && tabix -s 1 -b 2 -e 2 myresults.tsv.gz
    """
    # Format for LocusZoom
    formatted_df = format_gwas_output(
        gwas_df, 
        include_alpha=include_alpha,
        format_for_locuszoom=True
    )
    
    # Determine output filename
    if compress and not output_file.endswith('.gz'):
        final_output = output_file + '.gz'
    else:
        final_output = output_file
    
    # Save file
    if compress:
        formatted_df.to_csv(final_output, sep='\t', index=False, compression='gzip')
        print(f"Saved LocusZoom-formatted results to: {final_output}")
        print(f"\nFor optimal LocusZoom performance, run:")
        print(f"  gunzip {final_output}")
        print(f"  bgzip {final_output.replace('.gz', '')}")
        print(f"  tabix -s 1 -b 2 -e 2 {final_output}")
    else:
        formatted_df.to_csv(final_output, sep='\t', index=False)
        print(f"Saved LocusZoom-formatted results to: {final_output}")
        print(f"\nFor optimal LocusZoom performance, run:")
        print(f"  bgzip {final_output}")
        print(f"  tabix -s 1 -b 2 -e 2 {final_output}.gz")


def validate_locuszoom_format(gwas_df: pd.DataFrame) -> dict:
    """
    Validate that GWAS results meet LocusZoom format requirements.
    
    Args:
        gwas_df: GWAS results DataFrame
        
    Returns:
        Dictionary with validation results
    """
    validation = {
        'valid': True,
        'errors': [],
        'warnings': [],
        'info': []
    }
    
    # Check required columns
    required_cols = ['chrom', 'pos', 'pval']
    for col in required_cols:
        if col not in gwas_df.columns:
            validation['valid'] = False
            validation['errors'].append(f"Missing required column: {col}")
    
    # Check recommended columns
    recommended_cols = ['ref_allele', 'alt_allele', 'coef', 'std_err', 'eaf']
    for col in recommended_cols:
        if col not in gwas_df.columns:
            validation['warnings'].append(f"Missing recommended column: {col}")
    
    # Check chromosome format
    if 'chrom' in gwas_df.columns:
        valid_chroms = set([str(i) for i in range(1, 26)] + ['X', 'Y', 'M', 'MT'])
        unique_chroms = set(gwas_df['chrom'].astype(str).unique())
        invalid_chroms = unique_chroms - valid_chroms
        if invalid_chroms:
            validation['warnings'].append(
                f"Non-standard chromosome names found: {invalid_chroms}. "
                f"LocusZoom expects: 1-25, X, Y, M, or MT"
            )
    
    # Check if sorted
    if 'chrom' in gwas_df.columns and 'pos' in gwas_df.columns:
        df_check = gwas_df.copy()
        df_check['chrom_num'] = df_check['chrom'].apply(
            lambda x: 23 if str(x).upper() == 'X' else 
                     24 if str(x).upper() == 'Y' else 
                     25 if str(x).upper() in ['M', 'MT'] else 
                     int(x) if str(x).isdigit() else 99
        )
        is_sorted = (df_check[['chrom_num', 'pos']] == 
                    df_check[['chrom_num', 'pos']].sort_values(['chrom_num', 'pos'])).all().all()
        
        if not is_sorted:
            validation['warnings'].append(
                "Data not sorted by chromosome and position. "
                "LocusZoom requires sorted data."
            )
    
    # Summary
    validation['info'].append(f"Total variants: {len(gwas_df)}")
    if 'pval' in gwas_df.columns:
        sig_variants = (gwas_df['pval'] < 5e-8).sum()
        validation['info'].append(f"Genome-wide significant variants (p<5e-8): {sig_variants}")
    
    return validation


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
                vid_val = row.get('variant_id', idx)
                coef_val = row.get('coef', np.nan)
                report_lines.append(f"    {vid_val}: p={row['pval']:.2e}, coef={coef_val:.4f}")
        
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
            if len(valid_alpha) > 0:
                report_lines.append(f"  Mean alpha: {valid_alpha.mean():.4f}")
                report_lines.append(f"  Median alpha: {valid_alpha.median():.4f}")
                report_lines.append(f"  Alpha range: [{valid_alpha.min():.4f}, {valid_alpha.max():.4f}]")
            
            # Distribution of alpha values
            alpha_ranges = {
                'Near additive (0.4-0.6)': ((valid_alpha >= 0.4) & (valid_alpha <= 0.6)).sum(),
                'Dominant-like (> 0.6)': (valid_alpha > 0.6).sum(),
                'Recessive-like (< 0.4)': (valid_alpha < 0.4).sum()
            }
            
            if len(valid_alpha) > 0:
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
