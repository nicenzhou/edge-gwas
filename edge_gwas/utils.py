"""
Utility functions for EDGE GWAS analysis.
"""

import pandas as pd
import numpy as np
from pandas_plink import read_plink1_bin
from sklearn.model_selection import train_test_split
from typing import Tuple, List, Optional
import logging

logger = logging.getLogger(__name__)


def load_plink_data(
    bed_file: str,
    bim_file: str,
    fam_file: str,
    verbose: bool = True
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load PLINK binary format data.
    
    Args:
        bed_file: Path to .bed file
        bim_file: Path to .bim file
        fam_file: Path to .fam file
        verbose: Print loading information
        
    Returns:
        Tuple of (genotype_df, variant_info_df)
        - genotype_df: DataFrame with samples as index, variants as columns
        - variant_info_df: DataFrame with variant information
    """
    if verbose:
        logger.info(f"Loading PLINK data from {bed_file}")
    
    # Read PLINK data
    G = read_plink1_bin(bed_file, bim_file, fam_file, verbose=verbose)
    
    # Extract genotype matrix
    genotypes = G.values.T  # Transpose to get samples x variants
    sample_ids = G.coords['sample'].values
    variant_ids = G.coords['snp'].values
    
    # Create genotype DataFrame
    genotype_df = pd.DataFrame(
        genotypes,
        index=sample_ids,
        columns=variant_ids
    )
    genotype_df.index.name = 'sample_id'
    
    # Create variant info DataFrame
    variant_info_df = pd.DataFrame({
        'variant_id': variant_ids,
        'chrom': G.coords['chrom'].values,
        'pos': G.coords['pos'].values,
        'ref_allele': G.a1.values,
        'alt_allele': G.a0.values,
    })
    variant_info_df.set_index('variant_id', inplace=True)
    
    if verbose:
        logger.info(f"Loaded {len(sample_ids)} samples and {len(variant_ids)} variants")
    
    return genotype_df, variant_info_df


def prepare_phenotype_data(
    phenotype_file: str,
    outcome_col: str,
    covariate_cols: List[str],
    sample_id_col: str = 'IID',
    sep: str = '\t',
    log_transform_outcome: bool = False
) -> pd.DataFrame:
    """
    Load and prepare phenotype data.
    
    Args:
        phenotype_file: Path to phenotype file
        outcome_col: Name of outcome column
        covariate_cols: List of covariate column names
        sample_id_col: Name of sample ID column
        sep: File separator
        log_transform_outcome: Apply log10(x+1) transformation to outcome
        
    Returns:
        DataFrame with sample IDs as index, outcome and covariates as columns
    """
    logger.info(f"Loading phenotype data from {phenotype_file}")
    
    # Read phenotype file
    pheno_df = pd.read_csv(phenotype_file, sep=sep)
    
    # Check required columns exist
    required_cols = [sample_id_col, outcome_col] + covariate_cols
    missing_cols = [col for col in required_cols if col not in pheno_df.columns]
    if missing_cols:
        raise ValueError(f"Missing columns in phenotype file: {missing_cols}")
    
    # Set sample ID as index
    pheno_df.set_index(sample_id_col, inplace=True)
    
    # Select outcome and covariates
    pheno_df = pheno_df[[outcome_col] + covariate_cols]
    
    # Log transform outcome if requested
    if log_transform_outcome:
        logger.info(f"Applying log10(x+1) transformation to {outcome_col}")
        pheno_df[f'log10_{outcome_col}'] = np.log10(pheno_df[outcome_col] + 1)
        outcome_col = f'log10_{outcome_col}'
    
    # Remove missing values
    n_before = len(pheno_df)
    pheno_df = pheno_df.dropna()
    n_after = len(pheno_df)
    
    if n_before > n_after:
        logger.info(f"Removed {n_before - n_after} samples with missing data")
    
    logger.info(f"Prepared phenotype data for {n_after} samples")
    
    return pheno_df


def stratified_train_test_split(
    genotype_df: pd.DataFrame,
    phenotype_df: pd.DataFrame,
    outcome_col: str,
    test_size: float = 0.5,
    random_state: int = 42,
    is_binary: bool = True
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Split data into training and test sets with stratification.
    
    Args:
        genotype_df: Genotype DataFrame
        phenotype_df: Phenotype DataFrame
        outcome_col: Name of outcome column for stratification
        test_size: Proportion of data for test set
        random_state: Random seed for reproducibility
        is_binary: Whether outcome is binary (for stratification)
        
    Returns:
        Tuple of (train_geno, test_geno, train_pheno, test_pheno)
    """
    logger.info(f"Splitting data into train/test ({1-test_size:.0%}/{test_size:.0%})")
    
    # Get common samples
    common_samples = genotype_df.index.intersection(phenotype_df.index)
    
    if len(common_samples) == 0:
        raise ValueError("No common samples found between genotype and phenotype data")
    
    genotype_df = genotype_df.loc[common_samples]
    phenotype_df = phenotype_df.loc[common_samples]
    
    # Stratify if binary outcome
    if is_binary:
        stratify = phenotype_df[outcome_col]
    else:
        stratify = None
    
    # Split samples
    train_idx, test_idx = train_test_split(
        common_samples,
        test_size=test_size,
        random_state=random_state,
        stratify=stratify
    )
    
    # Split data
    train_geno = genotype_df.loc[train_idx]
    test_geno = genotype_df.loc[test_idx]
    train_pheno = phenotype_df.loc[train_idx]
    test_pheno = phenotype_df.loc[test_idx]
    
    logger.info(f"Training set: {len(train_idx)} samples")
    logger.info(f"Test set: {len(test_idx)} samples")
    
    if is_binary and stratify is not None:
        train_cases = train_pheno[outcome_col].sum()
        test_cases = test_pheno[outcome_col].sum()
        logger.info(f"Training cases/controls: {train_cases}/{len(train_idx)-train_cases}")
        logger.info(f"Test cases/controls: {test_cases}/{len(test_idx)-test_cases}")
    
    return train_geno, test_geno, train_pheno, test_pheno


def filter_variants_by_maf(
    genotype_df: pd.DataFrame,
    min_maf: float = 0.01,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Filter variants by minor allele frequency.
    
    Args:
        genotype_df: Genotype DataFrame
        min_maf: Minimum minor allele frequency
        verbose: Print filtering information
        
    Returns:
        Filtered genotype DataFrame
    """
    n_before = genotype_df.shape[1]
    
    # Calculate MAF for each variant
    mafs = []
    for col in genotype_df.columns:
        geno = genotype_df[col].dropna()
        if len(geno) == 0:
            mafs.append(0)
            continue
        
        # Calculate allele frequency
        af = (2 * (geno == 0).sum() + (geno == 1).sum()) / (2 * len(geno))
        # MAF is minimum of af and 1-af
        maf = min(af, 1 - af)
        mafs.append(maf)
    
    # Filter variants
    maf_series = pd.Series(mafs, index=genotype_df.columns)
    variants_to_keep = maf_series[maf_series >= min_maf].index
    
    genotype_df = genotype_df[variants_to_keep]
    
    n_after = genotype_df.shape[1]
    
    if verbose:
        logger.info(f"Filtered variants by MAF >= {min_maf}")
        logger.info(f"Kept {n_after}/{n_before} variants ({n_after/n_before*100:.1f}%)")
    
    return genotype_df


def filter_variants_by_missing(
    genotype_df: pd.DataFrame,
    max_missing: float = 0.1,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Filter variants by missing genotype rate.
    
    Args:
        genotype_df: Genotype DataFrame
        max_missing: Maximum proportion of missing genotypes allowed
        verbose: Print filtering information
        
    Returns:
        Filtered genotype DataFrame
    """
    n_before = genotype_df.shape[1]
    
    # Calculate missing rate for each variant
    missing_rates = genotype_df.isna().mean()
    
    # Filter variants
    variants_to_keep = missing_rates[missing_rates <= max_missing].index
    genotype_df = genotype_df[variants_to_keep]
    
    n_after = genotype_df.shape[1]
    
    if verbose:
        logger.info(f"Filtered variants by missing rate <= {max_missing}")
        logger.info(f"Kept {n_after}/{n_before} variants ({n_after/n_before*100:.1f}%)")
    
    return genotype_df


def merge_alpha_with_gwas(
    gwas_df: pd.DataFrame,
    alpha_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Merge GWAS results with alpha values and additional variant information.
    
    Args:
        gwas_df: GWAS results DataFrame
        alpha_df: Alpha values DataFrame
        
    Returns:
        Merged DataFrame
    """
    # Select relevant columns from alpha_df
    alpha_cols = ['variant_id', 'alpha_value', 'ref_allele', 'alt_allele', 
                  'eaf', 'coef_het', 'coef_hom', 'pval_het', 'pval_hom']
    
    # Only keep columns that exist
    alpha_cols = [col for col in alpha_cols if col in alpha_df.columns]
    
    merged_df = pd.merge(
        gwas_df,
        alpha_df[alpha_cols],
        on='variant_id',
        how='left',
        suffixes=('', '_alpha')
    )
    
    return merged_df


def calculate_genomic_inflation(pvalues: pd.Series) -> float:
    """
    Calculate genomic inflation factor (lambda).
    
    Args:
        pvalues: Series of p-values
        
    Returns:
        Genomic inflation factor
    """
    from scipy import stats
    
    # Remove NaN values
    pvalues = pvalues.dropna()
    
    if len(pvalues) == 0:
        return np.nan
    
    # Convert p-values to chi-square statistics (df=1)
    chi2_stats = stats.chi2.ppf(1 - pvalues, df=1)
    
    # Calculate lambda as median(chi2) / expected_median
    # Expected median for chi2(1) is qchisq(0.5, 1) = 0.4549364
    lambda_gc = np.median(chi2_stats) / 0.4549364
    
    return lambda_gc


def qq_plot_data(pvalues: pd.Series) -> Tuple[np.ndarray, np.ndarray]:
    """
    Prepare data for QQ plot.
    
    Args:
        pvalues: Series of p-values
        
    Returns:
        Tuple of (expected -log10 p-values, observed -log10 p-values)
    """
    # Remove NaN and zeros
    pvalues = pvalues.dropna()
    pvalues = pvalues[pvalues > 0]
    
    # Sort p-values
    pvalues_sorted = np.sort(pvalues)
    
    # Calculate expected p-values under null
    n = len(pvalues_sorted)
    expected_pvals = np.arange(1, n + 1) / (n + 1)
    
    # Convert to -log10 scale
    observed = -np.log10(pvalues_sorted)
    expected = -np.log10(expected_pvals)
    
    return expected, observed
