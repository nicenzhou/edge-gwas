"""
Utility functions for EDGE GWAS analysis.
"""

import os
import sys
import tempfile
import subprocess
import logging
import shlex
from typing import Tuple, Optional, Union, List, TypedDict
import pandas as pd
import numpy as np
from pandas_plink import read_plink1_bin
from scipy import stats
from sklearn.model_selection import train_test_split, KFold

logger = logging.getLogger(__name__)

# Constants
MIN_SAMPLES_FOR_ANALYSIS = 10
MIN_SAMPLES_FOR_HWE = 10
MIN_EXPECTED_COUNT_HWE = 5
LARGE_DATASET_THRESHOLD = 1e9  # 1 billion elements


# Type hints for return values
class CaseControlBalance(TypedDict):
    """Return type for case/control balance check."""
    case_count: int
    control_count: int
    ratio: float


def validate_genotype_df(genotype_df: pd.DataFrame, name: str = "genotype_df") -> None:
    """
    Validate genotype DataFrame format.
    
    Args:
        genotype_df: Genotype DataFrame to validate
        name: Name of the DataFrame for error messages
        
    Raises:
        TypeError: If not a pandas DataFrame
        ValueError: If DataFrame is empty or has duplicate sample IDs
    """
    if not isinstance(genotype_df, pd.DataFrame):
        raise TypeError(f"{name} must be a pandas DataFrame")
    
    if genotype_df.empty:
        raise ValueError(f"{name} is empty")
    
    if genotype_df.index.duplicated().any():
        n_duplicates = genotype_df.index.duplicated().sum()
        raise ValueError(
            f"{name} has {n_duplicates} duplicate sample IDs in index. "
            f"Example duplicates: {genotype_df.index[genotype_df.index.duplicated()].unique()[:5].tolist()}"
        )


def validate_phenotype_df(
    phenotype_df: pd.DataFrame,
    outcome_col: str,
    covariate_cols: List[str]
) -> None:
    """
    Validate phenotype DataFrame format.
    
    Args:
        phenotype_df: Phenotype DataFrame to validate
        outcome_col: Name of outcome column
        covariate_cols: List of covariate column names
        
    Raises:
        TypeError: If not a pandas DataFrame
        ValueError: If DataFrame is empty, has duplicates, or missing columns
    """
    validate_genotype_df(phenotype_df, "phenotype_df")
    
    required_cols = [outcome_col] + covariate_cols
    missing_cols = [col for col in required_cols if col not in phenotype_df.columns]
    if missing_cols:
        raise ValueError(
            f"Missing columns in phenotype_df: {missing_cols}. "
            f"Available columns: {phenotype_df.columns.tolist()}"
        )


def load_plink_data(
    bed_file: str,
    bim_file: str,
    fam_file: str,
    verbose: bool = True
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load PLINK binary format data (.bed/.bim/.fam).
    
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
    genotypes = G.values 
    sample_ids = G.coords['sample'].values
    variant_ids = G.coords['snp'].values
    
    # Convert to numpy array if it's a dask/xarray
    if hasattr(genotypes, 'compute'):
        genotypes = genotypes.compute()
    genotypes = np.asarray(genotypes)
    
    # Create DataFrame: variants as rows, samples as columns
    # Then transpose to get samples as rows, variants as columns
    genotype_df = pd.DataFrame(
        genotypes,  
        index=pd.Index(sample_ids, name='sample_id'),
        columns=pd.Index(variant_ids, name='variant_id')
    )
    
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


def load_pgen_data(
    pgen_file: str,
    pvar_file: str,
    psam_file: str,
    verbose: bool = True
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load PLINK 2 binary format data (.pgen/.pvar/.psam).
    
    Args:
        pgen_file: Path to .pgen file
        pvar_file: Path to .pvar file
        psam_file: Path to .psam file
        verbose: Print loading information
        
    Returns:
        Tuple of (genotype_df, variant_info_df)
        - genotype_df: DataFrame with samples as index, variants as columns
        - variant_info_df: DataFrame with variant information
        
    Note:
        Requires pgenlib package: pip install pgenlib
    """
    try:
        from pgenlib import PgenReader
    except ImportError:
        raise ImportError(
            "pgenlib is required to load PGEN files. "
            "Install with: pip install pgenlib"
        )
    
    if verbose:
        logger.info(f"Loading PGEN data from {pgen_file}")
    
    # Read .psam file (sample information)
    psam_df = pd.read_csv(psam_file, sep='\t', comment='#')
    sample_ids = psam_df['IID'].values if 'IID' in psam_df.columns else psam_df.iloc[:, 1].values
    
    # Read .pvar file (variant information)
    pvar_df = pd.read_csv(pvar_file, sep='\t', comment='#')
    variant_ids = pvar_df['ID'].values if 'ID' in pvar_df.columns else pvar_df.iloc[:, 2].values
    chromosomes = pvar_df['CHROM'].values if 'CHROM' in pvar_df.columns else pvar_df.iloc[:, 0].values
    positions = pvar_df['POS'].values if 'POS' in pvar_df.columns else pvar_df.iloc[:, 1].values
    ref_alleles = pvar_df['REF'].values if 'REF' in pvar_df.columns else pvar_df.iloc[:, 3].values
    alt_alleles = pvar_df['ALT'].values if 'ALT' in pvar_df.columns else pvar_df.iloc[:, 4].values
    
    # Open .pgen file
    pgen_reader = PgenReader(pgen_file.encode())
    
    # Initialize genotype matrix
    n_samples = len(sample_ids)
    n_variants = len(variant_ids)
    
    # Use memory mapping for large datasets
    if n_samples * n_variants > LARGE_DATASET_THRESHOLD:
        if verbose:
            logger.info("Large dataset detected. Using memory-mapped array.")
        temp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.dat')
        genotypes = np.memmap(
            temp_file.name,
            dtype=np.float32,
            mode='w+',
            shape=(n_samples, n_variants)
        )
    else:
        genotypes = np.zeros((n_samples, n_variants), dtype=np.float32)
    
    # Read genotypes
    variant_buffer = np.empty(n_samples, dtype=np.int32)
    for i in range(n_variants):
        pgen_reader.read(i, variant_buffer)
        genotypes[:, i] = variant_buffer
    
    pgen_reader.close()
    
    # Convert missing values (-9) to NaN
    genotypes[genotypes == -9] = np.nan
    
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
        'chrom': chromosomes,
        'pos': positions,
        'ref_allele': ref_alleles,
        'alt_allele': alt_alleles,
    })
    variant_info_df.set_index('variant_id', inplace=True)
    
    if verbose:
        logger.info(f"Loaded {n_samples} samples and {n_variants} variants")
    
    return genotype_df, variant_info_df


def load_bgen_data(
    bgen_file: str,
    sample_file: Optional[str] = None,
    verbose: bool = True
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load BGEN format data.
    
    Args:
        bgen_file: Path to .bgen file
        sample_file: Path to .sample file (optional, can be embedded in BGEN)
        verbose: Print loading information
        
    Returns:
        Tuple of (genotype_df, variant_info_df)
        - genotype_df: DataFrame with samples as index, variants as columns (dosages)
        - variant_info_df: DataFrame with variant information
        
    Note:
        Requires bgen_reader package: pip install bgen-reader
    """
    try:
        from bgen_reader import open_bgen
    except ImportError:
        raise ImportError(
            "bgen_reader is required to load BGEN files. "
            "Install with: pip install bgen-reader"
        )
    
    if verbose:
        logger.info(f"Loading BGEN data from {bgen_file}")
    
    # Open BGEN file
    bgen = open_bgen(bgen_file, verbose=verbose)
    
    # Read genotype probabilities
    # This returns probabilities for each genotype (0/0, 0/1, 1/1)
    probs = bgen.read()
    
    # Convert probabilities to dosages (expected allele count)
    # Dosage = P(0/1) * 1 + P(1/1) * 2
    dosages = probs[:, :, 1] + 2 * probs[:, :, 2]
    
    # Get sample IDs
    if sample_file is not None:
        # Read sample file (skip first two lines which are headers)
        sample_df = pd.read_csv(sample_file, sep=' ', skiprows=2, header=None)
        sample_ids = sample_df[0].values  # First column is ID_1
    else:
        # Use sample IDs from BGEN file
        sample_ids = bgen.samples
    
    # Get variant information
    variant_ids = bgen.ids
    chromosomes = bgen.chromosomes
    positions = bgen.positions
    ref_alleles = bgen.allele_ids[:, 0]
    alt_alleles = bgen.allele_ids[:, 1]
    
    # Create genotype DataFrame (samples x variants)
    genotype_df = pd.DataFrame(
        dosages.T,  # Transpose to get samples x variants
        index=sample_ids,
        columns=variant_ids
    )
    genotype_df.index.name = 'sample_id'
    
    # Create variant info DataFrame
    variant_info_df = pd.DataFrame({
        'variant_id': variant_ids,
        'chrom': chromosomes,
        'pos': positions,
        'ref_allele': ref_alleles,
        'alt_allele': alt_alleles,
    })
    variant_info_df.set_index('variant_id', inplace=True)
    
    if verbose:
        logger.info(f"Loaded {len(sample_ids)} samples and {len(variant_ids)} variants")
        logger.info("Note: BGEN genotypes are dosages (0-2 continuous values)")
    
    return genotype_df, variant_info_df


def load_vcf_data(
    vcf_file: str,
    dosage: bool = True,
    verbose: bool = True
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load VCF format data.
    
    Args:
        vcf_file: Path to .vcf or .vcf.gz file
        dosage: If True, use dosages (DS field); if False, use hard calls (GT field)
        verbose: Print loading information
        
    Returns:
        Tuple of (genotype_df, variant_info_df)
        - genotype_df: DataFrame with samples as index, variants as columns
        - variant_info_df: DataFrame with variant information
        
    Note:
        Requires cyvcf2 package: pip install cyvcf2
    """
    try:
        from cyvcf2 import VCF
    except ImportError:
        raise ImportError(
            "cyvcf2 is required to load VCF files. "
            "Install with: pip install cyvcf2"
        )
    
    if verbose:
        logger.info(f"Loading VCF data from {vcf_file}")
    
    # Open VCF file
    vcf = VCF(vcf_file)
    
    # Get sample IDs
    sample_ids = vcf.samples
    
    # Initialize lists to store data
    genotypes_list = []
    variant_ids_list = []
    chromosomes_list = []
    positions_list = []
    ref_alleles_list = []
    alt_alleles_list = []
    
    # Read variants
    for variant in vcf:
        # Get variant information
        variant_id = variant.ID if variant.ID else f"{variant.CHROM}:{variant.POS}"
        variant_ids_list.append(variant_id)
        chromosomes_list.append(variant.CHROM)
        positions_list.append(variant.POS)
        ref_alleles_list.append(variant.REF)
        alt_alleles_list.append(variant.ALT[0] if variant.ALT else '.')
        
        # Get genotypes
        if dosage and 'DS' in variant.FORMAT:
            # Use dosages if available
            geno = variant.format('DS')[:, 0]  # Get dosage for first alt allele
        else:
            # Use hard calls (convert GT to 0/1/2)
            gt = variant.genotypes
            # gt is list of [allele1, allele2, phased] for each sample
            geno = np.array([sum(g[:2]) if -1 not in g[:2] else np.nan 
                           for g in gt])
        
        genotypes_list.append(geno)
    
    # Close VCF
    vcf.close()
    
    # Convert to arrays
    genotypes = np.array(genotypes_list).T  # Transpose to samples x variants
    
    # Create genotype DataFrame
    genotype_df = pd.DataFrame(
        genotypes,
        index=sample_ids,
        columns=variant_ids_list
    )
    genotype_df.index.name = 'sample_id'
    
    # Create variant info DataFrame
    variant_info_df = pd.DataFrame({
        'variant_id': variant_ids_list,
        'chrom': chromosomes_list,
        'pos': positions_list,
        'ref_allele': ref_alleles_list,
        'alt_allele': alt_alleles_list,
    })
    variant_info_df.set_index('variant_id', inplace=True)
    
    if verbose:
        logger.info(f"Loaded {len(sample_ids)} samples and {len(variant_ids_list)} variants")
        if dosage:
            logger.info("Using dosages (DS field)")
        else:
            logger.info("Using hard calls (GT field)")
    
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
    
    # Validate required columns
    validate_phenotype_df(pheno_df, outcome_col, covariate_cols)
    
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
    is_binary: bool = True,
    geno_id_col: Optional[str] = None,
    pheno_id_col: Optional[str] = None
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Split data into training and test sets with stratification.
    
    Args:
        genotype_df: Genotype DataFrame (samples x variants)
        phenotype_df: Phenotype DataFrame
        outcome_col: Name of outcome column for stratification
        test_size: Proportion of samples in test set (default: 0.5)
        random_state: Random seed for reproducibility
        is_binary: Whether outcome is binary (enables stratification)
        geno_id_col: Column name in genotype_df for sample IDs
                    If None, uses genotype_df.index
        pheno_id_col: Column name in phenotype_df for sample IDs
                     If None, uses phenotype_df.index
        
    Returns:
        Tuple of (train_geno, test_geno, train_pheno, test_pheno)
        
    Raises:
        ValueError: If no common samples found or stratification fails
    """
    logger.info(f"Splitting data into train/test ({1-test_size:.0%}/{test_size:.0%})")
    
    # Get sample identifiers from genotype DataFrame
    if geno_id_col is None:
        geno_samples = genotype_df.index
        logger.debug("Using genotype DataFrame index for sample IDs")
    elif geno_id_col in genotype_df.columns:
        geno_samples = genotype_df[geno_id_col]
        logger.debug(f"Using genotype column '{geno_id_col}' for sample IDs")
    elif genotype_df.index.name == geno_id_col:
        geno_samples = genotype_df.index
        logger.debug(f"Using genotype index (name='{geno_id_col}') for sample IDs")
    else:
        raise ValueError(
            f"Sample ID column '{geno_id_col}' not found in genotype DataFrame. "
            f"Available columns: {genotype_df.columns.tolist()}, "
            f"Index name: '{genotype_df.index.name}'"
        )
    
    # Get sample identifiers from phenotype DataFrame
    if pheno_id_col is None:
        pheno_samples = phenotype_df.index
        logger.debug("Using phenotype DataFrame index for sample IDs")
    elif pheno_id_col in phenotype_df.columns:
        pheno_samples = phenotype_df[pheno_id_col]
        logger.debug(f"Using phenotype column '{pheno_id_col}' for sample IDs")
    elif phenotype_df.index.name == pheno_id_col:
        pheno_samples = phenotype_df.index
        logger.debug(f"Using phenotype index (name='{pheno_id_col}') for sample IDs")
    else:
        raise ValueError(
            f"Sample ID column '{pheno_id_col}' not found in phenotype DataFrame. "
            f"Available columns: {phenotype_df.columns.tolist()}, "
            f"Index name: '{phenotype_df.index.name}'"
        )
    
    # Convert both to strings for consistent comparison
    geno_samples_str = pd.Index(geno_samples.astype(str))
    pheno_samples_str = pd.Index(pheno_samples.astype(str))
    
    # Check for duplicate string conversions
    if len(geno_samples_str) != len(set(geno_samples_str)):
        n_duplicates = len(geno_samples_str) - len(set(geno_samples_str))
        logger.warning(
            f"{n_duplicates} genotype sample IDs map to the same string representation. "
            f"This may cause issues with sample matching."
        )
    if len(pheno_samples_str) != len(set(pheno_samples_str)):
        n_duplicates = len(pheno_samples_str) - len(set(pheno_samples_str))
        logger.warning(
            f"{n_duplicates} phenotype sample IDs map to the same string representation. "
            f"This may cause issues with sample matching."
        )
    
    # Find common samples (using string comparison)
    common_samples_str = geno_samples_str.intersection(pheno_samples_str)
    
    if len(common_samples_str) == 0:
        raise ValueError(
            f"No common samples found between genotype and phenotype data.\n"
            f"Genotype has {len(geno_samples)} samples (type: {type(geno_samples[0]).__name__}): "
            f"{list(geno_samples[:5])}...\n"
            f"Phenotype has {len(pheno_samples)} samples (type: {type(pheno_samples[0]).__name__}): "
            f"{list(pheno_samples[:5])}...\n"
            f"Hint: Check if sample IDs match between files. "
            f"Use pheno.set_index('IID') to set sample IDs as index."
        )
    
    # Convert back to original type for indexing
    # Create mapping from string to original values
    geno_str_to_orig = dict(zip(geno_samples_str, geno_samples))
    pheno_str_to_orig = dict(zip(pheno_samples_str, pheno_samples))
    
    # Get common samples in original types
    common_samples_geno = pd.Index([geno_str_to_orig[s] for s in common_samples_str])
    common_samples_pheno = pd.Index([pheno_str_to_orig[s] for s in common_samples_str])
    
    n_dropped_geno = len(geno_samples) - len(common_samples_str)
    n_dropped_pheno = len(pheno_samples) - len(common_samples_str)
    
    logger.info(
        f"Found {len(common_samples_str)} common samples "
        f"(dropped {n_dropped_geno} from genotype, {n_dropped_pheno} from phenotype)"
    )
    
    # Subset to common samples and align
    if geno_id_col is None or genotype_df.index.name == geno_id_col:
        genotype_df_subset = genotype_df.loc[common_samples_geno]
    else:
        genotype_df_subset = genotype_df[genotype_df[geno_id_col].isin(common_samples_geno)].copy()
        genotype_df_subset = genotype_df_subset.set_index(geno_id_col)
        genotype_df_subset = genotype_df_subset.loc[common_samples_geno]
    
    if pheno_id_col is None or phenotype_df.index.name == pheno_id_col:
        phenotype_df_subset = phenotype_df.loc[common_samples_pheno]
    else:
        phenotype_df_subset = phenotype_df[phenotype_df[pheno_id_col].isin(common_samples_pheno)].copy()
        phenotype_df_subset = phenotype_df_subset.set_index(pheno_id_col)
        phenotype_df_subset = phenotype_df_subset.loc[common_samples_pheno]
    
    # Ensure indices match (convert pheno index to match geno type)
    phenotype_df_subset.index = genotype_df_subset.index
    
    # Prepare stratification
    if is_binary:
        stratify = phenotype_df_subset[outcome_col]
        
        # Validate stratification is possible
        value_counts = stratify.value_counts()
        if len(value_counts) < 2:
            logger.warning(
                f"Outcome has only {len(value_counts)} unique value(s). "
                f"Disabling stratification."
            )
            stratify = None
        elif any(value_counts < 2):
            logger.warning(
                f"Some outcome classes have fewer than 2 samples: {value_counts.to_dict()}. "
                f"Stratification may fail."
            )
    else:
        stratify = None
    
    # Perform train-test split (use genotype index for consistency)
    try:
        train_idx, test_idx = train_test_split(
            genotype_df_subset.index,
            test_size=test_size,
            random_state=random_state,
            stratify=stratify
        )
    except ValueError as e:
        if stratify is not None:
            logger.warning(
                f"Stratification failed: {e}. "
                f"Falling back to random split without stratification."
            )
            train_idx, test_idx = train_test_split(
                genotype_df_subset.index,
                test_size=test_size,
                random_state=random_state,
                stratify=None
            )
        else:
            raise
    
    # Split data
    train_geno = genotype_df_subset.loc[train_idx]
    test_geno = genotype_df_subset.loc[test_idx]
    train_pheno = phenotype_df_subset.loc[train_idx]
    test_pheno = phenotype_df_subset.loc[test_idx]
    
    logger.info(f"Training set: {len(train_idx)} samples")
    logger.info(f"Test set: {len(test_idx)} samples")
    
    # Log class distribution for binary outcomes
    if is_binary and stratify is not None:
        train_outcome = train_pheno[outcome_col]
        test_outcome = test_pheno[outcome_col]
        
        # Handle both binary (0/1) and case/control (1/2) encoding
        train_cases = int(train_outcome.sum())
        test_cases = int(test_outcome.sum())
        
        logger.info(
            f"Training cases/controls: {train_cases}/{len(train_idx)-train_cases} "
            f"({train_cases/len(train_idx):.1%} cases)"
        )
        logger.info(
            f"Test cases/controls: {test_cases}/{len(test_idx)-test_cases} "
            f"({test_cases/len(test_idx):.1%} cases)"
        )
    
    return train_geno, test_geno, train_pheno, test_pheno


def filter_variants_by_maf(
    genotype_df: pd.DataFrame,
    min_maf: float = 0.01,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Filter variants by minor allele frequency (optimized vectorized version).
    
    Args:
        genotype_df: Genotype DataFrame (works with both hard calls and dosages)
        min_maf: Minimum minor allele frequency
        verbose: Print filtering information
        
    Returns:
        Filtered genotype DataFrame
    """
    n_before = genotype_df.shape[1]
    
    # Vectorized MAF calculation
    af = genotype_df.mean(axis=0) / 2  # Allele frequency
    
    # Handle NaN values (variants with all missing data)
    af = af.fillna(0)
    
    # MAF is min(af, 1-af)
    maf = np.minimum(af, 1 - af)
    
    # Filter variants
    variants_to_keep = maf[maf >= min_maf].index
    genotype_df_filtered = genotype_df[variants_to_keep]
    
    n_after = genotype_df_filtered.shape[1]
    
    if verbose:
        logger.info(f"Filtered variants by MAF >= {min_maf}")
        logger.info(f"Kept {n_after}/{n_before} variants ({n_after/n_before*100:.1f}%)")
    
    return genotype_df_filtered


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
    missing_rates = genotype_df.isna().mean(axis=0)
    
    # Filter variants
    variants_to_keep = missing_rates[missing_rates <= max_missing].index
    genotype_df_filtered = genotype_df[variants_to_keep]
    
    n_after = genotype_df_filtered.shape[1]
    
    if verbose:
        logger.info(f"Filtered variants by missing rate <= {max_missing}")
        logger.info(f"Kept {n_after}/{n_before} variants ({n_after/n_before*100:.1f}%)")
    
    return genotype_df_filtered


def filter_samples_by_call_rate(
    genotype_df: pd.DataFrame,
    phenotype_df: pd.DataFrame,
    min_call_rate: float = 0.95,
    verbose: bool = True
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Filter samples by genotype call rate.
    
    Args:
        genotype_df: Genotype DataFrame
        phenotype_df: Phenotype DataFrame
        min_call_rate: Minimum call rate (proportion of non-missing genotypes)
        verbose: Print filtering information
        
    Returns:
        Tuple of (filtered_genotype_df, filtered_phenotype_df)
    """
    n_before = genotype_df.shape[0]
    
    # Calculate call rate for each sample
    sample_call_rate = genotype_df.notna().mean(axis=1)
    
    # Filter samples
    good_samples = sample_call_rate[sample_call_rate >= min_call_rate].index
    
    genotype_df_filtered = genotype_df.loc[good_samples]
    phenotype_df_filtered = phenotype_df.loc[phenotype_df.index.intersection(good_samples)]
    
    n_after = genotype_df_filtered.shape[0]
    
    if verbose:
        logger.info(f"Filtered samples by call rate >= {min_call_rate}")
        logger.info(f"Kept {n_after}/{n_before} samples ({n_after/n_before*100:.1f}%)")
    
    return genotype_df_filtered, phenotype_df_filtered


def check_case_control_balance(
    phenotype_df: pd.DataFrame,
    outcome_col: str,
    verbose: bool = True
) -> CaseControlBalance:
    """
    Check case/control balance in binary outcome.
    
    Args:
        phenotype_df: Phenotype DataFrame
        outcome_col: Name of outcome column
        verbose: Print balance information
        
    Returns:
        Dictionary with case_count, control_count, and ratio
    """
    case_count = int(phenotype_df[outcome_col].sum())
    control_count = len(phenotype_df) - case_count
    
    if control_count > 0:
        ratio = case_count / control_count
    else:
        ratio = np.inf
    
    if verbose:
        logger.info(f"Cases: {case_count}, Controls: {control_count}")
        logger.info(f"Case/control ratio: {ratio:.2f}")
    
    return {
        'case_count': case_count,
        'control_count': control_count,
        'ratio': ratio
    }


def calculate_hwe_pvalues(
    genotype_df: pd.DataFrame,
    verbose: bool = True
) -> pd.Series:
    """
    Calculate Hardy-Weinberg Equilibrium p-values for each variant.
    
    Args:
        genotype_df: Genotype DataFrame
        verbose: Print calculation information
        
    Returns:
        Series of HWE p-values for each variant
    """
    try:
        from tqdm import tqdm
        use_tqdm = True
    except ImportError:
        use_tqdm = False
    
    if verbose:
        logger.info("Calculating Hardy-Weinberg Equilibrium p-values...")
    
    hwe_pvals = {}
    
    variants = genotype_df.columns
    iterator = tqdm(variants, desc="Calculating HWE") if (verbose and use_tqdm) else variants
    
    for variant in iterator:
        geno = genotype_df[variant].dropna()
        
        # Only use hard calls (0, 1, 2)
        geno = geno[geno.isin([0, 1, 2])]
        
        if len(geno) < MIN_SAMPLES_FOR_HWE:  # Skip if too few samples
            hwe_pvals[variant] = np.nan
            continue
        
        # Count genotypes
        n_aa = (geno == 0).sum()  # Homozygous reference
        n_ab = (geno == 1).sum()  # Heterozygous
        n_bb = (geno == 2).sum()  # Homozygous alternate
        
        n_total = n_aa + n_ab + n_bb
        
        # Calculate allele frequencies
        p = (2 * n_aa + n_ab) / (2 * n_total)  # Reference allele frequency
        q = 1 - p  # Alternate allele frequency
        
        # Expected genotype counts under HWE
        exp_aa = n_total * p * p
        exp_ab = n_total * 2 * p * q
        exp_bb = n_total * q * q
        
        # Chi-square test (only if expected counts are sufficient)
        if exp_aa >= MIN_EXPECTED_COUNT_HWE and exp_ab >= MIN_EXPECTED_COUNT_HWE and exp_bb >= MIN_EXPECTED_COUNT_HWE:
            obs = np.array([n_aa, n_ab, n_bb])
            exp = np.array([exp_aa, exp_ab, exp_bb])
            
            chi2_stat = np.sum((obs - exp) ** 2 / exp)
            pval = 1 - stats.chi2.cdf(chi2_stat, df=1)
            
            hwe_pvals[variant] = pval
        else:
            # Too few expected counts for valid chi-square test
            hwe_pvals[variant] = np.nan
    
    hwe_series = pd.Series(hwe_pvals)
    
    if verbose:
        n_valid = hwe_series.notna().sum()
        n_total = len(hwe_series)
        logger.info(f"Calculated HWE p-values for {n_valid}/{n_total} variants")
        if n_valid < n_total:
            logger.info(f"  {n_total - n_valid} variants skipped (insufficient data)")
    
    return hwe_series


def filter_variants_by_hwe(
    genotype_df: pd.DataFrame,
    hwe_threshold: float = 1e-6,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Filter variants by Hardy-Weinberg Equilibrium p-value.
    
    Args:
        genotype_df: Genotype DataFrame
        hwe_threshold: Minimum HWE p-value threshold
        verbose: Print filtering information
        
    Returns:
        Filtered genotype DataFrame
    """
    n_before = genotype_df.shape[1]
    
    # Calculate HWE p-values
    hwe_pvals = calculate_hwe_pvalues(genotype_df, verbose=False)
    
    # Filter variants (keep those with p-value >= threshold or NaN)
    variants_to_keep = hwe_pvals[(hwe_pvals >= hwe_threshold) | hwe_pvals.isna()].index
    genotype_df_filtered = genotype_df[variants_to_keep]
    
    n_after = genotype_df_filtered.shape[1]
    
    if verbose:
        n_filtered_hwe = (hwe_pvals < hwe_threshold).sum()
        n_na = hwe_pvals.isna().sum()
        logger.info(f"Filtered variants by HWE p-value >= {hwe_threshold}")
        logger.info(f"  Removed {n_filtered_hwe} variants failing HWE")
        logger.info(f"  Kept {n_na} variants with insufficient data for HWE test")
        logger.info(f"Kept {n_after}/{n_before} variants ({n_after/n_before*100:.1f}%)")
    
    return genotype_df_filtered


def standard_gwas(
    genotype_df: pd.DataFrame,
    phenotype_df: pd.DataFrame,
    outcome: str,
    covariates: List[str],
    outcome_type: str = 'binary'
) -> pd.DataFrame:
    """
    Perform standard additive GWAS for comparison with EDGE.
    
    Args:
        genotype_df: Genotype DataFrame
        phenotype_df: Phenotype DataFrame
        outcome: Name of outcome column
        covariates: List of covariate column names
        outcome_type: 'binary' for logistic regression, 'continuous' for linear regression
        
    Returns:
        DataFrame with variant_id, coef, pval, std_err
    """
    from sklearn.linear_model import LogisticRegression, LinearRegression
    
    logger.info(f"Running standard GWAS ({outcome_type} outcome)...")
    
    results = []
    
    try:
        from tqdm import tqdm
        use_tqdm = True
    except ImportError:
        use_tqdm = False
    
    variants = genotype_df.columns
    iterator = tqdm(variants, desc="Standard GWAS") if use_tqdm else variants
    
    for variant in iterator:
        try:
            # Prepare data
            X = phenotype_df[covariates].copy()
            X['genotype'] = genotype_df[variant]
            X = X.dropna()
            
            if len(X) < MIN_SAMPLES_FOR_ANALYSIS:  # Skip if too few samples
                continue
            
            y = phenotype_df.loc[X.index, outcome]
            
            # Fit model
            if outcome_type == 'binary':
                model = LogisticRegression(max_iter=1000, penalty=None)
                model.fit(X, y)
                
                # Get coefficient for genotype (last column)
                coef = model.coef_[0][-1]
                
                # Calculate standard error using Hessian approximation
                # This is simplified - full calculation would use statsmodels
                predictions = model.predict_proba(X)[:, 1]
                W = predictions * (1 - predictions)
                
                # Wald test approximation
                std_err = np.sqrt(1 / (np.sum(W * X['genotype']**2) + 1e-10))
                z_score = coef / std_err
                pval = 2 * (1 - stats.norm.cdf(abs(z_score)))
                
            else:  # continuous
                model = LinearRegression()
                model.fit(X, y)
                
                # Get coefficient
                coef = model.coef_[-1]
                
                # Calculate standard error
                y_pred = model.predict(X)
                residuals = y - y_pred
                mse = np.sum(residuals**2) / (len(X) - X.shape[1])
                
                # Standard error for last coefficient
                X_matrix = X.values
                var_coef = mse * np.linalg.inv(X_matrix.T @ X_matrix)[-1, -1]
                std_err = np.sqrt(var_coef)
                
                # T-test
                t_stat = coef / std_err
                pval = 2 * (1 - stats.t.cdf(abs(t_stat), df=len(X) - X.shape[1]))
            
            results.append({
                'variant_id': variant,
                'coef': coef,
                'std_err': std_err,
                'pval': pval
            })
            
        except Exception as e:
            logger.warning(f"Failed to analyze variant {variant}: {str(e)}")
            continue
    
    results_df = pd.DataFrame(results)
    
    logger.info(f"Completed standard GWAS for {len(results_df)} variants")
    
    return results_df


def cross_validated_edge_analysis(
    genotype_df: pd.DataFrame,
    phenotype_df: pd.DataFrame,
    outcome: str,
    covariates: List[str],
    outcome_type: str = 'binary',
    n_folds: int = 5,
    n_jobs: int = 8,
    random_state: int = 42
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Perform k-fold cross-validation for EDGE analysis.
    
    Args:
        genotype_df: Genotype DataFrame
        phenotype_df: Phenotype DataFrame
        outcome: Name of outcome column
        covariates: List of covariate column names
        outcome_type: 'binary' or 'continuous'
        n_folds: Number of cross-validation folds
        n_jobs: Number of parallel jobs for EDGE analysis
        random_state: Random seed for reproducibility
        
    Returns:
        Tuple of (avg_alpha, meta_gwas_df, combined_alpha, combined_gwas)
        - avg_alpha: Averaged alpha values across folds
        - meta_gwas_df: Meta-analyzed GWAS results
        - combined_alpha: All alpha values from all folds
        - combined_gwas: All GWAS results from all folds
    """
    from edge_gwas.core import EDGEAnalysis
    from scipy.stats import combine_pvalues
    
    kf = KFold(n_splits=n_folds, shuffle=True, random_state=random_state)
    edge = EDGEAnalysis(outcome_type=outcome_type, n_jobs=n_jobs)
    
    all_alpha_values = []
    all_gwas_results = []
    
    for fold, (train_idx, test_idx) in enumerate(kf.split(genotype_df), 1):
        logger.info(f"Processing fold {fold}/{n_folds}...")
        
        # Split data
        train_samples = genotype_df.index[train_idx]
        test_samples = genotype_df.index[test_idx]
        
        train_g = genotype_df.loc[train_samples]
        test_g = genotype_df.loc[test_samples]
        train_p = phenotype_df.loc[train_samples]
        test_p = phenotype_df.loc[test_samples]
        
        # Run EDGE
        alpha_df, gwas_df = edge.run_full_analysis(
            train_g, train_p, test_g, test_p,
            outcome=outcome,
            covariates=covariates
        )
        
        # Store results
        alpha_df['fold'] = fold
        gwas_df['fold'] = fold
        all_alpha_values.append(alpha_df)
        all_gwas_results.append(gwas_df)
    
    # Combine results
    combined_alpha = pd.concat(all_alpha_values, ignore_index=True)
    combined_gwas = pd.concat(all_gwas_results, ignore_index=True)
    
    # Average alpha values across folds
    avg_alpha = combined_alpha.groupby('variant_id').agg({
        'alpha_value': ['mean', 'std'],
        'eaf': 'mean'
    }).reset_index()
    avg_alpha.columns = ['variant_id', 'alpha_mean', 'alpha_std', 'eaf']
    
    # Meta-analysis of p-values across folds (Fisher's method)
    meta_gwas = []
    for variant in combined_gwas['variant_id'].unique():
        variant_data = combined_gwas[combined_gwas['variant_id'] == variant]
        
        # Combine p-values using Fisher's method
        _, combined_pval = combine_pvalues(variant_data['pval'], method='fisher')
        
        meta_gwas.append({
            'variant_id': variant,
            'pval': combined_pval,
            'mean_coef': variant_data['coef'].mean(),
            'std_coef': variant_data['coef'].std()
        })
    
    meta_gwas_df = pd.DataFrame(meta_gwas)
    
    logger.info("Cross-validation complete")
    logger.info(f"Mean alpha std across variants: {avg_alpha['alpha_std'].mean():.3f}")
    
    return avg_alpha, meta_gwas_df, combined_alpha, combined_gwas


def calculate_pca_sklearn(
    genotype_df: pd.DataFrame,
    n_pcs: int = 10,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Calculate principal components using scikit-learn (basic PCA without relatedness correction).
    
    Args:
        genotype_df: Genotype DataFrame (samples x variants)
        n_pcs: Number of principal components to calculate
        verbose: Print progress information
        
    Returns:
        DataFrame with 'IID' column and PC1, PC2, ..., PCn columns
        Index is set to IID
        
    Note:
        This is a basic PCA without correction for relatedness.
        For more robust PCA accounting for relatedness, use calculate_pca_plink().
    """
    from sklearn.decomposition import PCA
    from sklearn.impute import SimpleImputer
    
    if verbose:
        logger.info(f"Calculating {n_pcs} principal components using scikit-learn...")
    
    # Validate input
    validate_genotype_df(genotype_df)
    
    # Impute missing values with mean
    imputer = SimpleImputer(strategy='mean')
    genotype_imputed = imputer.fit_transform(genotype_df.values)
    
    # Standardize genotypes (mean=0, std=1)
    genotype_std = (genotype_imputed - genotype_imputed.mean(axis=0)) / (genotype_imputed.std(axis=0) + 1e-10)
    
    # Calculate PCA
    pca = PCA(n_components=n_pcs)
    pcs = pca.fit_transform(genotype_std)
    
    # Create DataFrame with IID column
    pc_cols = [f'PC{i+1}' for i in range(n_pcs)]
    pca_df = pd.DataFrame(
        pcs,
        columns=pc_cols
    )
    pca_df['IID'] = genotype_df.index.astype(str)
    
    # Reorder columns to put IID first
    pca_df = pca_df[['IID'] + pc_cols]
    
    # Set IID as index
    pca_df.set_index('IID', inplace=True)
    
    if verbose:
        logger.info(f"Explained variance ratio: {pca.explained_variance_ratio_[:5]}")
        logger.info(f"Total variance explained by {n_pcs} PCs: {pca.explained_variance_ratio_.sum():.3f}")
        logger.info(f"PCA complete for {len(pca_df)} samples")
    
    return pca_df


def calculate_pca_plink(
    plink_prefix: str,
    n_pcs: int = 10,
    output_prefix: Optional[str] = None,
    maf_threshold: Optional[float] = 0.01,
    ld_window: Optional[int] = 50,
    ld_step: Optional[int] = 5,
    ld_r2: Optional[float] = 0.2,
    approx: bool = False,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Calculate principal components using PLINK2.
    
    Args:
        plink_prefix: Prefix for PLINK binary files (.bed/.bim/.fam)
        n_pcs: Number of principal components to calculate
        output_prefix: Prefix for output files (default: temp directory)
        maf_threshold: MAF threshold for variant filtering (default: 0.01, None to skip)
        ld_window: Window size for LD pruning in variant count (default: 50, None to skip LD pruning)
        ld_step: Step size for LD pruning in variant count (default: 5, None to skip LD pruning)
        ld_r2: R² threshold for LD pruning (default: 0.2, None to skip LD pruning)
        approx: Use approximate PCA for large cohorts (faster, recommended for >5000 samples)
        verbose: Print progress information
        
    Returns:
        DataFrame with 'IID' column and PC1, PC2, ..., PCn columns
        Index is set to IID
        
    Note:
        Requires PLINK2 to be installed and available in PATH.
        Download from: https://www.cog-genomics.org/plink/2.0/
        
        For large cohorts (>5000 samples), use approx=True for faster computation.
        Set maf_threshold, ld_window, ld_step, or ld_r2 to None to skip those filters.
    """
    # Create temporary directory if no output prefix specified
    if output_prefix is None:
        temp_dir = tempfile.mkdtemp()
        output_prefix = os.path.join(temp_dir, 'pca')
        cleanup = True
    else:
        temp_dir = None
        cleanup = False
    
    try:
        if verbose:
            method = "approximate" if approx else "exact"
            logger.info(f"Calculating {n_pcs} PCs using PLINK2 ({method} method)...")
            if maf_threshold is not None:
                logger.info(f"MAF threshold: {maf_threshold}")
            if ld_window is not None and ld_r2 is not None:
                logger.info(f"LD pruning: window={ld_window}, step={ld_step}, r²<{ld_r2}")
        
        # Step 1: LD pruning (optional)
        # Only perform LD pruning if parameters are provided
        if ld_window is not None and ld_step is not None and ld_r2 is not None:
            prune_prefix = f"{output_prefix}_pruned"
            
            # Build LD pruning command
            cmd_prune = [
                'plink2',
                '--bfile', plink_prefix,
            ]
            
            # Add MAF filter if specified
            if maf_threshold is not None:
                cmd_prune.extend(['--maf', str(maf_threshold)])
            
            # Add LD pruning parameters (using variant count, not kb)
            cmd_prune.extend([
                '--indep-pairwise', str(ld_window), str(ld_step), str(ld_r2),
                '--out', prune_prefix
            ])
            
            if verbose:
                logger.info("Step 1: LD pruning...")
            
            result = subprocess.run(cmd_prune, capture_output=True, text=True)
            if result.returncode != 0:
                raise RuntimeError(f"PLINK2 LD pruning failed:\n{result.stderr}")
            
            # Use pruned variants for PCA
            extract_file = f'{prune_prefix}.prune.in'
        else:
            # Skip LD pruning
            if verbose:
                logger.info("Skipping LD pruning (parameters not provided)")
            extract_file = None
        
        # Step 2: Calculate PCA
        cmd_pca = [
            'plink2',
            '--bfile', plink_prefix,
        ]
        
        # Add MAF filter if specified and no LD pruning was done
        if maf_threshold is not None and extract_file is None:
            cmd_pca.extend(['--maf', str(maf_threshold)])
        
        # Add extract file if LD pruning was performed
        if extract_file is not None:
            cmd_pca.extend(['--extract', extract_file])
        
        # Add PCA command
        if approx:
            # Use approximate PCA (PLINK2 will use default approximation)
            cmd_pca.extend([
                '--pca', 'approx', str(n_pcs),
                '--out', output_prefix
            ])
        else:
            # Use exact PCA
            cmd_pca.extend([
                '--pca', str(n_pcs),
                '--out', output_prefix
            ])
        
        if verbose:
            logger.info("Step 2: Calculating PCA...")
        
        result = subprocess.run(cmd_pca, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"PLINK2 PCA calculation failed:\n{result.stderr}")
        
        # Read eigenvec file
        eigenvec_file = f"{output_prefix}.eigenvec"
        pca_df = pd.read_csv(eigenvec_file, sep='\s+', header=None)
        
        # Set column names
        # PLINK2 format: #FID IID PC1 PC2 ... PCn
        pc_cols = [f'PC{i+1}' for i in range(n_pcs)]
        pca_df.columns = ['FID', 'IID'] + pc_cols
        
        # Keep only IID and PC columns
        pca_df = pca_df[['IID'] + pc_cols]
        
        # Convert IID to string
        pca_df['IID'] = pca_df['IID'].astype(str)
        
        # Set IID as index
        pca_df.set_index('IID', inplace=True)
        
        # Read eigenval file for variance explained
        if verbose:
            eigenval_file = f"{output_prefix}.eigenval"
            if os.path.exists(eigenval_file):
                eigenvals = pd.read_csv(eigenval_file, header=None).values.flatten()
                total_var = eigenvals.sum()
                var_explained = eigenvals / total_var
                logger.info(f"Variance explained by first 5 PCs: {var_explained[:5]}")
                logger.info(f"Total variance explained by {n_pcs} PCs: {var_explained.sum():.3f}")
            
            logger.info(f"PCA calculation complete. Found {len(pca_df)} samples.")
        
        return pca_df
        
    finally:
        # Clean up temporary directory
        if cleanup and temp_dir is not None:
            import shutil
            shutil.rmtree(temp_dir, ignore_errors=True)


def calculate_pca_pcair(
    plink_prefix: str,
    n_pcs: int = 10,
    kinship_matrix: Optional[str] = None,
    divergence_matrix: Optional[str] = None,
    output_prefix: Optional[str] = None,
    kin_threshold: float = 0.0884,
    div_threshold: float = -0.0884,
    maf_threshold: float = 0.01,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Calculate PC-AiR (Principal Components - Analysis in Related samples).
    
    PC-AiR is a method for computing principal components that accounts for 
    relatedness and ancestry in genetic data. It identifies an ancestry 
    representative subset of unrelated samples and uses this subset to 
    compute ancestry informative PCs.
    
    Args:
        plink_prefix: Prefix for PLINK binary files (.bed/.bim/.fam)
        n_pcs: Number of principal components to calculate
        kinship_matrix: Path to kinship matrix file (GCTA GRM format prefix)
                       If None, will compute using calculate_grm_gcta()
        divergence_matrix: Path to divergence matrix file (optional)
        output_prefix: Prefix for output files (default: temp directory)
        kin_threshold: Kinship threshold for defining relatedness (default: 0.0884, ~ 2nd degree)
        div_threshold: Divergence threshold (default: -0.0884)
        maf_threshold: MAF threshold for GRM calculation (if kinship_matrix is None)
        verbose: Print progress information
        
    Returns:
        DataFrame with 'IID' column and PC1, PC2, ..., PCn columns
        Index is set to IID
        
    Note:
        Requires R with GENESIS package installed:
        ```R
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
        BiocManager::install("GENESIS")
        BiocManager::install("SNPRelate")
        BiocManager::install("gdsfmt")
        ```
        
    Reference:
        Conomos et al. (2015) Genetic Epidemiology
        https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4608645/
    """
    # Check if R and Rscript are available
    try:
        result = subprocess.run(['Rscript', '--version'], capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError("Rscript not found. Please install R.")
    except FileNotFoundError:
        raise RuntimeError(
            "Rscript not found in PATH. Please install R and ensure it's available.\n"
            "Download R from: https://www.r-project.org/"
        )
    
    # Create temporary directory if no output prefix specified
    if output_prefix is None:
        temp_dir = tempfile.mkdtemp()
        output_prefix = os.path.join(temp_dir, 'pcair')
        cleanup = True
    else:
        temp_dir = None
        cleanup = False
    
    try:
        if verbose:
            logger.info(f"Calculating {n_pcs} PCs using PC-AiR...")
            logger.info(f"Kinship threshold: {kin_threshold}, Divergence threshold: {div_threshold}")
        
        # If kinship matrix not provided, calculate using GCTA
        if kinship_matrix is None:
            if verbose:
                logger.info("Calculating kinship matrix using GCTA...")
            kinship_matrix = calculate_grm_gcta(
                plink_prefix=plink_prefix,
                output_prefix=f"{output_prefix}_grm",
                maf_threshold=maf_threshold,
                verbose=verbose
            )
        
        # Create R script for PC-AiR with proper escaping
        bed_file = shlex.quote(f"{plink_prefix}.bed")
        bim_file = shlex.quote(f"{plink_prefix}.bim")
        fam_file = shlex.quote(f"{plink_prefix}.fam")
        gds_file = shlex.quote(f"{output_prefix}.gds")
        grm_bin_file = shlex.quote(f"{kinship_matrix}.grm.bin")
        grm_n_file = shlex.quote(f"{kinship_matrix}.grm.N.bin")
        grm_id_file = shlex.quote(f"{kinship_matrix}.grm.id")
        output_file = shlex.quote(f"{output_prefix}_pcair.txt")
        variance_file = shlex.quote(f"{output_prefix}_variance.txt")
        unrelated_file = shlex.quote(f"{output_prefix}_unrelated.txt")
        
        r_script = f"""
# Load required libraries
suppressPackageStartupMessages({{
    library(GENESIS)
    library(SNPRelate)
    library(gdsfmt)
}})

# Convert PLINK to GDS format
snpgdsBED2GDS(
    bed.fn = {bed_file},
    bim.fn = {bim_file}, 
    fam.fn = {fam_file},
    out.gdsfn = {gds_file}
)

# Open GDS file
gds <- snpgdsOpen({gds_file})

# Load kinship matrix from GCTA format
# Read GRM binary file
grm_bin <- file({grm_bin_file}, "rb")
grm_n_file_handle <- file({grm_n_file}, "rb")

# Read sample IDs
sample_ids <- read.table({grm_id_file}, header=FALSE, stringsAsFactors=FALSE)
n_samples <- nrow(sample_ids)

# Read number of values (lower triangle including diagonal)
n_values <- n_samples * (n_samples + 1) / 2

# Read GRM values
grm_values <- readBin(grm_bin, what="numeric", n=n_values, size=4)
close(grm_bin)
close(grm_n_file_handle)

# Reconstruct full symmetric matrix
kin_matrix <- matrix(0, nrow=n_samples, ncol=n_samples)
idx <- 1
for(i in 1:n_samples) {{
    for(j in 1:i) {{
        kin_matrix[i,j] <- grm_values[idx]
        kin_matrix[j,i] <- grm_values[idx]
        idx <- idx + 1
    }}
}}

rownames(kin_matrix) <- sample_ids$V2
colnames(kin_matrix) <- sample_ids$V2

# Run PC-AiR
cat("Running PC-AiR analysis...\\n")
pc_air <- pcair(
    gds = gds,
    kinobj = kin_matrix,
    kin.thresh = {kin_threshold},
    divobj = NULL,
    div.thresh = {div_threshold},
    num.cores = 1
)

# Extract PCs
pcs <- pc_air$vectors[, 1:{n_pcs}, drop=FALSE]
colnames(pcs) <- paste0("PC", 1:{n_pcs})

# Save results with IID column
output_df <- data.frame(
    IID = as.character(rownames(pcs)),
    pcs,
    stringsAsFactors = FALSE
)

write.table(
    output_df,
    file = {output_file},
    quote = FALSE,
    row.names = FALSE,
    sep = "\\t"
)

# Save variance explained
write.table(
    data.frame(variance = pc_air$values[1:{n_pcs}]),
    file = {variance_file},
    quote = FALSE,
    row.names = FALSE,
    sep = "\\t"
)

# Save unrelated set
write.table(
    data.frame(IID = as.character(pc_air$unrels)),
    file = {unrelated_file},
    quote = FALSE,
    row.names = FALSE,
    sep = "\\t"
)

cat("PC-AiR complete.\\n")
cat(paste("Unrelated samples:", length(pc_air$unrels), "\\n"))
cat(paste("Related samples:", length(pc_air$rels), "\\n"))

# Close GDS
snpgdsClose(gds)
"""
        
        # Write R script to file
        r_script_file = f"{output_prefix}_pcair.R"
        with open(r_script_file, 'w') as f:
            f.write(r_script)
        
        # Run R script
        if verbose:
            logger.info("Running PC-AiR in R (this may take a while)...")
        
        result = subprocess.run(
            ['Rscript', r_script_file],
            capture_output=True,
            text=True
        )
        
        if result.returncode != 0:
            raise RuntimeError(f"PC-AiR calculation failed:\n{result.stderr}\n{result.stdout}")
        
        # Read results
        pca_file = f"{output_prefix}_pcair.txt"
        if not os.path.exists(pca_file):
            raise RuntimeError(f"PC-AiR output file not found: {pca_file}")
        
        pca_df = pd.read_csv(pca_file, sep='\t')
        
        # Convert IID to string and set as index
        pca_df['IID'] = pca_df['IID'].astype(str)
        pca_df.set_index('IID', inplace=True)
        
        # Read variance explained
        if verbose:
            variance_file = f"{output_prefix}_variance.txt"
            if os.path.exists(variance_file):
                variance = pd.read_csv(variance_file, sep='\t')
                total_var = variance['variance'].sum()
                var_explained = variance['variance'] / total_var
                logger.info(f"Variance explained by first 5 PCs: {var_explained[:5].values}")
                logger.info(f"Total variance explained by {n_pcs} PCs: {var_explained.sum():.3f}")
            
            # Read unrelated samples info
            unrelated_file = f"{output_prefix}_unrelated.txt"
            if os.path.exists(unrelated_file):
                unrelated_df = pd.read_csv(unrelated_file, sep='\t')
                logger.info(f"Number of unrelated samples: {len(unrelated_df)}")
                logger.info(f"Total samples: {len(pca_df)}")
        
        return pca_df
        
    finally:
        # Clean up temporary directory
        if cleanup and temp_dir is not None:
            import shutil
            shutil.rmtree(temp_dir, ignore_errors=True)


def calculate_grm_gcta(
    plink_prefix: str,
    output_prefix: Optional[str] = None,
    maf_threshold: float = 0.01,
    method: str = 'grm',
    max_threads: int = 1,
    verbose: bool = True
) -> str:
    """
    Calculate genetic relationship matrix (GRM) using GCTA.
    
    Args:
        plink_prefix: Prefix for PLINK binary files (.bed/.bim/.fam)
        output_prefix: Prefix for output GRM files (default: temp directory)
        maf_threshold: MAF threshold for variant filtering
        method: GRM calculation method ('grm' for full, 'grm-sparse' for sparse)
        max_threads: Maximum number of threads to use
        verbose: Print progress information
        
    Returns:
        Path to output GRM prefix (files will be prefix.grm.bin, prefix.grm.N.bin, prefix.grm.id)
        
    Note:
        Requires GCTA to be installed and available in PATH.
        Download from: https://yanglab.westlake.edu.cn/software/gcta/
        
    Output files:
        - prefix.grm.bin: GRM values (lower triangle, binary format)
        - prefix.grm.N.bin: Number of SNPs used for each pair
        - prefix.grm.id: Sample IDs (FID and IID)
    """
    # Determine GCTA command (gcta64 or gcta)
    gcta_cmd = 'gcta64'
    test_result = subprocess.run(['which', 'gcta64'], capture_output=True)
    if test_result.returncode != 0:
        gcta_cmd = 'gcta'
        # Test if gcta is available
        test_result = subprocess.run(['which', 'gcta'], capture_output=True)
        if test_result.returncode != 0:
            raise RuntimeError(
                "GCTA not found in PATH. Please install GCTA.\n"
                "Download from: https://yanglab.westlake.edu.cn/software/gcta/"
            )
    
    # Create temporary directory if no output prefix specified
    if output_prefix is None:
        temp_dir = tempfile.mkdtemp()
        output_prefix = os.path.join(temp_dir, 'grm')
    
    if verbose:
        logger.info(f"Calculating GRM using GCTA (method: {method})...")
        logger.info(f"MAF threshold: {maf_threshold}, Threads: {max_threads}")
    
    # Build GCTA command
    cmd = [
        gcta_cmd,
        '--bfile', plink_prefix,
        '--make-' + method,
        '--maf', str(maf_threshold),
        '--thread-num', str(max_threads),
        '--out', output_prefix
    ]
    
    # Run GCTA
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"GCTA GRM calculation failed:\n{result.stderr}")
    
    if verbose:
        logger.info(f"GRM calculation complete. Output: {output_prefix}.grm.*")
    
    return output_prefix


def load_grm_gcta(grm_prefix: str, verbose: bool = True) -> Tuple[np.ndarray, pd.DataFrame]:
    """
    Load GRM calculated by GCTA.
    
    Args:
        grm_prefix: Prefix for GRM files (without .grm.bin extension)
        verbose: Print loading information
        
    Returns:
        Tuple of (grm_matrix, sample_ids_df)
        - grm_matrix: n_samples x n_samples symmetric GRM matrix
        - sample_ids_df: DataFrame with FID and IID columns
        
    Raises:
        FileNotFoundError: If GRM files are not found
    """
    grm_bin_file = f"{grm_prefix}.grm.bin"
    grm_n_file = f"{grm_prefix}.grm.N.bin"
    grm_id_file = f"{grm_prefix}.grm.id"
    
    # Check files exist
    for f in [grm_bin_file, grm_n_file, grm_id_file]:
        if not os.path.exists(f):
            raise FileNotFoundError(f"GRM file not found: {f}")
    
    if verbose:
        logger.info(f"Loading GRM from {grm_prefix}...")
    
    # Load sample IDs
    sample_ids = pd.read_csv(grm_id_file, sep='\t', header=None, names=['FID', 'IID'])
    
    # Convert IID to string for consistency
    sample_ids['IID'] = sample_ids['IID'].astype(str)
    
    n_samples = len(sample_ids)
    
    # Load GRM values (lower triangle)
    grm_values = np.fromfile(grm_bin_file, dtype=np.float32)
    
    # Expected number of values in lower triangle (including diagonal)
    expected_n = n_samples * (n_samples + 1) // 2
    
    if len(grm_values) != expected_n:
        raise ValueError(
            f"GRM file has {len(grm_values)} values, "
            f"expected {expected_n} for {n_samples} samples"
        )
    
    # Reconstruct full symmetric matrix
    grm_matrix = np.zeros((n_samples, n_samples), dtype=np.float32)
    
    idx = 0
    for i in range(n_samples):
        for j in range(i + 1):
            grm_matrix[i, j] = grm_values[idx]
            grm_matrix[j, i] = grm_values[idx]  # Symmetric
            idx += 1
    
    if verbose:
        logger.info(f"Loaded {n_samples} x {n_samples} GRM matrix")
        logger.info(f"Mean diagonal: {np.diag(grm_matrix).mean():.3f}")
        off_diag_mean = (grm_matrix.sum() - np.trace(grm_matrix)) / (n_samples * (n_samples - 1))
        logger.info(f"Mean off-diagonal: {off_diag_mean:.6f}")
    
    return grm_matrix, sample_ids


def attach_pcs_to_phenotype(
    phenotype_df: pd.DataFrame,
    pca_df: pd.DataFrame,
    n_pcs: int = 10,
    pc_prefix: str = 'PC',
    sample_id_col: Optional[str] = None,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Attach principal components to phenotype DataFrame.
    
    Args:
        phenotype_df: Phenotype DataFrame
        pca_df: PCA DataFrame with IID as index and PC columns
        n_pcs: Number of PCs to attach (will use PC1 to PCn)
        pc_prefix: Prefix for PC column names (default: 'PC')
        sample_id_col: Column name in phenotype_df to use for matching with PCA IIDs
                      If None, uses phenotype_df.index
        verbose: Print information about merging
        
    Returns:
        Phenotype DataFrame with PC columns added
        
    Raises:
        ValueError: If requested PCs are not available in pca_df
        
    Example:
        >>> # When phenotype has IID as index
        >>> pheno_df = pd.DataFrame({'age': [25, 30], 'disease': [0, 1]}, 
        ...                         index=['sample1', 'sample2'])
        >>> pca_df = calculate_pca_plink('genotypes')
        >>> pheno_with_pcs = attach_pcs_to_phenotype(pheno_df, pca_df, n_pcs=10)
        
        >>> # When phenotype has IID as a column
        >>> pheno_df = pd.DataFrame({'IID': ['sample1', 'sample2'], 
        ...                          'age': [25, 30], 'disease': [0, 1]})
        >>> pheno_with_pcs = attach_pcs_to_phenotype(pheno_df, pca_df, n_pcs=10, 
        ...                                          sample_id_col='IID')
    """
    # Check if requested PCs are available
    available_pcs = [col for col in pca_df.columns if col.startswith(pc_prefix)]
    max_available = len(available_pcs)
    
    if n_pcs > max_available:
        raise ValueError(
            f"Requested {n_pcs} PCs, but only {max_available} are available in pca_df. "
            f"Available PC columns: {available_pcs}"
        )
    
    # Select requested PCs
    pc_cols = [f'{pc_prefix}{i+1}' for i in range(n_pcs)]
    
    # Check all requested PCs exist
    missing_pcs = [pc for pc in pc_cols if pc not in pca_df.columns]
    if missing_pcs:
        raise ValueError(f"Missing PC columns in pca_df: {missing_pcs}")
    
    pcs_to_add = pca_df[pc_cols].copy()
    
    # Reset PCA index to have IID as a column for merging
    pcs_to_add['IID'] = pcs_to_add.index.astype(str)
    
    # Prepare phenotype dataframe for merging
    result_df = phenotype_df.copy()
    
    if sample_id_col is None:
        # Use index for matching
        result_df['_merge_id'] = result_df.index.astype(str)
        original_index = result_df.index
        use_index = True
    else:
        # Use specified column for matching
        if sample_id_col not in result_df.columns:
            raise ValueError(f"Column '{sample_id_col}' not found in phenotype_df")
        result_df['_merge_id'] = result_df[sample_id_col].astype(str)
        use_index = False
    
    # Merge PCs with phenotype data
    merged_df = result_df.merge(
        pcs_to_add,
        left_on='_merge_id',
        right_on='IID',
        how='left',
        suffixes=('', '_pc')
    )
    
    # Drop merge columns
    merged_df.drop(['_merge_id', 'IID'], axis=1, inplace=True, errors='ignore')
    
    # Restore original index if it was used
    if use_index:
        merged_df.index = original_index
    
    if verbose:
        n_pheno_samples = len(phenotype_df)
        n_pca_samples = len(pca_df)
        n_with_pcs = merged_df[pc_cols[0]].notna().sum()
        
        logger.info(f"Attaching {n_pcs} PCs to phenotype data")
        logger.info(f"Phenotype samples: {n_pheno_samples}")
        logger.info(f"PCA samples: {n_pca_samples}")
        logger.info(f"Samples with PCs after merge: {n_with_pcs}")
        
        if n_with_pcs < n_pheno_samples:
            n_missing = n_pheno_samples - n_with_pcs
            logger.warning(f"{n_missing} samples in phenotype_df have no PCs (will have NA values)")
        
        if n_with_pcs < n_pca_samples:
            n_extra = n_pca_samples - n_with_pcs
            logger.info(f"{n_extra} samples in pca_df are not in phenotype_df")
        
        # Check for missing values
        for pc_col in pc_cols:
            n_missing_pc = merged_df[pc_col].isna().sum()
            if n_missing_pc > 0:
                logger.warning(f"{pc_col}: {n_missing_pc} samples with missing values")
    
    return merged_df


def get_pc_covariate_list(n_pcs: int, pc_prefix: str = 'PC') -> List[str]:
    """
    Generate list of PC covariate names for use in EDGE analysis.
    
    Args:
        n_pcs: Number of PCs
        pc_prefix: Prefix for PC column names (default: 'PC')
        
    Returns:
        List of PC column names ['PC1', 'PC2', ..., 'PCn']
        
    Example:
        >>> pc_list = get_pc_covariate_list(5)
        >>> print(pc_list)
        ['PC1', 'PC2', 'PC3', 'PC4', 'PC5']
        
        >>> # Use in EDGE analysis
        >>> covariates = ['age', 'sex'] + get_pc_covariate_list(10)
        >>> alpha_df, gwas_df = edge.run_full_analysis(
        ...     train_g, train_p, test_g, test_p,
        ...     outcome='disease',
        ...     covariates=covariates
        ... )
    """
    return [f'{pc_prefix}{i+1}' for i in range(n_pcs)]


def identify_related_samples(
    grm_matrix: np.ndarray,
    sample_ids: pd.DataFrame,
    threshold: float = 0.0884,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Identify pairs of related samples based on GRM threshold.
    
    Args:
        grm_matrix: n_samples x n_samples GRM matrix
        sample_ids: DataFrame with sample IDs (from load_grm_gcta)
        threshold: Relatedness threshold (default: 0.0884 ~ 3rd degree relatives)
                  Common thresholds:
                  - 0.354: 1st degree (parent-offspring, full siblings)
                  - 0.177: 2nd degree (half-siblings, grandparent-grandchild)
                  - 0.0884: 3rd degree (first cousins)
        verbose: Print summary statistics
        
    Returns:
        DataFrame with columns: IID1, IID2, kinship
        Sorted by kinship (descending)
        
    Example:
        >>> grm_matrix, sample_ids = load_grm_gcta('output/grm')
        >>> related_pairs = identify_related_samples(grm_matrix, sample_ids, threshold=0.177)
        >>> print(f"Found {len(related_pairs)} related pairs")
    """
    n_samples = len(sample_ids)
    related_pairs = []
    
    # Convert IIDs to strings for consistency
    iids = sample_ids['IID'].astype(str).values
    
    # Find pairs above threshold (exclude diagonal)
    for i in range(n_samples):
        for j in range(i + 1, n_samples):
            kinship = grm_matrix[i, j]
            if kinship >= threshold:
                related_pairs.append({
                    'IID1': iids[i],
                    'IID2': iids[j],
                    'kinship': kinship
                })
    
    # Create DataFrame and sort by kinship
    related_df = pd.DataFrame(related_pairs)
    
    if len(related_df) > 0:
        related_df = related_df.sort_values('kinship', ascending=False)
    
    if verbose:
        logger.info(f"Relatedness threshold: {threshold}")
        logger.info(f"Found {len(related_df)} related pairs out of {n_samples * (n_samples - 1) // 2} possible pairs")
        
        if len(related_df) > 0:
            logger.info(f"Maximum kinship: {related_df['kinship'].max():.4f}")
            logger.info(f"Mean kinship (related pairs): {related_df['kinship'].mean():.4f}")
            
            # Count by degree of relatedness
            first_degree = (related_df['kinship'] >= 0.354).sum()
            second_degree = ((related_df['kinship'] >= 0.177) & (related_df['kinship'] < 0.354)).sum()
            third_degree = ((related_df['kinship'] >= 0.0884) & (related_df['kinship'] < 0.177)).sum()
            
            if first_degree > 0:
                logger.info(f"  1st degree relatives: {first_degree} pairs")
            if second_degree > 0:
                logger.info(f"  2nd degree relatives: {second_degree} pairs")
            if third_degree > 0:
                logger.info(f"  3rd degree relatives: {third_degree} pairs")
    
    return related_df


def filter_related_samples(
    phenotype_df: pd.DataFrame,
    grm_matrix: np.ndarray,
    sample_ids: pd.DataFrame,
    threshold: float = 0.0884,
    method: str = 'greedy',
    sample_id_col: Optional[str] = None,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Filter out related samples to create an unrelated subset.
    
    Args:
        phenotype_df: Phenotype DataFrame
        grm_matrix: n_samples x n_samples GRM matrix
        sample_ids: DataFrame with sample IDs (from load_grm_gcta)
        threshold: Relatedness threshold
        method: Method for selecting unrelated samples
                'greedy': Iteratively remove sample with most relatives
                'random': Randomly remove one from each related pair
        sample_id_col: Column name in phenotype_df for sample IDs
                      If None, uses phenotype_df.index
        verbose: Print filtering information
        
    Returns:
        Filtered phenotype DataFrame with unrelated samples only
        
    Example:
        >>> grm_matrix, sample_ids = load_grm_gcta('output/grm')
        >>> unrelated_pheno = filter_related_samples(
        ...     phenotype_df, grm_matrix, sample_ids, threshold=0.0884
        ... )
    """
    # Get related pairs
    related_df = identify_related_samples(
        grm_matrix, sample_ids, threshold=threshold, verbose=False
    )
    
    n_original = len(phenotype_df)
    
    if len(related_df) == 0:
        if verbose:
            logger.info("No related samples found. Returning original phenotype data.")
        return phenotype_df.copy()
    
    # Get set of all samples involved in relatedness
    all_related = set(related_df['IID1']).union(set(related_df['IID2']))
    
    if method == 'greedy':
        # Greedy algorithm: iteratively remove sample with most relatives
        samples_to_remove = set()
        remaining_pairs = related_df.copy()
        
        while len(remaining_pairs) > 0:
            # Count how many times each sample appears
            sample_counts = pd.concat([
                remaining_pairs['IID1'],
                remaining_pairs['IID2']
            ]).value_counts()
            
            # Remove sample with most relatives
            sample_to_remove = sample_counts.index[0]
            samples_to_remove.add(sample_to_remove)
            
            # Remove all pairs involving this sample
            remaining_pairs = remaining_pairs[
                (remaining_pairs['IID1'] != sample_to_remove) &
                (remaining_pairs['IID2'] != sample_to_remove)
            ]
        
    elif method == 'random':
        # Random method: randomly choose one from each pair
        samples_to_remove = set()
        
        for _, row in related_df.iterrows():
            iid1, iid2 = row['IID1'], row['IID2']
            
            # Skip if already removed
            if iid1 in samples_to_remove or iid2 in samples_to_remove:
                continue
            
            # Randomly choose one to remove
            to_remove = np.random.choice([iid1, iid2])
            samples_to_remove.add(to_remove)
    
    else:
        raise ValueError(f"Unknown method: {method}. Use 'greedy' or 'random'.")
    
    # Filter phenotype data
    if sample_id_col is None:
        # Use index
        pheno_sample_ids = phenotype_df.index.astype(str)
        samples_to_keep = [s for s in pheno_sample_ids if s not in samples_to_remove]
        filtered_df = phenotype_df.loc[samples_to_keep]
    else:
        # Use specified column
        if sample_id_col not in phenotype_df.columns:
            raise ValueError(f"Column '{sample_id_col}' not found in phenotype_df")
        
        pheno_sample_ids = phenotype_df[sample_id_col].astype(str)
        mask = ~pheno_sample_ids.isin(samples_to_remove)
        filtered_df = phenotype_df[mask].copy()
    
    if verbose:
        n_filtered = len(filtered_df)
        n_removed = n_original - n_filtered
        logger.info(f"Filtered related samples using '{method}' method")
        logger.info(f"Original samples: {n_original}")
        logger.info(f"Related samples removed: {n_removed}")
        logger.info(f"Unrelated samples remaining: {n_filtered} ({n_filtered/n_original*100:.1f}%)")
    
    return filtered_df

# Backward compatibility: keep old function name
additive_gwas = standard_gwas






