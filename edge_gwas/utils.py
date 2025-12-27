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


def _ensure_minor_allele_is_alt(
    genotype_df: pd.DataFrame,
    variant_info_df: pd.DataFrame,
    verbose: bool = True
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Ensure minor allele is coded as ALT (2 in 0/1/2 encoding).
    
    For variants where major allele is ALT, flip the encoding:
    - 0 (REF/REF) → 2 (ALT/ALT)
    - 1 (REF/ALT) → 1 (stays same)
    - 2 (ALT/ALT) → 0 (REF/REF)
    
    Also swap REF and ALT allele labels in variant_info.
    
    Args:
        genotype_df: Genotype DataFrame (samples x variants)
        variant_info_df: Variant information DataFrame
        verbose: Print recoding information
        
    Returns:
        Tuple of (recoded_genotype_df, updated_variant_info_df)
    """
    genotype_df = genotype_df.copy()
    variant_info_df = variant_info_df.copy()
    
    flipped_variants = []
    
    for variant_id in genotype_df.columns:
        # Calculate allele frequency for ALT allele
        geno = genotype_df[variant_id]
        valid_geno = geno.dropna()
        
        if len(valid_geno) == 0:
            continue
        
        # ALT allele frequency (currently coded as 2)
        alt_freq = (2 * (valid_geno == 2).sum() + (valid_geno == 1).sum()) / (2 * len(valid_geno))
        
        # If ALT allele frequency > 0.5, it's actually the major allele
        # We need to flip it so minor allele becomes ALT
        if alt_freq > 0.5:
            # Flip genotypes: 0→2, 1→1, 2→0
            genotype_df[variant_id] = geno.map({0: 2, 1: 1, 2: 0})
            
            # Swap REF and ALT in variant info
            if variant_id in variant_info_df.index:
                ref_allele = variant_info_df.loc[variant_id, 'ref_allele']
                alt_allele = variant_info_df.loc[variant_id, 'alt_allele']
                
                variant_info_df.loc[variant_id, 'ref_allele'] = alt_allele
                variant_info_df.loc[variant_id, 'alt_allele'] = ref_allele
            
            flipped_variants.append(variant_id)
    
    if verbose and flipped_variants:
        logger.info(f"Flipped {len(flipped_variants)} variants to ensure minor allele is ALT")
        logger.info(f"Example flipped variants: {flipped_variants[:5]}")
    
    return genotype_df, variant_info_df


def load_plink_data(
    bed_file: str,
    bim_file: str,
    fam_file: str,
    minor_allele_as_alt: bool = True,
    verbose: bool = True
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load PLINK binary format data (.bed/.bim/.fam).
    
    Args:
        bed_file: Path to .bed file
        bim_file: Path to .bim file
        fam_file: Path to .fam file
        minor_allele_as_alt: If True, ensure minor allele is coded as ALT (2)
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
    
    # Create DataFrame: samples as rows, variants as columns
    genotype_df = pd.DataFrame(
        genotypes,  
        index=pd.Index(sample_ids, name='sample_id'),
        columns=pd.Index(variant_ids, name='variant_id')
    )
    
    # Create variant info DataFrame
    # Note: In PLINK BIM format:
    # - a1 is typically the minor/effect allele (should be our ALT)
    # - a0 is typically the major/reference allele (should be our REF)
    variant_info_df = pd.DataFrame({
        'variant_id': variant_ids,
        'chrom': G.coords['chrom'].values,
        'pos': G.coords['pos'].values,
        'ref_allele': G.a0.values,  # Major allele
        'alt_allele': G.a1.values,  # Minor allele
    })
    variant_info_df.set_index('variant_id', inplace=True)
    
    if verbose:
        logger.info(f"Loaded {len(sample_ids)} samples and {len(variant_ids)} variants")
    
    # Ensure minor allele is coded as ALT (2)
    if minor_allele_as_alt:
        if verbose:
            logger.info("Checking and recoding to ensure minor allele is ALT...")
        genotype_df, variant_info_df = _ensure_minor_allele_is_alt(
            genotype_df, variant_info_df, verbose
        )
    
    # Add MAF to variant info
    maf_list = []
    for variant_id in genotype_df.columns:
        geno = genotype_df[variant_id]
        valid_geno = geno.dropna()
        if len(valid_geno) > 0:
            # ALT allele frequency (now guaranteed to be minor)
            alt_freq = (2 * (valid_geno == 2).sum() + (valid_geno == 1).sum()) / (2 * len(valid_geno))
            maf = min(alt_freq, 1 - alt_freq)  # Should equal alt_freq after recoding
        else:
            maf = np.nan
        maf_list.append(maf)
    
    variant_info_df['MAF'] = maf_list
    
    if verbose:
        logger.info(f"MAF range: {variant_info_df['MAF'].min():.4f} - {variant_info_df['MAF'].max():.4f}")
    
    return genotype_df, variant_info_df


def load_pgen_data(
    pgen_file: str,
    pvar_file: str,
    psam_file: str,
    minor_allele_as_alt: bool = True,
    verbose: bool = True
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load PLINK 2 binary format data (.pgen/.pvar/.psam).
    
    Args:
        pgen_file: Path to .pgen file
        pvar_file: Path to .pvar file
        psam_file: Path to .psam file
        minor_allele_as_alt: If True, ensure minor allele is coded as ALT (2)
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
    
    # Ensure minor allele is coded as ALT (2)
    if minor_allele_as_alt:
        if verbose:
            logger.info("Checking and recoding to ensure minor allele is ALT...")
        genotype_df, variant_info_df = _ensure_minor_allele_is_alt(
            genotype_df, variant_info_df, verbose
        )
    
    # Add MAF to variant info
    maf_list = []
    for variant_id in genotype_df.columns:
        geno = genotype_df[variant_id]
        valid_geno = geno.dropna()
        if len(valid_geno) > 0:
            alt_freq = (2 * (valid_geno == 2).sum() + (valid_geno == 1).sum()) / (2 * len(valid_geno))
            maf = min(alt_freq, 1 - alt_freq)
        else:
            maf = np.nan
        maf_list.append(maf)
    
    variant_info_df['MAF'] = maf_list
    
    if verbose:
        logger.info(f"MAF range: {variant_info_df['MAF'].min():.4f} - {variant_info_df['MAF'].max():.4f}")
    
    return genotype_df, variant_info_df


def load_bgen_data(
    bgen_file: str,
    sample_file: Optional[str] = None,
    minor_allele_as_alt: bool = True,
    verbose: bool = True
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load BGEN format data.
    
    Args:
        bgen_file: Path to .bgen file
        sample_file: Path to .sample file (optional, can be embedded in BGEN)
        minor_allele_as_alt: If True, ensure minor allele is coded as ALT (2)
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
    probs = bgen.read()
    
    # Convert probabilities to dosages (expected allele count)
    # Dosage = P(0/1) * 1 + P(1/1) * 2
    dosages = probs[:, :, 1] + 2 * probs[:, :, 2]
    
    # Get sample IDs
    if sample_file is not None:
        sample_df = pd.read_csv(sample_file, sep=' ', skiprows=2, header=None)
        sample_ids = sample_df[0].values
    else:
        sample_ids = bgen.samples
    
    # Get variant information
    variant_ids = bgen.ids
    chromosomes = bgen.chromosomes
    positions = bgen.positions
    ref_alleles = bgen.allele_ids[:, 0]
    alt_alleles = bgen.allele_ids[:, 1]
    
    # Create genotype DataFrame (samples x variants)
    genotype_df = pd.DataFrame(
        dosages.T,
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
    
    # Ensure minor allele is coded as ALT (2)
    if minor_allele_as_alt:
        if verbose:
            logger.info("Checking and recoding to ensure minor allele is ALT...")
        
        # For dosage data, we need to flip dosages: dosage_new = 2 - dosage_old
        flipped_variants = []
        
        for variant_id in genotype_df.columns:
            geno = genotype_df[variant_id]
            valid_geno = geno.dropna()
            
            if len(valid_geno) == 0:
                continue
            
            # Calculate mean dosage (= 2 * ALT allele frequency)
            mean_dosage = valid_geno.mean()
            alt_freq = mean_dosage / 2
            
            # If ALT allele frequency > 0.5, flip
            if alt_freq > 0.5:
                # Flip dosages: new_dosage = 2 - old_dosage
                genotype_df[variant_id] = 2 - geno
                
                # Swap REF and ALT in variant info
                if variant_id in variant_info_df.index:
                    ref_allele = variant_info_df.loc[variant_id, 'ref_allele']
                    alt_allele = variant_info_df.loc[variant_id, 'alt_allele']
                    
                    variant_info_df.loc[variant_id, 'ref_allele'] = alt_allele
                    variant_info_df.loc[variant_id, 'alt_allele'] = ref_allele
                
                flipped_variants.append(variant_id)
        
        if verbose and flipped_variants:
            logger.info(f"Flipped {len(flipped_variants)} variants to ensure minor allele is ALT")
            logger.info(f"Example flipped variants: {flipped_variants[:5]}")
    
    # Add MAF to variant info
    maf_list = []
    for variant_id in genotype_df.columns:
        geno = genotype_df[variant_id]
        valid_geno = geno.dropna()
        if len(valid_geno) > 0:
            # For dosage data, ALT freq = mean(dosage) / 2
            alt_freq = valid_geno.mean() / 2
            maf = min(alt_freq, 1 - alt_freq)
        else:
            maf = np.nan
        maf_list.append(maf)
    
    variant_info_df['MAF'] = maf_list
    
    if verbose:
        logger.info(f"MAF range: {variant_info_df['MAF'].min():.4f} - {variant_info_df['MAF'].max():.4f}")
    
    return genotype_df, variant_info_df


def load_vcf_data(
    vcf_file: str,
    dosage: bool = True,
    minor_allele_as_alt: bool = True,
    verbose: bool = True
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load VCF format data.
    
    Args:
        vcf_file: Path to .vcf or .vcf.gz file
        dosage: If True, use dosages (DS field); if False, use hard calls (GT field)
        minor_allele_as_alt: If True, ensure minor allele is coded as ALT (2)
        verbose: Print loading information
        
    Returns:
        Tuple of (genotype_df, variant_info_df)
        
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
            ds_values = variant.format('DS')
            # Handle both 1D and 2D arrays
            if ds_values.ndim == 2:
                geno = ds_values[:, 0]  # Get dosage for first alt allele
            else:
                geno = ds_values  # Already 1D
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
    
    # Ensure minor allele is coded as ALT (2)
    if minor_allele_as_alt:
        if verbose:
            logger.info("Checking and recoding to ensure minor allele is ALT...")
        
        if dosage:
            # For dosage data, flip: dosage_new = 2 - dosage_old
            flipped_variants = []
            
            for variant_id in genotype_df.columns:
                geno = genotype_df[variant_id]
                valid_geno = geno.dropna()
                
                if len(valid_geno) == 0:
                    continue
                
                # Calculate mean dosage - ensure it's a scalar
                mean_dosage = np.mean(valid_geno.values)  # FIXED: Use np.mean on .values
                alt_freq = mean_dosage / 2
                
                # If ALT allele frequency > 0.5, flip
                if alt_freq > 0.5:
                    genotype_df[variant_id] = 2 - geno
                    
                    # Swap REF and ALT in variant info
                    if variant_id in variant_info_df.index:
                        ref_allele = variant_info_df.loc[variant_id, 'ref_allele']
                        alt_allele = variant_info_df.loc[variant_id, 'alt_allele']
                        
                        variant_info_df.loc[variant_id, 'ref_allele'] = alt_allele
                        variant_info_df.loc[variant_id, 'alt_allele'] = ref_allele
                    
                    flipped_variants.append(variant_id)
            
            if verbose and flipped_variants:
                logger.info(f"Flipped {len(flipped_variants)} variants (dosage)")
                logger.info(f"Example flipped variants: {flipped_variants[:5]}")
        else:
            # For hard calls, use the standard function
            genotype_df, variant_info_df = _ensure_minor_allele_is_alt(
                genotype_df, variant_info_df, verbose
            )
    
    # Add MAF to variant info
    maf_list = []
    for variant_id in genotype_df.columns:
        geno = genotype_df[variant_id]
        valid_geno = geno.dropna()
        if len(valid_geno) > 0:
            if dosage:
                # For dosage data - ensure scalar
                alt_freq = np.mean(valid_geno.values) / 2  # FIXED: Use np.mean on .values
            else:
                # For hard calls
                alt_freq = (2 * (valid_geno == 2).sum() + (valid_geno == 1).sum()) / (2 * len(valid_geno))
            maf = min(alt_freq, 1 - alt_freq)
        else:
            maf = np.nan
        maf_list.append(maf)
    
    variant_info_df['MAF'] = maf_list
    
    if verbose:
        logger.info(f"MAF range: {variant_info_df['MAF'].min():.4f} - {variant_info_df['MAF'].max():.4f}")
    
    return genotype_df, variant_info_df


def validate_genotype_df(
    genotype_df: pd.DataFrame,
    variant_info_df: Optional[pd.DataFrame] = None,
    name: str = "genotype_df",
    check_encoding: bool = True,
    verbose: bool = True,
    return_details: bool = False
) -> Union[None, bool, Tuple[bool, pd.DataFrame]]:
    """
    Validate genotype DataFrame format and encoding.
    
    Performs:
    1. Basic format validation (type, not empty, no duplicates)
    2. Encoding validation (minor allele as ALT, valid values) - if variant_info provided
    
    Args:
        genotype_df: Genotype DataFrame (samples x variants)
        variant_info_df: Optional variant information DataFrame
        name: Name for error messages
        check_encoding: If True, validate encoding (requires variant_info_df)
        verbose: Print validation results
        return_details: If True, return (passed, report_df)
        
    Returns:
        None (raises errors if invalid) OR
        bool (validation passed) OR  
        Tuple[bool, pd.DataFrame] (if return_details=True)
        
    Raises:
        TypeError: If not a pandas DataFrame
        ValueError: If empty, duplicates, or invalid format
    """
    # ========================================================================
    # PART 1: Basic Format Validation
    # ========================================================================
    
    if not isinstance(genotype_df, pd.DataFrame):
        raise TypeError(f"{name} must be a pandas DataFrame, got {type(genotype_df)}")
    
    if genotype_df.empty:
        raise ValueError(f"{name} is empty")
    
    # Check duplicate sample IDs
    if genotype_df.index.duplicated().any():
        n_dup = genotype_df.index.duplicated().sum()
        examples = genotype_df.index[genotype_df.index.duplicated()].unique()[:5].tolist()
        raise ValueError(f"{name} has {n_dup} duplicate sample IDs. Examples: {examples}")
    
    # Check duplicate variant IDs
    if genotype_df.columns.duplicated().any():
        n_dup = genotype_df.columns.duplicated().sum()
        examples = genotype_df.columns[genotype_df.columns.duplicated()].unique()[:5].tolist()
        raise ValueError(f"{name} has {n_dup} duplicate variant IDs. Examples: {examples}")
    
    n_samples = len(genotype_df)
    n_variants = len(genotype_df.columns)
    
    if verbose:
        logger.info(f"✓ {name} basic validation passed: {n_samples} samples, {n_variants} variants")
    
    # If no encoding check needed, return
    if not check_encoding or variant_info_df is None:
        return None
    
    # ========================================================================
    # PART 2: Encoding Validation
    # ========================================================================
    
    if verbose:
        logger.info("Checking encoding (minor allele as ALT)...")
    
    # Detect data type
    all_vals = genotype_df.values.flatten()
    all_vals = all_vals[~np.isnan(all_vals)]
    is_dosage = not np.all(np.isin(all_vals, [0, 1, 2]))
    
    validation_results = []
    n_fail = 0
    n_warning = 0
    
    for variant_id in genotype_df.columns:
        geno = genotype_df[variant_id].dropna()
        
        result = {'variant_id': variant_id, 'n_valid': len(geno)}
        
        if len(geno) == 0:
            result['status'] = 'WARNING'
            result['issue'] = 'No valid genotypes'
            n_warning += 1
            validation_results.append(result)
            continue
        
        # Calculate frequencies
        if is_dosage:
            alt_freq = geno.mean() / 2
        else:
            alt_freq = (2 * (geno == 2).sum() + (geno == 1).sum()) / (2 * len(geno))
        
        maf = min(alt_freq, 1 - alt_freq)
        result['alt_freq'] = alt_freq
        result['maf'] = maf
        result['minor_is_alt'] = alt_freq <= 0.5
        
        # Validation
        if alt_freq > 0.5:
            result['status'] = 'FAIL'
            result['issue'] = f'Minor allele is REF (ALT freq={alt_freq:.3f})'
            n_fail += 1
        else:
            result['status'] = 'PASS'
            result['issue'] = None
        
        validation_results.append(result)
    
    results_df = pd.DataFrame(validation_results)
    validation_passed = n_fail == 0
    
    if verbose:
        n_pass = (results_df['status'] == 'PASS').sum()
        logger.info(f"Encoding validation: PASS={n_pass}, FAIL={n_fail}, WARNING={n_warning}")
        logger.info(f"Overall: {'✓ PASSED' if validation_passed else '✗ FAILED'}")
        
        if not results_df.empty and 'maf' in results_df.columns:
            valid_maf = results_df['maf'].dropna()
            if len(valid_maf) > 0:
                logger.info(f"MAF range: {valid_maf.min():.4f} - {valid_maf.max():.4f}")
    
    if return_details:
        return validation_passed, results_df
    else:
        return validation_passed


def validate_and_fix_encoding(
    genotype_df: pd.DataFrame,
    variant_info_df: pd.DataFrame,
    verbose: bool = True
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Validate and automatically fix genotype encoding.
    
    Args:
        genotype_df: Genotype DataFrame
        variant_info_df: Variant info DataFrame
        verbose: Print progress
        
    Returns:
        Tuple of (fixed_genotype_df, fixed_variant_info_df, report_df)
    """
    # Validate current encoding
    is_valid, initial_report = validate_genotype_df(
        genotype_df, variant_info_df, 
        check_encoding=True, verbose=verbose, return_details=True
    )
    
    if is_valid:
        if verbose:
            logger.info("✓ No fixes needed")
        initial_report['was_fixed'] = False
        return genotype_df.copy(), variant_info_df.copy(), initial_report
    
    # Fix encoding
    if verbose:
        logger.info("Applying fixes...")
    
    geno_fixed, info_fixed = _ensure_minor_allele_is_alt(
        genotype_df, variant_info_df, verbose=verbose
    )
    
    # Re-validate
    is_valid_after, final_report = validate_genotype_df(
        geno_fixed, info_fixed,
        check_encoding=True, verbose=verbose, return_details=True
    )
    
    # Mark fixed variants
    failed_vars = set(initial_report[initial_report['status'] == 'FAIL']['variant_id'])
    final_report['was_fixed'] = final_report['variant_id'].isin(failed_vars)
    
    if verbose:
        n_fixed = final_report['was_fixed'].sum()
        logger.info(f"✓ Fixed {n_fixed} variants")
        if is_valid_after:
            logger.info("✓ All checks now pass!")
    
    return geno_fixed, info_fixed, final_report


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
        sample_id_col: Name of sample ID column (will become index)
        sep: File separator (default: tab)
        log_transform_outcome: Apply log10(x+1) transformation to outcome
        
    Returns:
        DataFrame with sample IDs as index, outcome and covariates as columns
        
    Example:
        >>> pheno = prepare_phenotype_data(
        ...     'pheno.txt',
        ...     outcome_col='disease',
        ...     covariate_cols=['age', 'sex'],
        ...     sep=' '
        ... )
    """
    logger.info(f"Loading phenotype data from {phenotype_file}")
    
    # Read phenotype file
    pheno_df = pd.read_csv(phenotype_file, sep=sep)
    
    logger.info(f"Loaded {len(pheno_df)} samples with columns: {list(pheno_df.columns)}")
    
    # Check sample ID column exists
    if sample_id_col not in pheno_df.columns:
        raise ValueError(
            f"Sample ID column '{sample_id_col}' not found. "
            f"Available columns: {list(pheno_df.columns)}"
        )
    
    # Validate required columns BEFORE setting index
    validate_phenotype_df(pheno_df, outcome_col, covariate_cols)
    
    # Set sample ID as index
    pheno_df.set_index(sample_id_col, inplace=True)
    
    # Check for duplicate sample IDs
    if pheno_df.index.duplicated().any():
        n_dup = pheno_df.index.duplicated().sum()
        examples = pheno_df.index[pheno_df.index.duplicated()].unique()[:5].tolist()
        raise ValueError(
            f"Phenotype data has {n_dup} duplicate sample IDs. "
            f"Examples: {examples}"
        )
    
    # Select outcome and covariates
    selected_cols = [outcome_col] + covariate_cols
    pheno_df = pheno_df[selected_cols]
    
    # Log transform outcome if requested
    if log_transform_outcome:
        logger.info(f"Applying log10(x+1) transformation to {outcome_col}")
        pheno_df[outcome_col] = np.log10(pheno_df[outcome_col] + 1)
    
    # Remove missing values
    n_before = len(pheno_df)
    pheno_df = pheno_df.dropna()
    n_after = len(pheno_df)
    
    if n_before > n_after:
        logger.info(f"Removed {n_before - n_after} samples with missing data")
    
    logger.info(f"✓ Prepared phenotype data for {n_after} samples")
    
    return pheno_df


def validate_phenotype_df(
    phenotype_df: pd.DataFrame,
    outcome_col: str,
    covariate_cols: List[str],
    name: str = "phenotype_df"
) -> None:
    """
    Validate phenotype DataFrame format.
    
    Args:
        phenotype_df: Phenotype DataFrame to validate
        outcome_col: Name of outcome column
        covariate_cols: List of covariate column names
        name: Name of the DataFrame for error messages
        
    Raises:
        TypeError: If not a pandas DataFrame
        ValueError: If required columns are missing or DataFrame is invalid
    """
    if not isinstance(phenotype_df, pd.DataFrame):
        raise TypeError(f"{name} must be a pandas DataFrame, got {type(phenotype_df)}")
    
    if phenotype_df.empty:
        raise ValueError(f"{name} is empty")
    
    # Check outcome column exists
    if outcome_col not in phenotype_df.columns:
        raise ValueError(
            f"Outcome column '{outcome_col}' not found in {name}. "
            f"Available columns: {list(phenotype_df.columns)}"
        )
    
    # Check covariate columns exist
    missing_covariates = [col for col in covariate_cols if col not in phenotype_df.columns]
    if missing_covariates:
        raise ValueError(
            f"Covariate columns {missing_covariates} not found in {name}. "
            f"Available columns: {list(phenotype_df.columns)}"
        )
    
    logger.info(f"✓ {name} validation passed")


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
    # For hard calls: mean/2 = allele frequency
    # For dosages: mean/2 = allele frequency
    af = genotype_df.mean(axis=0, skipna=True) / 2  # Allele frequency
    
    # Handle NaN values (variants with all missing data)
    af = af.fillna(0)
    
    # MAF is min(af, 1-af)
    maf = np.minimum(af, 1 - af)
    
    # Filter variants
    variants_to_keep = maf >= min_maf
    genotype_df_filtered = genotype_df.loc[:, variants_to_keep]
    
    n_after = genotype_df_filtered.shape[1]
    
    if verbose:
        logger.info(f"Filtered variants by MAF >= {min_maf}")
        logger.info(f"Kept {n_after}/{n_before} variants ({n_after/n_before*100:.1f}%)")
        if n_after > 0:
            logger.info(f"MAF range in filtered data: {maf[variants_to_keep].min():.4f} - {maf[variants_to_keep].max():.4f}")
    
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
        max_missing: Maximum proportion of missing genotypes allowed (0-1)
        verbose: Print filtering information
        
    Returns:
        Filtered genotype DataFrame
    """
    n_before = genotype_df.shape[1]
    
    # Calculate missing rate for each variant
    missing_rates = genotype_df.isna().mean(axis=0)
    
    # Filter variants
    variants_to_keep = missing_rates <= max_missing
    genotype_df_filtered = genotype_df.loc[:, variants_to_keep]
    
    n_after = genotype_df_filtered.shape[1]
    
    if verbose:
        logger.info(f"Filtered variants by missing rate <= {max_missing}")
        logger.info(f"Kept {n_after}/{n_before} variants ({n_after/n_before*100:.1f}%)")
        if n_after > 0:
            remaining_missing = missing_rates[variants_to_keep]
            logger.info(f"Missing rate range in filtered data: {remaining_missing.min():.4f} - {remaining_missing.max():.4f}")
    
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
        genotype_df: Genotype DataFrame (samples as index)
        phenotype_df: Phenotype DataFrame (sample IDs as index)
        min_call_rate: Minimum call rate (proportion of non-missing genotypes, 0-1)
        verbose: Print filtering information
        
    Returns:
        Tuple of (filtered_genotype_df, filtered_phenotype_df)
    """
    n_samples_before = genotype_df.shape[0]
    n_pheno_before = phenotype_df.shape[0]
    
    # Calculate call rate for each sample (proportion of non-missing)
    sample_call_rate = genotype_df.notna().mean(axis=1)
    
    # Filter samples
    samples_to_keep = sample_call_rate >= min_call_rate
    good_samples = genotype_df.index[samples_to_keep]
    
    # Filter genotype data
    genotype_df_filtered = genotype_df.loc[samples_to_keep, :]
    
    # Filter phenotype data - keep samples that are in good_samples
    # Convert to strings for matching
    good_samples_str = set(good_samples.astype(str))
    pheno_index_str = phenotype_df.index.astype(str)
    pheno_mask = pheno_index_str.isin(good_samples_str)
    phenotype_df_filtered = phenotype_df[pheno_mask].copy()
    
    n_samples_after = genotype_df_filtered.shape[0]
    n_pheno_after = phenotype_df_filtered.shape[0]
    
    if verbose:
        logger.info(f"Filtered samples by call rate >= {min_call_rate}")
        logger.info(f"Genotype samples: kept {n_samples_after}/{n_samples_before} ({n_samples_after/n_samples_before*100:.1f}%)")
        logger.info(f"Phenotype samples: kept {n_pheno_after}/{n_pheno_before} ({n_pheno_after/n_pheno_before*100:.1f}%)")
        if n_samples_after > 0:
            remaining_call_rates = sample_call_rate[samples_to_keep]
            logger.info(f"Call rate range: {remaining_call_rates.min():.4f} - {remaining_call_rates.max():.4f}")
    
    return genotype_df_filtered, phenotype_df_filtered


def filter_genotype_data(
    genotype_df: pd.DataFrame,
    phenotype_df: Optional[pd.DataFrame] = None,
    min_maf: Optional[float] = None,
    max_missing_per_variant: Optional[float] = None,
    min_call_rate_per_sample: Optional[float] = None,
    verbose: bool = True
) -> Union[pd.DataFrame, Tuple[pd.DataFrame, pd.DataFrame]]:
    """
    Comprehensive genotype data filtering with multiple QC criteria.
    
    Applies filters in this order:
    1. Variant MAF filter (if min_maf specified)
    2. Variant missing rate filter (if max_missing_per_variant specified)
    3. Sample call rate filter (if min_call_rate_per_sample specified)
    
    Args:
        genotype_df: Genotype DataFrame (samples x variants)
        phenotype_df: Optional phenotype DataFrame (required if filtering samples)
        min_maf: Minimum minor allele frequency (e.g., 0.01 for 1%)
                If None, no MAF filtering
        max_missing_per_variant: Maximum missing rate per variant (e.g., 0.1 for 10%)
                                If None, no variant missing rate filtering
        min_call_rate_per_sample: Minimum call rate per sample (e.g., 0.95 for 95%)
                                 If None, no sample filtering
        verbose: Print filtering information
        
    Returns:
        If min_call_rate_per_sample is None: filtered_genotype_df
        If min_call_rate_per_sample is specified: (filtered_genotype_df, filtered_phenotype_df)
        
    Examples:
        >>> # Filter by MAF only
        >>> geno_filtered = filter_genotype_data(geno, min_maf=0.01)
        
        >>> # Filter by MAF and missing rate
        >>> geno_filtered = filter_genotype_data(
        ...     geno, min_maf=0.01, max_missing_per_variant=0.1
        ... )
        
        >>> # Filter variants and samples
        >>> geno_filtered, pheno_filtered = filter_genotype_data(
        ...     geno, pheno, 
        ...     min_maf=0.01, 
        ...     max_missing_per_variant=0.1,
        ...     min_call_rate_per_sample=0.95
        ... )
        
        >>> # No filtering (just returns copy)
        >>> geno_copy = filter_genotype_data(geno)
    """
    if verbose:
        logger.info("="*70)
        logger.info("Starting genotype QC filtering")
        logger.info("="*70)
        logger.info(f"Input: {genotype_df.shape[0]} samples x {genotype_df.shape[1]} variants")
    
    # Start with copy of input data
    geno_filtered = genotype_df.copy()
    pheno_filtered = phenotype_df.copy() if phenotype_df is not None else None
    
    n_variants_initial = geno_filtered.shape[1]
    n_samples_initial = geno_filtered.shape[0]
    
    # Track filtering steps
    filter_steps = []
    
    # Step 1: MAF filter
    if min_maf is not None:
        if verbose:
            logger.info(f"\n[1/3] Filtering variants by MAF >= {min_maf}")
        
        n_before = geno_filtered.shape[1]
        geno_filtered = filter_variants_by_maf(geno_filtered, min_maf=min_maf, verbose=verbose)
        n_after = geno_filtered.shape[1]
        
        filter_steps.append({
            'step': 'MAF filter',
            'criterion': f'>= {min_maf}',
            'before': n_before,
            'after': n_after,
            'removed': n_before - n_after
        })
    else:
        if verbose:
            logger.info("\n[1/3] Skipping MAF filter (not specified)")
    
    # Step 2: Variant missing rate filter
    if max_missing_per_variant is not None:
        if verbose:
            logger.info(f"\n[2/3] Filtering variants by missing rate <= {max_missing_per_variant}")
        
        n_before = geno_filtered.shape[1]
        geno_filtered = filter_variants_by_missing(
            geno_filtered, 
            max_missing=max_missing_per_variant, 
            verbose=verbose
        )
        n_after = geno_filtered.shape[1]
        
        filter_steps.append({
            'step': 'Variant missing',
            'criterion': f'<= {max_missing_per_variant}',
            'before': n_before,
            'after': n_after,
            'removed': n_before - n_after
        })
    else:
        if verbose:
            logger.info("\n[2/3] Skipping variant missing rate filter (not specified)")
    
    # Step 3: Sample call rate filter
    if min_call_rate_per_sample is not None:
        if phenotype_df is None:
            raise ValueError(
                "phenotype_df must be provided when filtering by sample call rate"
            )
        
        if verbose:
            logger.info(f"\n[3/3] Filtering samples by call rate >= {min_call_rate_per_sample}")
        
        n_samples_before = geno_filtered.shape[0]
        n_pheno_before = pheno_filtered.shape[0]
        
        geno_filtered, pheno_filtered = filter_samples_by_call_rate(
            geno_filtered,
            pheno_filtered,
            min_call_rate=min_call_rate_per_sample,
            verbose=verbose
        )
        
        n_samples_after = geno_filtered.shape[0]
        n_pheno_after = pheno_filtered.shape[0]
        
        filter_steps.append({
            'step': 'Sample call rate',
            'criterion': f'>= {min_call_rate_per_sample}',
            'before': n_samples_before,
            'after': n_samples_after,
            'removed': n_samples_before - n_samples_after
        })
    else:
        if verbose:
            logger.info("\n[3/3] Skipping sample call rate filter (not specified)")
    
    # Summary
    if verbose:
        logger.info("\n" + "="*70)
        logger.info("FILTERING SUMMARY")
        logger.info("="*70)
        
        if filter_steps:
            summary_df = pd.DataFrame(filter_steps)
            print(summary_df.to_string(index=False))
            print()
        
        n_variants_final = geno_filtered.shape[1]
        n_samples_final = geno_filtered.shape[0]
        
        logger.info(f"Variants: {n_variants_initial} → {n_variants_final} "
                   f"({n_variants_final/n_variants_initial*100:.1f}% retained)")
        logger.info(f"Samples:  {n_samples_initial} → {n_samples_final} "
                   f"({n_samples_final/n_samples_initial*100:.1f}% retained)")
        
        if pheno_filtered is not None:
            logger.info(f"Phenotype samples: {pheno_filtered.shape[0]}")
        
        logger.info("="*70)
    
    # Return based on whether sample filtering was done
    if min_call_rate_per_sample is not None:
        return geno_filtered, pheno_filtered
    else:
        return geno_filtered
        

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
    file_prefix: str,
    n_pcs: int = 10,
    file_format: str = 'bfile',
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
        file_prefix: Prefix for input files
        n_pcs: Number of principal components to calculate
        file_format: Input file format ('bfile', 'pfile', 'vcf', 'bgen'). Default: 'bfile'
        output_prefix: Prefix for output files (default: temp directory)
        maf_threshold: MAF threshold for variant filtering (default: 0.01, None to skip)
        ld_window: Window size for LD pruning in variant count (default: 50, None to skip)
        ld_step: Step size for LD pruning in variant count (default: 5, None to skip)
        ld_r2: R² threshold for LD pruning (default: 0.2, None to skip)
        approx: Use approximate PCA for large cohorts
        verbose: Print progress information
        
    Returns:
        DataFrame with IID as index and PC1, PC2, ..., PCn columns
    """
    # Validate file format
    valid_formats = ['bfile', 'pfile', 'vcf', 'bgen']
    if file_format not in valid_formats:
        raise ValueError(f"file_format must be one of {valid_formats}, got '{file_format}'")
    
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
            logger.info(f"Input format: {file_format}")
            if maf_threshold is not None:
                logger.info(f"MAF threshold: {maf_threshold}")
            if ld_window is not None and ld_r2 is not None:
                logger.info(f"LD pruning: window={ld_window}, step={ld_step}, r²<{ld_r2}")
        
        # Build file input flag based on format
        if file_format == 'bfile':
            file_flag = ['--bfile', file_prefix]
        elif file_format == 'pfile':
            file_flag = ['--pfile', file_prefix]
        elif file_format == 'vcf':
            file_flag = ['--vcf', file_prefix]
        elif file_format == 'bgen':
            file_flag = ['--bgen', file_prefix]
        
        # Step 1: LD pruning (optional)
        if ld_window is not None and ld_step is not None and ld_r2 is not None:
            prune_prefix = f"{output_prefix}_pruned"
            cmd_prune = ['plink2'] + file_flag
            
            if maf_threshold is not None:
                cmd_prune.extend(['--maf', str(maf_threshold)])
            
            cmd_prune.extend([
                '--indep-pairwise', str(ld_window), str(ld_step), str(ld_r2),
                '--out', prune_prefix
            ])
            
            if verbose:
                logger.info("Step 1: LD pruning...")
            
            result = subprocess.run(cmd_prune, capture_output=True, text=True)
            if result.returncode != 0:
                raise RuntimeError(f"PLINK2 LD pruning failed:\n{result.stderr}")
            
            extract_file = f'{prune_prefix}.prune.in'
        else:
            if verbose:
                logger.info("Skipping LD pruning")
            extract_file = None
        
        # Step 2: Calculate PCA
        cmd_pca = ['plink2'] + file_flag
        
        if maf_threshold is not None and extract_file is None:
            cmd_pca.extend(['--maf', str(maf_threshold)])
        
        if extract_file is not None:
            cmd_pca.extend(['--extract', extract_file])
        
        if approx:
            cmd_pca.extend(['--pca', 'approx', str(n_pcs), '--out', output_prefix])
        else:
            cmd_pca.extend(['--pca', str(n_pcs), '--out', output_prefix])
        
        if verbose:
            logger.info("Step 2: Calculating PCA...")
        
        result = subprocess.run(cmd_pca, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"PLINK2 PCA calculation failed:\n{result.stderr}")
        
        # Read results - FIXED: skip comment lines starting with #
        eigenvec_file = f"{output_prefix}.eigenvec"
        pca_df = pd.read_csv(eigenvec_file, sep='\s+', comment='#', header=None, engine='python')
        
        # Set column names
        pc_cols = [f'PC{i+1}' for i in range(n_pcs)]
        pca_df.columns = ['FID', 'IID'] + pc_cols
        pca_df = pca_df[['IID'] + pc_cols]
        pca_df['IID'] = pca_df['IID'].astype(str)
        pca_df.set_index('IID', inplace=True)
        
        if verbose:
            eigenval_file = f"{output_prefix}.eigenval"
            if os.path.exists(eigenval_file):
                eigenvals = pd.read_csv(eigenval_file, header=None).values.flatten()
                var_explained = eigenvals / eigenvals.sum()
                logger.info(f"Variance explained by first 5 PCs: {var_explained[:5]}")
                logger.info(f"Total variance explained: {var_explained.sum():.3f}")
            logger.info(f"PCA complete. Found {len(pca_df)} samples.")
        
        return pca_df
        
    finally:
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
    """Calculate PC-AiR (Principal Components - Analysis in Related samples)."""
    
    try:
        result = subprocess.run(['Rscript', '--version'], capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError("Rscript not found. Please install R.")
    except FileNotFoundError:
        raise RuntimeError("Rscript not found in PATH. Please install R.")
    
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
        
        if kinship_matrix is None:
            if verbose:
                logger.info("Calculating kinship matrix using GCTA...")
            kinship_matrix = calculate_grm_gcta(
                plink_prefix=plink_prefix,
                output_prefix=f"{output_prefix}_grm",
                maf_threshold=maf_threshold,
                verbose=verbose
            )
        
        bed_file = os.path.abspath(f"{plink_prefix}.bed")
        bim_file = os.path.abspath(f"{plink_prefix}.bim")
        fam_file = os.path.abspath(f"{plink_prefix}.fam")
        gds_file = os.path.abspath(f"{output_prefix}.gds")
        grm_bin_file = os.path.abspath(f"{kinship_matrix}.grm.bin")
        grm_n_file = os.path.abspath(f"{kinship_matrix}.grm.N.bin")
        grm_id_file = os.path.abspath(f"{kinship_matrix}.grm.id")
        output_file = os.path.abspath(f"{output_prefix}_pcair.txt")
        variance_file = os.path.abspath(f"{output_prefix}_variance.txt")
        unrelated_file = os.path.abspath(f"{output_prefix}_unrelated.txt")
        
        r_script = f"""
suppressPackageStartupMessages({{
    library(GENESIS)
    library(SNPRelate)
    library(gdsfmt)
}})

snpgdsBED2GDS(bed.fn = "{bed_file}", bim.fn = "{bim_file}", fam.fn = "{fam_file}", out.gdsfn = "{gds_file}")
gds <- snpgdsOpen("{gds_file}")

grm_bin <- file("{grm_bin_file}", "rb")
grm_n_file_handle <- file("{grm_n_file}", "rb")
sample_ids <- read.table("{grm_id_file}", header=FALSE, stringsAsFactors=FALSE)
n_samples <- nrow(sample_ids)
n_values <- n_samples * (n_samples + 1) / 2
grm_values <- readBin(grm_bin, what="numeric", n=n_values, size=4)
close(grm_bin)
close(grm_n_file_handle)

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

pc_air <- pcair(gds = gds, kinobj = kin_matrix, kin.thresh = {kin_threshold}, divobj = NULL, div.thresh = {div_threshold}, num.cores = 1)

pcs <- pc_air$vectors[, 1:{n_pcs}, drop=FALSE]
colnames(pcs) <- paste0("PC", 1:{n_pcs})
output_df <- data.frame(IID = as.character(rownames(pcs)), pcs, stringsAsFactors = FALSE)

write.table(output_df, file = "{output_file}", quote = FALSE, row.names = FALSE, sep = "\\t")
write.table(data.frame(variance = pc_air$values[1:{n_pcs}]), file = "{variance_file}", quote = FALSE, row.names = FALSE, sep = "\\t")
write.table(data.frame(IID = as.character(pc_air$unrels)), file = "{unrelated_file}", quote = FALSE, row.names = FALSE, sep = "\\t")

snpgdsClose(gds)
"""
        
        r_script_file = f"{output_prefix}_pcair.R"
        with open(r_script_file, 'w') as f:
            f.write(r_script)
        
        if verbose:
            logger.info("Running PC-AiR in R...")
        
        result = subprocess.run(['Rscript', r_script_file], capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"PC-AiR failed:\n{result.stderr}\n{result.stdout}")
        
        pca_file = f"{output_prefix}_pcair.txt"
        if not os.path.exists(pca_file):
            raise RuntimeError(f"PC-AiR output file not found: {pca_file}")
        
        pca_df = pd.read_csv(pca_file, sep='\t')
        pca_df['IID'] = pca_df['IID'].astype(str)
        pca_df.set_index('IID', inplace=True)
        
        if verbose and os.path.exists(variance_file):
            variance = pd.read_csv(variance_file, sep='\t')
            var_explained = variance['variance'] / variance['variance'].sum()
            logger.info(f"Variance explained by first 5 PCs: {var_explained[:5].values}")
        
        return pca_df
        
    finally:
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
    drop_na: bool = False,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Attach principal components to phenotype DataFrame.
    
    Args:
        phenotype_df: Phenotype DataFrame (IID as index or column)
        pca_df: PCA DataFrame with IID as index and PC columns
        n_pcs: Number of PCs to attach (will use PC1 to PCn)
        pc_prefix: Prefix for PC column names (default: 'PC')
        sample_id_col: Column name in phenotype_df to use for matching with PCA IIDs
                      If None, uses phenotype_df.index
        drop_na: If True, remove samples with missing PCs after merging (default: False)
        verbose: Print information about merging
        
    Returns:
        Phenotype DataFrame with PC columns added
        
    Raises:
        ValueError: If requested PCs are not available in pca_df
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
    
    # Prepare PCA data - convert index to column
    pcs_to_add = pca_df[pc_cols].copy()
    pcs_to_add = pcs_to_add.reset_index()  # This creates 'IID' column from index
    pcs_to_add['IID'] = pcs_to_add['IID'].astype(str)
    
    # Prepare phenotype dataframe
    result_df = phenotype_df.copy()
    
    if sample_id_col is None:
        # Use index for matching
        result_df = result_df.reset_index()  # Convert index to column
        merge_col = result_df.columns[0]  # First column is the index
        result_df[merge_col] = result_df[merge_col].astype(str)
        use_index = True
        original_index_name = phenotype_df.index.name
    else:
        # Use specified column for matching
        if sample_id_col not in result_df.columns:
            raise ValueError(f"Column '{sample_id_col}' not found in phenotype_df")
        merge_col = sample_id_col
        result_df[merge_col] = result_df[merge_col].astype(str)
        use_index = False
        original_index_name = None
    
    # Merge PCs with phenotype data
    merged_df = result_df.merge(
        pcs_to_add,
        left_on=merge_col,
        right_on='IID',
        how='left',
        suffixes=('', '_pc')
    )
    
    # Drop the extra IID column from PCA
    if 'IID' in merged_df.columns and merge_col != 'IID':
        merged_df.drop('IID', axis=1, inplace=True)
    
    # Restore index if it was used
    if use_index:
        merged_df.set_index(merge_col, inplace=True)
        if original_index_name:
            merged_df.index.name = original_index_name
    
    # Remove samples with missing PCs if requested
    if drop_na:
        n_before = len(merged_df)
        merged_df = merged_df.dropna(subset=pc_cols)
        n_after = len(merged_df)
        if verbose and n_before > n_after:
            logger.info(f"Dropped {n_before - n_after} samples with missing PCs")
    
    if verbose:
        n_pheno_samples = len(phenotype_df)
        n_pca_samples = len(pca_df)
        n_with_pcs = merged_df[pc_cols[0]].notna().sum()
        
        logger.info(f"Attaching {n_pcs} PCs to phenotype data")
        logger.info(f"Phenotype samples: {n_pheno_samples}")
        logger.info(f"PCA samples: {n_pca_samples}")
        logger.info(f"Samples with PCs after merge: {n_with_pcs}")
        
        if not drop_na and n_with_pcs < n_pheno_samples:
            n_missing = n_pheno_samples - n_with_pcs
            logger.warning(f"{n_missing} samples in phenotype_df have no PCs (will have NA values)")
            logger.info("Set drop_na=True to remove samples with missing PCs")
    
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

def impute_covariates(
    phenotype_df: pd.DataFrame,
    covariate_cols: List[str],
    method: str = 'median',
    drop_na: bool = False,
    verbose: bool = True
) -> pd.DataFrame:
    """
    Impute missing values in covariates.
    
    Args:
        phenotype_df: Phenotype DataFrame with covariates
        covariate_cols: List of covariate column names to impute
        method: Imputation method. Options:
                - 'drop': Remove all rows with any missing values
                - 'mean': Replace with mean (numeric only)
                - 'median': Replace with median (numeric only)
                - 'mode': Replace with mode (works for categorical)
                - 'missforest': MissForest algorithm (requires missingpy)
                - 'knn': K-Nearest Neighbors (k=5)
                - 'mice': Multiple Imputation by Chained Equations (requires missingpy)
        drop_na: If True, drop rows with missing outcome after imputation (default: False)
        verbose: Print imputation information
        
    Returns:
        DataFrame with imputed covariates
        
    Note:
        For 'missforest' and 'mice', install: pip install missingpy
    """
    result_df = phenotype_df.copy()
    
    # Validate covariate columns exist
    missing_cols = [col for col in covariate_cols if col not in result_df.columns]
    if missing_cols:
        raise ValueError(f"Covariate columns not found: {missing_cols}")
    
    # Count missing values before imputation
    if verbose:
        n_total = len(result_df)
        missing_before = result_df[covariate_cols].isna().sum()
        logger.info(f"Missing values before imputation:")
        for col in covariate_cols:
            n_missing = missing_before[col]
            if n_missing > 0:
                logger.info(f"  {col}: {n_missing}/{n_total} ({n_missing/n_total*100:.1f}%)")
    
    if method == 'drop':
        # Drop rows with any missing covariates
        n_before = len(result_df)
        result_df = result_df.dropna(subset=covariate_cols)
        n_after = len(result_df)
        if verbose:
            logger.info(f"Dropped {n_before - n_after} rows with missing covariates")
    
    elif method in ['mean', 'median', 'mode']:
        # Simple imputation
        for col in covariate_cols:
            if result_df[col].isna().any():
                if method == 'mean':
                    fill_value = result_df[col].mean()
                elif method == 'median':
                    fill_value = result_df[col].median()
                elif method == 'mode':
                    fill_value = result_df[col].mode()[0] if len(result_df[col].mode()) > 0 else result_df[col].median()
                
                result_df[col].fillna(fill_value, inplace=True)
                if verbose:
                    logger.info(f"Imputed {col} with {method}: {fill_value}")
    
    elif method == 'knn':
        # KNN imputation
        from sklearn.impute import KNNImputer
        
        imputer = KNNImputer(n_neighbors=5)
        result_df[covariate_cols] = imputer.fit_transform(result_df[covariate_cols])
        if verbose:
            logger.info(f"Imputed covariates using KNN (k=5)")
    
    elif method == 'missforest':
        # MissForest imputation
        try:
            from missingpy import MissForest
        except ImportError:
            raise ImportError(
                "missingpy is required for MissForest imputation. "
                "Install with: pip install missingpy"
            )
        
        imputer = MissForest(max_iter=10, random_state=42)
        result_df[covariate_cols] = imputer.fit_transform(result_df[covariate_cols])
        if verbose:
            logger.info(f"Imputed covariates using MissForest")
    
    elif method == 'mice':
        # MICE imputation
        try:
            from missingpy import MICEImputer
        except ImportError:
            raise ImportError(
                "missingpy is required for MICE imputation. "
                "Install with: pip install missingpy"
            )
        
        imputer = MICEImputer(random_state=42)
        result_df[covariate_cols] = imputer.fit_transform(result_df[covariate_cols])
        if verbose:
            logger.info(f"Imputed covariates using MICE")
    
    else:
        raise ValueError(
            f"Unknown imputation method: {method}. "
            f"Options: 'drop', 'mean', 'median', 'mode', 'knn', 'missforest', 'mice'"
        )
    
    # Drop rows with missing outcome if requested
    if drop_na:
        n_before = len(result_df)
        result_df = result_df.dropna()
        n_after = len(result_df)
        if verbose and n_before > n_after:
            logger.info(f"Dropped {n_before - n_after} rows with any remaining missing values")
    
    # Count missing values after imputation
    if verbose:
        missing_after = result_df[covariate_cols].isna().sum()
        n_still_missing = missing_after.sum()
        if n_still_missing > 0:
            logger.warning(f"Still have {n_still_missing} missing values after imputation")
        else:
            logger.info(f"✓ All missing values imputed successfully")
    
    return result_df


def validate_and_align_data(
    genotype_df: pd.DataFrame,
    phenotype_df: pd.DataFrame,
    outcome_col: Optional[str] = None,
    covariate_cols: Optional[List[str]] = None,
    geno_id_col: Optional[str] = None,
    pheno_id_col: Optional[str] = None,
    keep_only_common: bool = True,
    verbose: bool = True
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Validate and align genotype and phenotype data by sample IDs.
    
    Args:
        genotype_df: Genotype DataFrame (samples x variants)
        phenotype_df: Phenotype DataFrame
        outcome_col: Name of outcome column (optional, for validation)
        covariate_cols: List of covariate columns (optional, for validation)
        geno_id_col: Column name for sample IDs in genotype_df (None = use index)
        pheno_id_col: Column name for sample IDs in phenotype_df (None = use index)
        keep_only_common: If True, keep only samples present in both datasets (default: True)
                         If False, raise error if samples don't match
        verbose: Print validation information
        
    Returns:
        Tuple of (aligned_genotype_df, aligned_phenotype_df)
        Both DataFrames will have matching samples in the same order
        
    Raises:
        ValueError: If no common samples found or if keep_only_common=False and samples don't match
    """
    if verbose:
        logger.info("Validating and aligning genotype and phenotype data...")
    
    # Store original index names
    geno_index_name = genotype_df.index.name
    pheno_index_name = phenotype_df.index.name
    
    # Get sample IDs from genotype
    if geno_id_col is None:
        geno_samples = genotype_df.index
    elif geno_id_col in genotype_df.columns:
        geno_samples = genotype_df[geno_id_col]
    else:
        raise ValueError(f"Column '{geno_id_col}' not found in genotype_df")
    
    # Get sample IDs from phenotype
    if pheno_id_col is None:
        pheno_samples = phenotype_df.index
    elif pheno_id_col in phenotype_df.columns:
        pheno_samples = phenotype_df[pheno_id_col]
    else:
        raise ValueError(f"Column '{pheno_id_col}' not found in phenotype_df")
    
    # Convert to strings for comparison
    geno_samples_str = pd.Index(geno_samples.astype(str))
    pheno_samples_str = pd.Index(pheno_samples.astype(str))
    
    # Check for duplicates
    if geno_samples_str.duplicated().any():
        n_dup = geno_samples_str.duplicated().sum()
        raise ValueError(f"Genotype data has {n_dup} duplicate sample IDs")
    
    if pheno_samples_str.duplicated().any():
        n_dup = pheno_samples_str.duplicated().sum()
        raise ValueError(f"Phenotype data has {n_dup} duplicate sample IDs")
    
    # Find common samples
    common_samples_str = geno_samples_str.intersection(pheno_samples_str)
    
    if len(common_samples_str) == 0:
        raise ValueError(
            f"No common samples found!\n"
            f"Genotype samples (first 5): {list(geno_samples[:5])}\n"
            f"Phenotype samples (first 5): {list(pheno_samples[:5])}\n"
            f"Check if sample IDs match between files."
        )
    
    # Report sample overlap
    n_geno = len(geno_samples)
    n_pheno = len(pheno_samples)
    n_common = len(common_samples_str)
    n_geno_only = n_geno - n_common
    n_pheno_only = n_pheno - n_common
    
    if verbose:
        logger.info(f"Sample overlap:")
        logger.info(f"  Genotype samples: {n_geno}")
        logger.info(f"  Phenotype samples: {n_pheno}")
        logger.info(f"  Common samples: {n_common}")
        if n_geno_only > 0:
            logger.info(f"  Genotype only: {n_geno_only}")
        if n_pheno_only > 0:
            logger.info(f"  Phenotype only: {n_pheno_only}")
    
    # Check if we should keep only common samples
    if not keep_only_common and (n_geno_only > 0 or n_pheno_only > 0):
        raise ValueError(
            f"Samples don't match perfectly (keep_only_common=False):\n"
            f"  {n_geno_only} samples only in genotype\n"
            f"  {n_pheno_only} samples only in phenotype\n"
            f"Set keep_only_common=True to keep only common samples."
        )
    
    # Create mapping from string to original values
    geno_str_to_orig = dict(zip(geno_samples_str, geno_samples))
    pheno_str_to_orig = dict(zip(pheno_samples_str, pheno_samples))
    
    # Get common samples in original types
    common_geno = pd.Index([geno_str_to_orig[s] for s in common_samples_str])
    common_pheno = pd.Index([pheno_str_to_orig[s] for s in common_samples_str])
    
    # Subset and align data
    if geno_id_col is None:
        geno_aligned = genotype_df.loc[common_geno].copy()
        # Preserve original index name
        geno_aligned.index.name = geno_index_name
    else:
        geno_aligned = genotype_df[genotype_df[geno_id_col].isin(common_geno)].copy()
        geno_aligned = geno_aligned.set_index(geno_id_col).loc[common_geno]
        # Set index name from column name
        geno_aligned.index.name = geno_id_col
    
    if pheno_id_col is None:
        pheno_aligned = phenotype_df.loc[common_pheno].copy()
        # Preserve original index name
        pheno_aligned.index.name = pheno_index_name
    else:
        pheno_aligned = phenotype_df[phenotype_df[pheno_id_col].isin(common_pheno)].copy()
        pheno_aligned = pheno_aligned.set_index(pheno_id_col).loc[common_pheno]
        # Set index name from column name
        pheno_aligned.index.name = pheno_id_col
    
    # Ensure indices match (convert pheno to match geno values AND preserve geno index name)
    pheno_aligned.index = geno_aligned.index.copy()
    pheno_aligned.index.name = geno_aligned.index.name  # FIXED: Preserve genotype index name
    
    # Validate outcome and covariates if specified
    if outcome_col is not None or covariate_cols is not None:
        required_cols = []
        if outcome_col:
            required_cols.append(outcome_col)
        if covariate_cols:
            required_cols.extend(covariate_cols)
        
        missing_cols = [col for col in required_cols if col not in pheno_aligned.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns in phenotype: {missing_cols}")
        
        # Check for missing values
        if verbose:
            for col in required_cols:
                n_missing = pheno_aligned[col].isna().sum()
                if n_missing > 0:
                    logger.warning(f"{col}: {n_missing} samples with missing values")
    
    if verbose:
        logger.info(f"✓ Data validated and aligned: {len(geno_aligned)} samples")
        logger.info(f"  Genotype shape: {geno_aligned.shape}")
        logger.info(f"  Genotype index: '{geno_aligned.index.name}'")
        logger.info(f"  Phenotype shape: {pheno_aligned.shape}")
        logger.info(f"  Phenotype index: '{pheno_aligned.index.name}'")
    
    return geno_aligned, pheno_aligned
    
    
# Backward compatibility: keep old function name
additive_gwas = standard_gwas






