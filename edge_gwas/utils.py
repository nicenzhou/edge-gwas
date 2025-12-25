"""
Utility functions for EDGE GWAS analysis.
"""

import pandas as pd
import numpy as np
from pandas_plink import read_plink1_bin
from sklearn.model_selection import train_test_split, KFold
from typing import Tuple, List, Optional
import logging
from scipy import stats

logger = logging.getLogger(__name__)


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
        genotype_df: Genotype DataFrame (works with both hard calls and dosages)
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
        
        # Calculate allele frequency (works for both 0/1/2 and dosages)
        # For dosages, this gives the average dosage
        af = geno.mean() / 2
        
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
    
    genotype_df = genotype_df.loc[good_samples]
    phenotype_df = phenotype_df.loc[phenotype_df.index.intersection(good_samples)]
    
    n_after = genotype_df.shape[0]
    
    if verbose:
        logger.info(f"Filtered samples by call rate >= {min_call_rate}")
        logger.info(f"Kept {n_after}/{n_before} samples ({n_after/n_before*100:.1f}%)")
    
    return genotype_df, phenotype_df


def check_case_control_balance(
    phenotype_df: pd.DataFrame,
    outcome_col: str,
    verbose: bool = True
) -> dict:
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
    if verbose:
        logger.info("Calculating Hardy-Weinberg Equilibrium p-values...")
    
    hwe_pvals = {}
    
    for variant in genotype_df.columns:
        geno = genotype_df[variant].dropna()
        
        # Only use hard calls (0, 1, 2)
        geno = geno[geno.isin([0, 1, 2])]
        
        if len(geno) < 10:  # Skip if too few samples
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
        
        # Chi-square test
        if exp_aa > 0 and exp_ab > 0 and exp_bb > 0:
            obs = np.array([n_aa, n_ab, n_bb])
            exp = np.array([exp_aa, exp_ab, exp_bb])
            
            chi2_stat = np.sum((obs - exp) ** 2 / exp)
            pval = 1 - stats.chi2.cdf(chi2_stat, df=1)
            
            hwe_pvals[variant] = pval
        else:
            hwe_pvals[variant] = np.nan
    
    hwe_series = pd.Series(hwe_pvals)
    
    if verbose:
        n_valid = hwe_series.notna().sum()
        logger.info(f"Calculated HWE p-values for {n_valid} variants")
    
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
    
    # Filter variants
    variants_to_keep = hwe_pvals[hwe_pvals >= hwe_threshold].index
    genotype_df = genotype_df[variants_to_keep]
    
    n_after = genotype_df.shape[1]
    
    if verbose:
        logger.info(f"Filtered variants by HWE p-value >= {hwe_threshold}")
        logger.info(f"Kept {n_after}/{n_before} variants ({n_after/n_before*100:.1f}%)")
    
    return genotype_df


def additive_gwas(
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
    
    logger.info(f"Running additive GWAS ({outcome_type} outcome)...")
    
    results = []
    
    for variant in genotype_df.columns:
        try:
            # Prepare data
            X = phenotype_df[covariates].copy()
            X['genotype'] = genotype_df[variant]
            X = X.dropna()
            
            if len(X) < 10:  # Skip if too few samples
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
    
    logger.info(f"Completed additive GWAS for {len(results_df)} variants")
    
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
