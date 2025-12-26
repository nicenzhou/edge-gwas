"""
Core EDGE GWAS analysis functions.
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm
from joblib import Parallel, delayed
import warnings
from statsmodels.tools.sm_exceptions import ConvergenceWarning
from typing import Optional, List, Tuple, Dict
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class EDGEAnalysis:
    """
    Main class for EDGE GWAS analysis.
    
    Attributes:
        outcome_type (str): Type of outcome - 'binary' or 'continuous'
        n_jobs (int): Number of parallel jobs for computation
        max_iter (int): Maximum iterations for model fitting
        verbose (bool): Whether to print progress messages
    """
    
    def __init__(
        self,
        outcome_type: str = 'binary',
        n_jobs: int = -1,
        max_iter: int = 1000,
        verbose: bool = True
    ):
        """
        Initialize EDGE Analysis.
        
        Args:
            outcome_type: 'binary' for logistic regression, 'continuous' for linear regression
            n_jobs: Number of parallel jobs (-1 uses all available cores)
            max_iter: Maximum iterations for model convergence
            verbose: Print progress information
        """
        if outcome_type not in ['binary', 'continuous']:
            raise ValueError("outcome_type must be 'binary' or 'continuous'")
        
        self.outcome_type = outcome_type
        self.n_jobs = n_jobs
        self.max_iter = max_iter
        self.verbose = verbose
        self.alpha_values = None
        self.skipped_snps = []
        
    def _fit_codominant_model(
        self,
        het_data: pd.Series,
        hom_data: pd.Series,
        phenotype_df: pd.DataFrame,
        outcome: str,
        covariates: List[str]
    ) -> pd.DataFrame:
        """
        Fit codominant model (separate effects for het and hom).
        
        Args:
            het_data: Heterozygous genotype indicator
            hom_data: Homozygous alternative genotype indicator
            phenotype_df: DataFrame containing outcome and covariates
            outcome: Name of outcome variable
            covariates: List of covariate names
            
        Returns:
            DataFrame with model results
        """
        # Merge genotype data
        data = pd.merge(
            het_data.to_frame(),
            hom_data.to_frame(),
            left_index=True,
            right_index=True,
            suffixes=('_het', '_hom')
        )
        
        # Merge with phenotype data
        merged_df = pd.merge(data, phenotype_df, left_index=True, right_index=True)
        merged_df = merged_df.dropna()
        
        snp_name = het_data.name
        
        # Prepare variables
        X = merged_df[[f'{snp_name}_het', f'{snp_name}_hom'] + covariates]
        y = merged_df[outcome]
        
        # Add constant
        X = sm.add_constant(X)
        
        # Fit model
        try:
            if self.outcome_type == 'binary':
                model = sm.Logit(y, X)
                result = model.fit(method='bfgs', maxiter=self.max_iter, disp=False)
            else:
                model = sm.OLS(y, X)
                result = model.fit()
                
            # Extract results
            result_df = pd.DataFrame({
                'snp': [snp_name],
                'coef_het': [result.params[f'{snp_name}_het']],
                'coef_hom': [result.params[f'{snp_name}_hom']],
                'std_err_het': [result.bse[f'{snp_name}_het']],
                'std_err_hom': [result.bse[f'{snp_name}_hom']],
                'stat_het': [result.tvalues[f'{snp_name}_het']],
                'stat_hom': [result.tvalues[f'{snp_name}_hom']],
                'pval_het': [result.pvalues[f'{snp_name}_het']],
                'pval_hom': [result.pvalues[f'{snp_name}_hom']],
                'conf_int_low_het': [result.conf_int().loc[f'{snp_name}_het', 0]],
                'conf_int_high_het': [result.conf_int().loc[f'{snp_name}_het', 1]],
                'conf_int_low_hom': [result.conf_int().loc[f'{snp_name}_hom', 0]],
                'conf_int_high_hom': [result.conf_int().loc[f'{snp_name}_hom', 1]],
                'n_samples': [len(y)]
            })
            
            return result_df
            
        except (np.linalg.LinAlgError, ConvergenceWarning) as e:
            logger.warning(f"Model fitting failed for {snp_name}: {str(e)}")
            self.skipped_snps.append(snp_name)
            return pd.DataFrame()
    
    def _fit_edge_model(
        self,
        edge_data: pd.Series,
        phenotype_df: pd.DataFrame,
        outcome: str,
        covariates: List[str]
    ) -> pd.DataFrame:
        """
        Fit EDGE-encoded model.
        
        Args:
            edge_data: EDGE-encoded genotype data
            phenotype_df: DataFrame containing outcome and covariates
            outcome: Name of outcome variable
            covariates: List of covariate names
            
        Returns:
            DataFrame with model results
        """
        # Merge with phenotype data
        merged_df = pd.merge(
            edge_data.to_frame(),
            phenotype_df,
            left_index=True,
            right_index=True
        )
        merged_df = merged_df.dropna()
        
        snp_name = edge_data.name
        
        # Prepare variables
        X = merged_df[[snp_name] + covariates]
        y = merged_df[outcome]
        
        # Add constant
        X = sm.add_constant(X)
        
        # Fit model
        try:
            if self.outcome_type == 'binary':
                model = sm.Logit(y, X)
                result = model.fit(method='bfgs', maxiter=self.max_iter, disp=False)
            else:
                model = sm.OLS(y, X)
                result = model.fit()
            
            # Extract results
            result_df = pd.DataFrame({
                'snp': [snp_name],
                'coef': [result.params[snp_name]],
                'std_err': [result.bse[snp_name]],
                'stat': [result.tvalues[snp_name]],
                'pval': [result.pvalues[snp_name]],
                'conf_int_low': [result.conf_int().loc[snp_name, 0]],
                'conf_int_high': [result.conf_int().loc[snp_name, 1]],
                'n_samples': [len(y)]
            })
            
            return result_df
            
        except (np.linalg.LinAlgError, ConvergenceWarning) as e:
            logger.warning(f"Model fitting failed for {snp_name}: {str(e)}")
            self.skipped_snps.append(snp_name)
            return pd.DataFrame()
    
    def calculate_alpha(
        self,
        genotype_data: pd.DataFrame,
        phenotype_df: pd.DataFrame,
        outcome: str,
        covariates: List[str],
        variant_info: Optional[pd.DataFrame] = None
    ) -> pd.DataFrame:
        """
        Calculate EDGE alpha values from training data.
        
        Args:
            genotype_data: DataFrame with genotypes (0, 1, 2 encoding)
                          Index should be sample IDs
            phenotype_df: DataFrame with outcome and covariates
                         Index should be sample IDs
            outcome: Name of outcome variable in phenotype_df
            covariates: List of covariate names in phenotype_df
            variant_info: Optional DataFrame with variant information
                         (columns: variant_id, ref_allele, alt_allele)
            
        Returns:
            DataFrame with alpha values for each variant
            Columns: variant_id, alpha_value, ref_allele, alt_allele, eaf, 
                    coef_het, coef_hom, convergence_status
        """
        self.skipped_snps = []
        alpha_results = []
        
        n_variants = genotype_data.shape[1]
        
        for idx, variant_id in enumerate(genotype_data.columns):
            if self.verbose and (idx + 1) % 100 == 0:
                logger.info(f"Processing variant {idx + 1}/{n_variants}")
            
            # Get genotype column
            geno = genotype_data[variant_id]
            
            # Create het and hom indicators
            het = geno.replace(2, 0)  # 0->0, 1->1, 2->0
            hom = geno.replace({0: 1, 1: 0, 2: 0})  # 0->1, 1->0, 2->0
            
            het.name = variant_id
            hom.name = variant_id
            
            # Fit codominant model
            result_df = self._fit_codominant_model(
                het, hom, phenotype_df, outcome, covariates
            )
            
            if result_df.empty:
                continue
            
            # Calculate alpha
            coef_het = result_df['coef_het'].iloc[0]
            coef_hom = result_df['coef_hom'].iloc[0]
            
            # Avoid division by zero
            if abs(coef_hom) < 1e-10:
                alpha_value = np.nan
                logger.warning(f"Alpha calculation failed for {variant_id}: coef_hom too close to zero")
            else:
                alpha_value = coef_het / coef_hom
            
            # Calculate EAF (effect allele frequency)
            valid_geno = geno.dropna()
            if len(valid_geno) > 0:
                eaf = (2 * (valid_geno == 0).sum() + (valid_geno == 1).sum()) / (2 * len(valid_geno))
            else:
                eaf = np.nan
            
            # Get variant info if provided
            if variant_info is not None and variant_id in variant_info.index:
                ref_allele = variant_info.loc[variant_id, 'ref_allele']
                alt_allele = variant_info.loc[variant_id, 'alt_allele']
            else:
                ref_allele = 'REF'
                alt_allele = 'ALT'
            
            alpha_results.append({
                'variant_id': variant_id,
                'alpha_value': alpha_value,
                'ref_allele': ref_allele,
                'alt_allele': alt_allele,
                'eaf': eaf,
                'coef_het': coef_het,
                'coef_hom': coef_hom,
                'std_err_het': result_df['std_err_het'].iloc[0],
                'std_err_hom': result_df['std_err_hom'].iloc[0],
                'pval_het': result_df['pval_het'].iloc[0],
                'pval_hom': result_df['pval_hom'].iloc[0],
                'n_samples': result_df['n_samples'].iloc[0],
                'convergence_status': 'converged'
            })
        
        self.alpha_values = pd.DataFrame(alpha_results)
        
        if self.verbose:
            logger.info(f"Alpha calculation complete. Processed {len(alpha_results)} variants.")
            logger.info(f"Skipped {len(self.skipped_snps)} variants due to convergence issues.")
        
        return self.alpha_values
    
    def apply_alpha(
        self,
        genotype_data: pd.DataFrame,
        phenotype_df: pd.DataFrame,
        outcome: str,
        covariates: List[str],
        alpha_values: Optional[pd.DataFrame] = None
    ) -> pd.DataFrame:
        """
        Apply EDGE alpha values to test data and perform GWAS.
        
        Args:
            genotype_data: DataFrame with genotypes (0, 1, 2 encoding)
                          Index should be sample IDs
            phenotype_df: DataFrame with outcome and covariates
                         Index should be sample IDs
            outcome: Name of outcome variable in phenotype_df
            covariates: List of covariate names in phenotype_df
            alpha_values: DataFrame with alpha values (from calculate_alpha)
                         If None, uses self.alpha_values from previous calculation
            
        Returns:
            DataFrame with GWAS results
            Columns: variant_id, coef, std_err, stat, pval, conf_int_low, 
                    conf_int_high, alpha_value, n_samples
        """
        # Use provided alpha values or stored values
        if alpha_values is None:
            if self.alpha_values is None:
                raise ValueError("No alpha values provided. Run calculate_alpha first or provide alpha_values.")
            alpha_values = self.alpha_values
        
        # Check that alpha_values has required columns
        required_cols = ['variant_id', 'alpha_value']
        if not all(col in alpha_values.columns for col in required_cols):
            raise ValueError(f"alpha_values must contain columns: {required_cols}")
        
        # Create alpha lookup dictionary
        alpha_dict = dict(zip(alpha_values['variant_id'], alpha_values['alpha_value']))
        
        self.skipped_snps = []
        gwas_results = []
        
        n_variants = genotype_data.shape[1]
        
        for idx, variant_id in enumerate(genotype_data.columns):
            if self.verbose and (idx + 1) % 100 == 0:
                logger.info(f"Processing variant {idx + 1}/{n_variants}")
            
            # Check if alpha value exists for this variant
            if variant_id not in alpha_dict:
                logger.warning(f"No alpha value found for {variant_id}, skipping.")
                self.skipped_snps.append(variant_id)
                continue
            
            alpha_value = alpha_dict[variant_id]
            
            # Skip if alpha is NaN
            if pd.isna(alpha_value):
                logger.warning(f"Alpha value is NaN for {variant_id}, skipping.")
                self.skipped_snps.append(variant_id)
                continue
            
            # Get genotype column
            geno = genotype_data[variant_id].copy()
            
            # Apply EDGE encoding: 0->1, 1->alpha, 2->0
            edge_encoded = geno.replace({0: 1.0, 1: alpha_value, 2: 0.0})
            edge_encoded.name = variant_id
            
            # Fit EDGE model
            result_df = self._fit_edge_model(
                edge_encoded, phenotype_df, outcome, covariates
            )
            
            if result_df.empty:
                continue
            
            # Add alpha value to results
            result_df['alpha_value'] = alpha_value
            result_df['variant_id'] = variant_id
            
            gwas_results.append(result_df)
        
        if not gwas_results:
            logger.warning("No variants were successfully analyzed.")
            return pd.DataFrame()
        
        gwas_df = pd.concat(gwas_results, ignore_index=True)
        
        if self.verbose:
            logger.info(f"EDGE GWAS complete. Analyzed {len(gwas_results)} variants.")
            logger.info(f"Skipped {len(self.skipped_snps)} variants.")
        
        return gwas_df
    
    def run_full_analysis(
        self,
        train_genotype: pd.DataFrame,
        train_phenotype: pd.DataFrame,
        test_genotype: pd.DataFrame,
        test_phenotype: pd.DataFrame,
        outcome: str,
        covariates: List[str],
        variant_info: Optional[pd.DataFrame] = None,
        output_prefix: Optional[str] = None
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Run complete EDGE analysis: calculate alpha on training data,
        apply alpha on test data.
        
        Args:
            train_genotype: Training genotype data
            train_phenotype: Training phenotype data
            test_genotype: Test genotype data
            test_phenotype: Test phenotype data
            outcome: Name of outcome variable
            covariates: List of covariate names
            variant_info: Optional variant information
            output_prefix: Optional prefix for output files
            
        Returns:
            Tuple of (alpha_df, gwas_df)
        """
        logger.info("Starting EDGE analysis...")
        
        # Calculate alpha values on training data
        logger.info("Step 1: Calculating alpha values on training data...")
        alpha_df = self.calculate_alpha(
            train_genotype,
            train_phenotype,
            outcome,
            covariates,
            variant_info
        )
        
        # Save alpha values if output prefix provided
        if output_prefix:
            alpha_file = f"{output_prefix}_alpha_values.csv"
            alpha_df.to_csv(alpha_file, index=False)
            logger.info(f"Alpha values saved to {alpha_file}")
        
        # Apply alpha values on test data
        logger.info("Step 2: Applying alpha values on test data...")
        gwas_df = self.apply_alpha(
            test_genotype,
            test_phenotype,
            outcome,
            covariates,
            alpha_df
        )
        
        # Save GWAS results if output prefix provided
        if output_prefix:
            gwas_file = f"{output_prefix}_gwas_results.csv"
            gwas_df.to_csv(gwas_file, index=False)
            logger.info(f"GWAS results saved to {gwas_file}")
            
            # Save skipped SNPs log
            if self.skipped_snps:
                log_file = f"{output_prefix}_skipped_snps.log"
                with open(log_file, 'w') as f:
                    f.write('\n'.join(self.skipped_snps))
                logger.info(f"Skipped SNPs log saved to {log_file}")
        
        logger.info("EDGE analysis complete!")
        
        return alpha_df, gwas_df
    
    def get_skipped_snps(self) -> List[str]:
        """
        Get list of SNPs that were skipped due to convergence issues.
        
        Returns:
            List of skipped SNP IDs
        """
        return self.skipped_snps
