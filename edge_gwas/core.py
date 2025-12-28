"""
Core EDGE GWAS analysis functions.
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import linalg
from scipy.stats import norm, rankdata
from scipy.optimize import minimize
from joblib import Parallel, delayed
import warnings
from statsmodels.tools.sm_exceptions import ConvergenceWarning
from typing import Optional, List, Tuple, Dict
import logging
import os

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
        outcome_transform (str): Transformation to apply to continuous outcomes
        ols_method (str): Optimization method for OLS regression
    """
    
    def __init__(
        self,
        outcome_type: str = 'binary',
        outcome_transform: Optional[str] = None,
        ols_method: str = 'bfgs',
        n_jobs: int = -1,
        max_iter: int = 1000,
        verbose: bool = True
    ):
        """
        Initialize EDGE Analysis.
        
        Args:
            outcome_type: 'binary' for logistic regression, 'continuous' for linear regression
            outcome_transform: Transformation for continuous outcomes:
                - None: No transformation
                - 'log': Natural log transformation
                - 'log10': Log base 10 transformation
                - 'inverse_normal': Inverse normal transformation (parametric)
                - 'rank_inverse_normal': Rank-based inverse normal transformation
            ols_method: Optimization method for OLS regression (default: 'bfgs'):
                - 'newton': Newton-Raphson algorithm
                - 'bfgs': Broyden-Fletcher-Goldfarb-Shannon algorithm (default)
                - 'lbfgs': Limited-memory BFGS algorithm
                - 'nm': Nelder-Mead algorithm
                - 'cg': Conjugate Gradient algorithm
                - 'ncg': Newton Conjugate Gradient algorithm
                - 'powell': Powell's conjugate direction algorithm
                - 'basinhopping': Basin-hopping algorithm
            n_jobs: Number of parallel jobs (-1 uses all available cores)
            max_iter: Maximum iterations for model convergence
            verbose: Print progress information
        """
        if outcome_type not in ['binary', 'continuous']:
            raise ValueError("outcome_type must be 'binary' or 'continuous'")
        
        valid_transforms = [None, 'log', 'log10', 'inverse_normal', 'rank_inverse_normal']
        if outcome_transform not in valid_transforms:
            raise ValueError(f"outcome_transform must be one of {valid_transforms}")
        
        if outcome_type == 'binary' and outcome_transform is not None:
            raise ValueError("outcome_transform can only be used with continuous outcomes")
        
        valid_ols_methods = ['newton', 'bfgs', 'lbfgs', 'nm', 'cg', 'ncg', 'powell', 'basinhopping']
        if ols_method not in valid_ols_methods:
            raise ValueError(f"ols_method must be one of {valid_ols_methods}")
        
        self.outcome_type = outcome_type
        self.outcome_transform = outcome_transform
        self.ols_method = ols_method
        self.n_jobs = n_jobs
        self.max_iter = max_iter
        self.verbose = verbose
        self.alpha_values = None
        self.skipped_snps = []
        
        if self.verbose:
            logger.info(f"OLS optimization method: {self.ols_method}")
        
    def _transform_outcome(self, y: pd.Series) -> pd.Series:
        """
        Apply transformation to continuous outcome.
        
        Args:
            y: Outcome series
            
        Returns:
            Transformed outcome series
        """
        if self.outcome_transform is None:
            return y
        
        y_transformed = y.copy()
        
        if self.outcome_transform == 'log':
            # Natural log transformation
            if (y <= 0).any():
                raise ValueError("Log transformation requires all positive values")
            y_transformed = np.log(y)
            if self.verbose:
                logger.info(f"Applied natural log transformation to outcome")
                
        elif self.outcome_transform == 'log10':
            # Log base 10 transformation
            if (y <= 0).any():
                raise ValueError("Log10 transformation requires all positive values")
            y_transformed = np.log10(y)
            if self.verbose:
                logger.info(f"Applied log10 transformation to outcome")
                
        elif self.outcome_transform == 'inverse_normal':
            # Parametric inverse normal transformation
            # Assumes data follows a normal distribution
            mean = y.mean()
            std = y.std()
            
            # Standardize
            z = (y - mean) / std
            
            # Apply inverse normal CDF
            y_transformed = pd.Series(
                norm.ppf((rankdata(z) - 0.5) / len(z)),
                index=y.index
            )
            if self.verbose:
                logger.info(f"Applied inverse normal transformation to outcome")
                
        elif self.outcome_transform == 'rank_inverse_normal':
            # Rank-based inverse normal transformation (RINT)
            # More robust to outliers
            n = len(y)
            
            # Get ranks (average for ties)
            ranks = rankdata(y, method='average')
            
            # Apply Blom's formula: (rank - 3/8) / (n + 1/4)
            # Alternative formulas:
            # - Van der Waerden: rank / (n + 1)
            # - Blom: (rank - 3/8) / (n + 1/4)
            # - Tukey: (rank - 1/3) / (n + 1/3)
            quantiles = (ranks - 3/8) / (n + 1/4)
            
            # Apply inverse normal CDF
            y_transformed = pd.Series(
                norm.ppf(quantiles),
                index=y.index
            )
            if self.verbose:
                logger.info(f"Applied rank-based inverse normal transformation to outcome")
        
        # Check for invalid values
        if y_transformed.isna().any() or np.isinf(y_transformed).any():
            n_invalid = y_transformed.isna().sum() + np.isinf(y_transformed).sum()
            logger.warning(f"Transformation produced {n_invalid} invalid values (NA or Inf)")
            y_transformed = y_transformed.replace([np.inf, -np.inf], np.nan)
        
        if self.verbose:
            logger.info(f"Outcome statistics after transformation:")
            logger.info(f"  Mean: {y_transformed.mean():.4f}")
            logger.info(f"  Std: {y_transformed.std():.4f}")
            logger.info(f"  Min: {y_transformed.min():.4f}")
            logger.info(f"  Max: {y_transformed.max():.4f}")
        
        return y_transformed
    
    def _prepare_grm_for_samples(
        self,
        sample_ids: pd.Index,
        grm_matrix: np.ndarray,
        grm_sample_ids: pd.DataFrame
    ) -> Tuple[np.ndarray, pd.Index]:
        """
        Extract and align GRM for samples in the analysis.
        
        Args:
            sample_ids: Sample IDs from the analysis data
            grm_matrix: Full GRM matrix from GCTA
            grm_sample_ids: DataFrame with FID and IID from GRM
            
        Returns:
            Tuple of (aligned_grm, common_sample_ids)
        """
        # Create sample ID mapping
        grm_sample_ids['sample_id'] = grm_sample_ids['IID'].astype(str)
        sample_ids_str = sample_ids.astype(str)
        
        # Find common samples maintaining order
        common_samples = [s for s in sample_ids_str if s in grm_sample_ids['sample_id'].values]
        
        if len(common_samples) == 0:
            raise ValueError("No common samples found between analysis data and GRM")
        
        if self.verbose:
            logger.info(f"Found {len(common_samples)} common samples between data and GRM")
        
        # Get indices for common samples in GRM
        grm_id_to_idx = {sid: idx for idx, sid in enumerate(grm_sample_ids['sample_id'])}
        grm_indices = [grm_id_to_idx[s] for s in common_samples]
        
        # Extract GRM submatrix for common samples
        aligned_grm = grm_matrix[np.ix_(grm_indices, grm_indices)]
        
        return aligned_grm, pd.Index(common_samples)
    
    def _transform_with_grm_linear(
        self,
        y: pd.Series,
        X: pd.DataFrame,
        grm: np.ndarray,
        h2: float = 0.5
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Transform phenotype and covariates using GRM for linear mixed model.
        
        This implements the transformation for a linear mixed model:
        y = Xβ + Zu + e, where u ~ N(0, σ²_g*GRM) and e ~ N(0, σ²_e*I)
        
        Args:
            y: Phenotype vector
            X: Design matrix (genotypes + covariates)
            grm: Genetic relationship matrix
            h2: Assumed heritability (default: 0.5)
            
        Returns:
            Tuple of (transformed_y, transformed_X)
        """
        n = len(y)
        
        # Construct V = h2*GRM + (1-h2)*I
        V = h2 * grm + (1 - h2) * np.eye(n)
        
        # Cholesky decomposition of V
        try:
            L = linalg.cholesky(V, lower=True)
            L_inv = linalg.solve_triangular(L, np.eye(n), lower=True)
            
            # Transform y and X: V^(-1/2) * y and V^(-1/2) * X
            y_transformed = L_inv @ y.values
            X_transformed = L_inv @ X.values
            
            return y_transformed, X_transformed
            
        except linalg.LinAlgError:
            logger.warning("GRM matrix is singular, using regularization")
            # Add small regularization term
            V_reg = V + 1e-6 * np.eye(n)
            L = linalg.cholesky(V_reg, lower=True)
            L_inv = linalg.solve_triangular(L, np.eye(n), lower=True)
            
            y_transformed = L_inv @ y.values
            X_transformed = L_inv @ X.values
            
            return y_transformed, X_transformed
    
    def _fit_ols_with_method(
        self,
        y: np.ndarray,
        X: np.ndarray
    ) -> sm.regression.linear_model.RegressionResultsWrapper:
        """
        Fit OLS model with specified optimization method.
        
        Args:
            y: Outcome vector
            X: Design matrix
            
        Returns:
            Fitted OLS model result
        """
        model = sm.OLS(y, X)
        
        # Fit using the specified optimization method
        result = model.fit(
            method=self.ols_method,
            maxiter=self.max_iter,
            disp=False
        )
        
        return result
    
    def _fit_logistic_mixed_model(
        self,
        y: np.ndarray,
        X: np.ndarray,
        grm: np.ndarray,
        tau: float = 1.0
    ) -> Dict:
        """
        Fit logistic mixed model with GRM using penalized quasi-likelihood.
        
        This implements a simplified GMMAT approach for binary outcomes.
        
        Args:
            y: Binary outcome vector
            X: Design matrix (including intercept)
            grm: Genetic relationship matrix
            tau: Variance component ratio (default: 1.0)
            
        Returns:
            Dictionary with coefficients, standard errors, and p-values
        """
        n = X.shape[0]
        p = X.shape[1]
        
        # Initialize coefficients with standard logistic regression
        try:
            init_model = sm.Logit(y, X)
            init_result = init_model.fit(disp=False, maxiter=100)
            beta = init_result.params.values
        except:
            beta = np.zeros(p)
            beta[0] = np.log(y.mean() / (1 - y.mean() + 1e-10))
        
        # Construct covariance matrix
        # V = tau * GRM + I (on the logit scale, this is approximate)
        Sigma = tau * grm + np.eye(n)
        
        try:
            # Cholesky decomposition
            L = linalg.cholesky(Sigma, lower=True)
            Sigma_inv = linalg.cho_solve((L, True), np.eye(n))
        except:
            # Fallback to pseudo-inverse if singular
            logger.warning("Singular GRM, using pseudo-inverse")
            Sigma_inv = linalg.pinv(Sigma)
        
        # Iteratively reweighted least squares with penalty
        for iteration in range(self.max_iter):
            # Calculate fitted values
            eta = X @ beta
            mu = 1 / (1 + np.exp(-eta))
            mu = np.clip(mu, 1e-10, 1 - 1e-10)
            
            # Working weights
            W = mu * (1 - mu)
            W = np.clip(W, 1e-10, None)
            W_mat = np.diag(W)
            
            # Working response
            z = eta + (y - mu) / W
            
            # Update beta with penalty
            # (X'WX + Sigma_inv)^{-1} X'Wz
            try:
                XtWX = X.T @ W_mat @ X
                XtWz = X.T @ (W * z)
                
                # Add penalty term (simplified)
                penalty_strength = 0.01  # Small penalty for stability
                A = XtWX + penalty_strength * Sigma_inv[:p, :p]
                
                beta_new = linalg.solve(A, XtWz, assume_a='pos')
                
                # Check convergence
                if np.max(np.abs(beta_new - beta)) < 1e-6:
                    beta = beta_new
                    break
                
                beta = beta_new
                
            except linalg.LinAlgError:
                logger.warning("Matrix inversion failed in logistic mixed model")
                break
        
        # Calculate standard errors
        try:
            # Information matrix
            mu = 1 / (1 + np.exp(-X @ beta))
            mu = np.clip(mu, 1e-10, 1 - 1e-10)
            W = mu * (1 - mu)
            W_mat = np.diag(W)
            
            XtWX = X.T @ W_mat @ X
            vcov = linalg.inv(XtWX)
            
            se = np.sqrt(np.diag(vcov))
            
            # Test statistics
            z_stats = beta / se
            pvals = 2 * (1 - norm.cdf(np.abs(z_stats)))
            
            # Confidence intervals
            ci_lower = beta - 1.96 * se
            ci_upper = beta + 1.96 * se
            
        except:
            logger.warning("Could not calculate standard errors for logistic mixed model")
            se = np.full_like(beta, np.nan)
            z_stats = np.full_like(beta, np.nan)
            pvals = np.full_like(beta, np.nan)
            ci_lower = np.full_like(beta, np.nan)
            ci_upper = np.full_like(beta, np.nan)
        
        return {
            'params': beta,
            'bse': se,
            'tvalues': z_stats,
            'pvalues': pvals,
            'conf_int_lower': ci_lower,
            'conf_int_upper': ci_upper
        }
    
    def _fit_codominant_model(
        self,
        het_data: pd.Series,
        hom_data: pd.Series,
        phenotype_df: pd.DataFrame,
        outcome: str,
        covariates: List[str],
        grm: Optional[np.ndarray] = None,
        grm_sample_ids: Optional[pd.Index] = None,
        mean_centered: bool = False
    ) -> pd.DataFrame:
        """
        Fit codominant model (separate effects for het and hom).
        
        Args:
            mean_centered: If True, fit codominant regression without intercept.
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
        
        # If GRM is provided, subset to common samples
        if grm is not None and grm_sample_ids is not None:
            merged_df = merged_df.loc[merged_df.index.intersection(grm_sample_ids)]
            if len(merged_df) == 0:
                logger.warning(f"No samples remain after GRM alignment for {het_data.name}")
                return pd.DataFrame()
        
        snp_name = het_data.name
        
        # Prepare variables
        X = merged_df[[f'{snp_name}_het', f'{snp_name}_hom'] + covariates]
        y = merged_df[outcome]
        
        # Apply outcome transformation if specified (for continuous outcomes)
        if self.outcome_type == 'continuous' and self.outcome_transform is not None:
            try:
                y = self._transform_outcome(y)
            except Exception as e:
                logger.warning(f"Outcome transformation failed for {snp_name}: {str(e)}")
                self.skipped_snps.append(snp_name)
                return pd.DataFrame()
        
        # Add constant ONLY if NOT mean-centered
        if not mean_centered:
            X = sm.add_constant(X)
        
        # Apply GRM if provided
        if grm is not None and grm_sample_ids is not None:
            # Align GRM to current samples
            sample_indices = [list(grm_sample_ids).index(s) for s in merged_df.index]
            aligned_grm = grm[np.ix_(sample_indices, sample_indices)]
            
            if self.outcome_type == 'continuous' or mean_centered:
                # Transform data for linear mixed model
                # (Use linear model for mean-centered binary outcome too)
                try:
                    y_transformed, X_transformed = self._transform_with_grm_linear(y, X, aligned_grm)
                    
                    # Fit OLS on transformed data
                    model = sm.OLS(y_transformed, X_transformed)
                    result = model.fit()
                except Exception as e:
                    logger.warning(f"GRM-based linear model fitting failed for {snp_name}: {str(e)}")
                    self.skipped_snps.append(snp_name)
                    return pd.DataFrame()
                    
            else:  # binary outcome without mean-centering
                # Fit logistic mixed model
                try:
                    result_dict = self._fit_logistic_mixed_model(
                        y.values, X.values, aligned_grm
                    )
                    
                    # Create a result-like object
                    class MixedModelResult:
                        def __init__(self, res_dict, feature_names):
                            self.params = pd.Series(res_dict['params'], index=feature_names)
                            self.bse = pd.Series(res_dict['bse'], index=feature_names)
                            self.tvalues = pd.Series(res_dict['tvalues'], index=feature_names)
                            self.pvalues = pd.Series(res_dict['pvalues'], index=feature_names)
                            self._conf_int_lower = pd.Series(res_dict['conf_int_lower'], index=feature_names)
                            self._conf_int_upper = pd.Series(res_dict['conf_int_upper'], index=feature_names)
                        
                        def conf_int(self):
                            return pd.DataFrame({
                                0: self._conf_int_lower,
                                1: self._conf_int_upper
                            })
                    
                    result = MixedModelResult(result_dict, X.columns)
                    
                except Exception as e:
                    logger.warning(f"GRM-based logistic model fitting failed for {snp_name}: {str(e)}")
                    self.skipped_snps.append(snp_name)
                    return pd.DataFrame()
        else:
            # Fit standard model without GRM
            try:
                if self.outcome_type == 'binary' and not mean_centered:
                    # Binary outcome: use logistic regression with optimization
                    model = sm.Logit(y, X)
                    result = model.fit(method='bfgs', maxiter=self.max_iter, disp=False)
                else:
                    # Continuous outcome or mean-centered: use OLS
                    model = sm.OLS(y, X)
                    result = model.fit()
            except Exception as e:
                logger.warning(f"Model fitting failed for {snp_name}: {str(e)}")
                self.skipped_snps.append(snp_name)
                return pd.DataFrame()
        
        # Extract results
        try:
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
            
        except Exception as e:
            logger.warning(f"Result extraction failed for {snp_name}: {str(e)}")
            self.skipped_snps.append(snp_name)
            return pd.DataFrame()


    
    def _fit_edge_model(
        self,
        edge_data: pd.Series,
        phenotype_df: pd.DataFrame,
        outcome: str,
        covariates: List[str],
        grm: Optional[np.ndarray] = None,
        grm_sample_ids: Optional[pd.Index] = None
    ) -> pd.DataFrame:
        """
        Fit EDGE-encoded model.
        
        Args:
            edge_data: EDGE-encoded genotype data
            phenotype_df: DataFrame containing outcome and covariates
            outcome: Name of outcome variable
            covariates: List of covariate names
            grm: Optional aligned GRM matrix for mixed model
            grm_sample_ids: Sample IDs corresponding to GRM rows
            
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
        
        # If GRM is provided, subset to common samples
        if grm is not None and grm_sample_ids is not None:
            merged_df = merged_df.loc[merged_df.index.intersection(grm_sample_ids)]
            if len(merged_df) == 0:
                logger.warning(f"No samples remain after GRM alignment for {edge_data.name}")
                return pd.DataFrame()
        
        snp_name = edge_data.name
        
        # Prepare variables
        X = merged_df[[snp_name] + covariates]
        y = merged_df[outcome]
        
        # Apply outcome transformation if specified (for continuous outcomes)
        if self.outcome_type == 'continuous' and self.outcome_transform is not None:
            try:
                y = self._transform_outcome(y)
            except Exception as e:
                logger.warning(f"Outcome transformation failed for {snp_name}: {str(e)}")
                self.skipped_snps.append(snp_name)
                return pd.DataFrame()
        
        # Add constant
        X = sm.add_constant(X)
        
        # Apply GRM if provided
        if grm is not None and grm_sample_ids is not None:
            # Align GRM to current samples
            sample_indices = [list(grm_sample_ids).index(s) for s in merged_df.index]
            aligned_grm = grm[np.ix_(sample_indices, sample_indices)]
            
            if self.outcome_type == 'continuous':
                # Transform data for linear mixed model
                try:
                    y_transformed, X_transformed = self._transform_with_grm_linear(y, X, aligned_grm)
                    
                    # Fit OLS on transformed data - NO OPTIMIZATION METHOD NEEDED
                    model = sm.OLS(y_transformed, X_transformed)
                    result = model.fit()
                except Exception as e:
                    logger.warning(f"GRM-based linear model fitting failed for {snp_name}: {str(e)}")
                    self.skipped_snps.append(snp_name)
                    return pd.DataFrame()
                    
            else:  # binary outcome
                # Fit logistic mixed model
                try:
                    result_dict = self._fit_logistic_mixed_model(
                        y.values, X.values, aligned_grm
                    )
                    
                    # Create a result-like object
                    class MixedModelResult:
                        def __init__(self, res_dict, feature_names):
                            self.params = pd.Series(res_dict['params'], index=feature_names)
                            self.bse = pd.Series(res_dict['bse'], index=feature_names)
                            self.tvalues = pd.Series(res_dict['tvalues'], index=feature_names)
                            self.pvalues = pd.Series(res_dict['pvalues'], index=feature_names)
                            self._conf_int_lower = pd.Series(res_dict['conf_int_lower'], index=feature_names)
                            self._conf_int_upper = pd.Series(res_dict['conf_int_upper'], index=feature_names)
                        
                        def conf_int(self):
                            return pd.DataFrame({
                                0: self._conf_int_lower,
                                1: self._conf_int_upper
                            })
                    
                    result = MixedModelResult(result_dict, X.columns)
                    
                except Exception as e:
                    logger.warning(f"GRM-based logistic model fitting failed for {snp_name}: {str(e)}")
                    self.skipped_snps.append(snp_name)
                    return pd.DataFrame()
        else:
            # Fit standard model without GRM
            try:
                if self.outcome_type == 'binary':
                    # Binary outcome: use logistic regression with optimization
                    model = sm.Logit(y, X)
                    result = model.fit(method='bfgs', maxiter=self.max_iter, disp=False)
                else:
                    # Continuous outcome: use OLS with direct solution (NO optimization method)
                    model = sm.OLS(y, X)
                    result = model.fit()
            except Exception as e:
                logger.warning(f"Model fitting failed for {snp_name}: {str(e)}")
                self.skipped_snps.append(snp_name)
                return pd.DataFrame()
        
        # Extract results
        try:
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
            
        except Exception as e:
            logger.warning(f"Result extraction failed for {snp_name}: {str(e)}")
            self.skipped_snps.append(snp_name)
            return pd.DataFrame()
    
    def calculate_alpha(
        self,
        genotype_data: pd.DataFrame,
        phenotype_df: pd.DataFrame,
        outcome: str,
        covariates: List[str],
        variant_info: Optional[pd.DataFrame] = None,
        grm_matrix: Optional[np.ndarray] = None,
        grm_sample_ids: Optional[pd.DataFrame] = None,
        mean_centered: bool = False
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
                         Index: variant_id
                         Columns: chrom, pos, ref_allele, alt_allele, MAF
            grm_matrix: Optional GRM matrix from GCTA (for population structure control)
            grm_sample_ids: DataFrame with FID and ID corresponding to GRM rows
            mean_centered: If True, use mean-centered codominant model without intercept
                          (default: False)
            
        Returns:
            DataFrame with alpha values for each variant
            Columns: chrom, pos, variant_id, alpha_value, ref_allele, alt_allele, 
                    eaf, coef_het, coef_hom, std_err_het, std_err_hom, 
                    pval_het, pval_hom, n_samples, convergence_status
        """
        self.skipped_snps = []
        alpha_results = []
        
        # Prepare GRM if provided
        aligned_grm = None
        grm_ids = None
        if grm_matrix is not None and grm_sample_ids is not None:
            if self.verbose:
                logger.info(f"Incorporating GRM for population structure control")
            aligned_grm, grm_ids = self._prepare_grm_for_samples(
                genotype_data.index,
                grm_matrix,
                grm_sample_ids
            )
        
        if self.verbose and mean_centered:
            logger.info("Using mean-centered codominant model (no intercept)")
        
        n_variants = genotype_data.shape[1]
        
        for idx, variant_id in enumerate(genotype_data.columns):
            if self.verbose and (idx + 1) % 100 == 0:
                logger.info(f"Processing variant {idx + 1}/{n_variants}")
            
            # Get genotype column
            geno = genotype_data[variant_id]
            
            # Create het and hom indicators
            het = (geno == 1).astype(int)
            hom = (geno == 2).astype(int)
            
            het.name = variant_id
            hom.name = variant_id
            
            # Fit codominant model with optional GRM and mean-centering
            result_df = self._fit_codominant_model(
                het, hom, phenotype_df, outcome, covariates, 
                aligned_grm, grm_ids, mean_centered
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
            
            alpha_results.append({
                'variant_id': variant_id,
                'alpha_value': alpha_value,
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
        
        alpha_df = pd.DataFrame(alpha_results)
        
        # Merge with variant_info if provided
        if variant_info is not None:
            # Ensure variant_info index name is variant_id
            if variant_info.index.name != 'variant_id':
                variant_info = variant_info.copy()
                variant_info.index.name = 'variant_id'
            
            # Reset index of alpha_df to merge on variant_id
            alpha_df = alpha_df.set_index('variant_id')
            
            # Select relevant columns from variant_info
            info_cols = []
            if 'chrom' in variant_info.columns:
                info_cols.append('chrom')
            if 'pos' in variant_info.columns:
                info_cols.append('pos')
            if 'ref_allele' in variant_info.columns:
                info_cols.append('ref_allele')
            if 'alt_allele' in variant_info.columns:
                info_cols.append('alt_allele')
            
            if info_cols:
                alpha_df = alpha_df.join(variant_info[info_cols], how='left')
            
            # Reset index to get variant_id as column
            alpha_df = alpha_df.reset_index()
        
        # Reorder columns
        final_cols = ['chrom', 'pos', 'variant_id', 'alpha_value', 'ref_allele', 'alt_allele',
                      'eaf', 'coef_het', 'coef_hom', 'std_err_het', 'std_err_hom',
                      'pval_het', 'pval_hom', 'n_samples', 'convergence_status']
        
        # Only include columns that exist
        final_cols = [col for col in final_cols if col in alpha_df.columns]
        alpha_df = alpha_df[final_cols]
        
        self.alpha_values = alpha_df
        
        if self.verbose:
            logger.info(f"Alpha calculation complete. Processed {len(alpha_df)} variants.")
            logger.info(f"Skipped {len(self.skipped_snps)} variants due to convergence issues.")
            if mean_centered:
                logger.info("Mean-centered codominant model was used (no intercept)")
            if self.outcome_transform:
                logger.info(f"Outcome transformation applied: {self.outcome_transform}")
            logger.info(f"OLS optimization method used: {self.ols_method}")
        
        return self.alpha_values

    
    def apply_alpha(
        self,
        genotype_data: pd.DataFrame,
        phenotype_df: pd.DataFrame,
        outcome: str,
        covariates: List[str],
        alpha_values: Optional[pd.DataFrame] = None,
        grm_matrix: Optional[np.ndarray] = None,
        grm_sample_ids: Optional[pd.DataFrame] = None,
        variant_info: Optional[pd.DataFrame] = None
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
            grm_matrix: Optional GRM matrix from GCTA (for population structure control)
            grm_sample_ids: DataFrame with FID and IID corresponding to GRM rows
            variant_info: Optional DataFrame with variant information
                         Index: variant_id
                         Columns: chrom, pos, ref_allele, alt_allele, MAF
            
        Returns:
            DataFrame with GWAS results
            Columns: chrom, pos, ref, alt, variant_id, snp, alpha_value, 
                    coef, std_err, stat, pval, conf_int_low, conf_int_high, 
                    n_samples, maf, eaf
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
        
        # Prepare GRM if provided
        aligned_grm = None
        grm_ids = None
        if grm_matrix is not None and grm_sample_ids is not None:
            if self.verbose:
                logger.info(f"Incorporating GRM for population structure control")
            aligned_grm, grm_ids = self._prepare_grm_for_samples(
                genotype_data.index,
                grm_matrix,
                grm_sample_ids
            )
        
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
            
            # Apply EDGE encoding: 0->0, 1->alpha, 2->1
            edge_encoded = geno.replace({0: 0.0, 1: alpha_value, 2: 1.0})
            edge_encoded.name = variant_id
            
            # Fit EDGE model with optional GRM
            result_df = self._fit_edge_model(
                edge_encoded, phenotype_df, outcome, covariates, aligned_grm, grm_ids
            )
            
            if result_df.empty:
                continue
            
            # Calculate frequencies
            valid_geno = geno.dropna()
            if len(valid_geno) > 0:
                # EAF (Effect Allele Frequency) = ALT allele frequency
                eaf = (2 * (valid_geno == 2).sum() + (valid_geno == 1).sum()) / (2 * len(valid_geno))
                # MAF (Minor Allele Frequency)
                maf = min(eaf, 1 - eaf)
            else:
                eaf = np.nan
                maf = np.nan
            
            gwas_results.append({
                'variant_id': variant_id,
                'snp': variant_id,
                'alpha_value': alpha_value,
                'coef': result_df['coef'].iloc[0],
                'std_err': result_df['std_err'].iloc[0],
                'stat': result_df['stat'].iloc[0],
                'pval': result_df['pval'].iloc[0],
                'conf_int_low': result_df['conf_int_low'].iloc[0],
                'conf_int_high': result_df['conf_int_high'].iloc[0],
                'n_samples': result_df['n_samples'].iloc[0],
                'maf': maf,
                'eaf': eaf
            })
        
        if not gwas_results:
            logger.warning("No variants were successfully analyzed.")
            return pd.DataFrame()
        
        gwas_df = pd.DataFrame(gwas_results)
        
        # Merge with variant_info if provided
        if variant_info is not None:
            # Ensure variant_info index name is variant_id
            if variant_info.index.name != 'variant_id':
                variant_info = variant_info.copy()
                variant_info.index.name = 'variant_id'
            
            # Set index of gwas_df to merge on variant_id
            gwas_df = gwas_df.set_index('variant_id')
            
            # Select relevant columns from variant_info
            info_cols = []
            if 'chrom' in variant_info.columns:
                info_cols.append('chrom')
            if 'pos' in variant_info.columns:
                info_cols.append('pos')
            if 'ref_allele' in variant_info.columns:
                info_cols.append('ref_allele')
            if 'alt_allele' in variant_info.columns:
                info_cols.append('alt_allele')
            
            if info_cols:
                gwas_df = gwas_df.join(variant_info[info_cols], how='left')
            
            # Reset index to get variant_id as column
            gwas_df = gwas_df.reset_index()
        
        # Rename ref_allele and alt_allele to ref and alt for final output
        if 'ref_allele' in gwas_df.columns:
            gwas_df = gwas_df.rename(columns={'ref_allele': 'ref'})
        if 'alt_allele' in gwas_df.columns:
            gwas_df = gwas_df.rename(columns={'alt_allele': 'alt'})
        
        # Reorder columns
        final_cols = ['chrom', 'pos', 'ref', 'alt', 'variant_id', 'snp', 'alpha_value',
                      'coef', 'std_err', 'stat', 'pval', 'conf_int_low', 'conf_int_high',
                      'n_samples', 'maf', 'eaf']
        
        # Only include columns that exist
        final_cols = [col for col in final_cols if col in gwas_df.columns]
        gwas_df = gwas_df[final_cols]
        
        if self.verbose:
            logger.info(f"EDGE GWAS complete. Analyzed {len(gwas_df)} variants.")
            logger.info(f"Skipped {len(self.skipped_snps)} variants.")
            if self.outcome_transform:
                logger.info(f"Outcome transformation applied: {self.outcome_transform}")
            logger.info(f"OLS optimization method used: {self.ols_method}")
            
            # Report MAF/EAF statistics
            if 'maf' in gwas_df.columns:
                valid_maf = gwas_df['maf'].dropna()
                if len(valid_maf) > 0:
                    logger.info(f"MAF range: {valid_maf.min():.4f} - {valid_maf.max():.4f}")
        
        return gwas_df


    def _parse_variant_id(
        self,
        variant_id: str,
        pattern: str
    ) -> Tuple[str, int, str, str]:
        """
        Parse variant ID to extract chromosome, position, ref, and alt alleles.
        
        Args:
            variant_id: Variant identifier string
            pattern: Pattern describing the variant_id format
            
        Returns:
            Tuple of (chrom, pos, ref_allele, alt_allele)
            
        Raises:
            ValueError: If variant_id cannot be parsed with the given pattern
        """
        parts = None
        chrom = None
        pos = None
        ref_allele = None
        alt_allele = None
        
        if pattern == 'chr:pos':
            # Format: 1:12345
            parts = variant_id.split(':')
            if len(parts) >= 2:
                chrom = parts[0]
                pos = int(parts[1])
                ref_allele = None
                alt_allele = None
        
        elif pattern == 'chr:pos:ref:alt':
            # Format: 1:12345:A:T
            parts = variant_id.split(':')
            if len(parts) >= 4:
                chrom = parts[0]
                pos = int(parts[1])
                ref_allele = parts[2]
                alt_allele = parts[3]
        
        elif pattern == 'chr:pos_ref_alt':
            # Format: 1:12345_A_T
            if ':' in variant_id and '_' in variant_id:
                chr_pos, alleles = variant_id.split(':', 1)
                chrom = chr_pos
                if '_' in alleles:
                    parts = alleles.split('_')
                    if len(parts) >= 3:
                        pos = int(parts[0])
                        ref_allele = parts[1]
                        alt_allele = parts[2]
        
        elif pattern == 'chr_pos_ref_alt':
            # Format: 1_12345_A_T
            parts = variant_id.split('_')
            if len(parts) >= 4:
                chrom = parts[0]
                pos = int(parts[1])
                ref_allele = parts[2]
                alt_allele = parts[3]
        
        else:
            raise ValueError(f"Unknown variant_id pattern: {pattern}")
        
        if chrom is None or pos is None:
            raise ValueError(f"Could not parse variant_id '{variant_id}' with pattern '{pattern}'")
        
        return chrom, pos, ref_allele, alt_allele
    
    def run_full_analysis(
        self,
        train_genotype: pd.DataFrame,
        train_phenotype: pd.DataFrame,
        test_genotype: pd.DataFrame,
        test_phenotype: pd.DataFrame,
        outcome: str,
        covariates: List[str],
        variant_info: Optional[pd.DataFrame] = None,
        grm_matrix: Optional[np.ndarray] = None,
        grm_sample_ids: Optional[pd.DataFrame] = None,
        mean_centered: bool = False,
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
            grm_matrix: Optional GRM matrix from GCTA
            grm_sample_ids: Optional sample IDs for GRM
            mean_centered: If True, use mean-centered model without intercept (default: False)
            output_prefix: Optional prefix for output files
            
        Returns:
            Tuple of (alpha_df, gwas_df)
        """
        logger.info("Starting EDGE analysis...")
        
        if self.outcome_transform:
            logger.info(f"Outcome transformation: {self.outcome_transform}")
        
        if mean_centered:
            logger.info("Using mean-centered codominant model (no intercept)")
        
        logger.info(f"OLS optimization method: {self.ols_method}")
        
        # Calculate alpha values on training data
        logger.info("Step 1: Calculating alpha values on training data...")
        alpha_df = self.calculate_alpha(
            train_genotype,
            train_phenotype,
            outcome,
            covariates,
            variant_info,
            grm_matrix,
            grm_sample_ids,
            mean_centered
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
            alpha_df,
            grm_matrix,
            grm_sample_ids
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

