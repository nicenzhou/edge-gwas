"""
Unit tests for EDGE GWAS core functionality.
"""

import pytest
import numpy as np
import pandas as pd
from edge_gwas import EDGEAnalysis


@pytest.fixture
def synthetic_data():
    """Generate synthetic data for testing."""
    np.random.seed(42)
    
    n_samples = 1000
    n_variants = 100
    
    # Generate genotypes (0, 1, 2)
    genotypes = np.random.choice([0, 1, 2], size=(n_samples, n_variants), 
                                p=[0.25, 0.5, 0.25])
    
    sample_ids = [f"SAMPLE_{i}" for i in range(n_samples)]
    variant_ids = [f"SNP_{i}" for i in range(n_variants)]
    
    genotype_df = pd.DataFrame(genotypes, index=sample_ids, columns=variant_ids)
    genotype_df.index.name = 'sample_id'
    
    # Generate phenotypes
    age = np.random.normal(50, 10, n_samples)
    sex = np.random.choice([0, 1], n_samples)
    pc1 = np.random.normal(0, 1, n_samples)
    pc2 = np.random.normal(0, 1, n_samples)
    
    # Binary outcome (affected by first SNP)
    risk_score = 0.1 * genotypes[:, 0] + 0.01 * age + 0.5 * sex + 0.1 * pc1
    prob = 1 / (1 + np.exp(-risk_score + 2))
    binary_outcome = np.random.binomial(1, prob)
    
    # Continuous outcome
    continuous_outcome = 25 + 0.5 * genotypes[:, 0] + 0.1 * age + sex + 0.2 * pc1 + np.random.normal(0, 2, n_samples)
    
    phenotype_df = pd.DataFrame({
        'sample_id': sample_ids,
        'binary_outcome': binary_outcome,
        'continuous_outcome': continuous_outcome,
        'age': age,
        'sex': sex,
        'PC1': pc1,
        'PC2': pc2
    })
    phenotype_df.set_index('sample_id', inplace=True)
    
    return genotype_df, phenotype_df


def test_edge_analysis_binary(synthetic_data):
    """Test EDGE analysis with binary outcome."""
    genotype_df, phenotype_df = synthetic_data
    
    # Split data
    train_geno = genotype_df.iloc[:500]
    test_geno = genotype_df.iloc[500:]
    train_pheno = phenotype_df.iloc[:500]
    test_pheno = phenotype_df.iloc[500:]
    
    # Initialize EDGE analysis
    edge = EDGEAnalysis(outcome_type='binary', n_jobs=1, verbose=False)
    
    # Calculate alpha
    alpha_df = edge.calculate_alpha(
        train_geno.iloc[:, :10],  # Use first 10 SNPs for speed
        train_pheno,
        outcome='binary_outcome',
        covariates=['age', 'sex', 'PC1', 'PC2']
    )
    
    assert len(alpha_df) > 0
    assert 'variant_id' in alpha_df.columns
    assert 'alpha_value' in alpha_df.columns
    
    # Apply alpha
    gwas_df = edge.apply_alpha(
        test_geno.iloc[:, :10],
        test_pheno,
        outcome='binary_outcome',
        covariates=['age', 'sex', 'PC1', 'PC2'],
        alpha_values=alpha_df
    )
    
    assert len(gwas_df) > 0
    assert 'variant_id' in gwas_df.columns
    assert 'pval' in gwas_df.columns


def test_edge_analysis_continuous(synthetic_data):
    """Test EDGE analysis with continuous outcome."""
    genotype_df, phenotype_df = synthetic_data
    
    # Split data
    train_geno = genotype_df.iloc[:500]
    test_geno = genotype_df.iloc[500:]
    train_pheno = phenotype_df.iloc[:500]
    test_pheno = phenotype_df.iloc[500:]
    
    # Initialize EDGE analysis
    edge = EDGEAnalysis(outcome_type='continuous', n_jobs=1, verbose=False)
    
    # Run full analysis
    alpha_df, gwas_df = edge.run_full_analysis(
        train_geno.iloc[:, :10],
        train_pheno,
        test_geno.iloc[:, :10],
        test_pheno,
        outcome='continuous_outcome',
        covariates=['age', 'sex', 'PC1', 'PC2']
    )
    
    assert len(alpha_df) > 0
    assert len(gwas_df) > 0
    assert all(col in gwas_df.columns for col in ['variant_id', 'coef', 'pval'])


def test_invalid_outcome_type():
    """Test that invalid outcome type raises error."""
    with pytest.raises(ValueError):
        edge = EDGEAnalysis(outcome_type='invalid')


def test_missing_alpha_values(synthetic_data):
    """Test that applying alpha without calculation raises error."""
    genotype_df, phenotype_df = synthetic_data
    
    edge = EDGEAnalysis(outcome_type='binary', verbose=False)
    
    with pytest.raises(ValueError):
        edge.apply_alpha(
            genotype_df.iloc[:100, :10],
            phenotype_df.iloc[:100],
            outcome='binary_outcome',
            covariates=['age', 'sex']
        )
