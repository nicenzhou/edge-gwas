"""
Example: EDGE GWAS analysis for continuous outcome (BMI).
"""

import pandas as pd
import numpy as np
from edge_gwas import EDGEAnalysis
from edge_gwas.utils import (
    load_plink_data,
    prepare_phenotype_data,
    stratified_train_test_split,
    filter_variants_by_maf,
    filter_variants_by_missing
)
from edge_gwas.io_handlers import save_results, create_summary_report

# ========================================
# 1. Load genotype data
# ========================================
print("Loading genotype data...")
genotype_df, variant_info = load_plink_data(
    bed_file="ukb21007_c1_bmi_eur_hardcall8_qced.bed",
    bim_file="ukb21007_c1_bmi_eur_hardcall8_qced.bim",
    fam_file="ukb21007_c1_bmi_eur_hardcall8_qced.fam"
)

# ========================================
# 2. Filter variants
# ========================================
print("Filtering variants...")
genotype_df = filter_variants_by_maf(genotype_df, min_maf=0.01)
genotype_df = filter_variants_by_missing(genotype_df, max_missing=0.1)

# ========================================
# 3. Load and prepare phenotype data
# ========================================
print("Loading phenotype data...")
phenotype_df = prepare_phenotype_data(
    phenotype_file="train_df_bmi_eur.csv",
    outcome_col="BMI_Instance_0",
    covariate_cols=['age_baseline', 'sex', 'genotype_batch', 
                    'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 
                    'PC6', 'PC7', 'PC8', 'PC9', 'PC10'],
    sample_id_col='IID',
    sep=',',
    log_transform_outcome=True  # Apply log transformation
)

# Update outcome name after transformation
outcome_name = 'log10_BMI_Instance_0'

# ========================================
# 4. Split into training and test sets
# ========================================
print("Splitting data into training and test sets...")
train_geno, test_geno, train_pheno, test_pheno = stratified_train_test_split(
    genotype_df,
    phenotype_df,
    outcome_col=outcome_name,
    test_size=0.5,
    random_state=42,
    is_binary=False  # Continuous outcome
)

# ========================================
# 5. Run EDGE GWAS analysis
# ========================================
print("Running EDGE GWAS analysis...")
edge = EDGEAnalysis(
    outcome_type='continuous',
    n_jobs=-1,
    max_iter=1000,
    verbose=True
)

alpha_df, gwas_df = edge.run_full_analysis(
    train_genotype=train_geno,
    train_phenotype=train_pheno,
    test_genotype=test_geno,
    test_phenotype=test_pheno,
    outcome=outcome_name,
    covariates=['age_baseline', 'sex', 'genotype_batch', 
                'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 
                'PC6', 'PC7', 'PC8', 'PC9', 'PC10'],
    variant_info=variant_info,
    output_prefix='bmi_eur_chr1'
)

# ========================================
# 6. Generate summary report
# ========================================
print("Generating summary report...")
report = create_summary_report(
    gwas_df,
    alpha_df,
    significance_threshold=5e-8,
    output_file='bmi_eur_chr1_summary.txt'
)

print(report)

# ========================================
# 7. Display results
# ========================================
print("\nTop 10 variants by p-value:")
top_variants = gwas_df.nsmallest(10, 'pval')
print(top_variants[['variant_id', 'coef', 'pval', 'alpha_value']])

print("\nAlpha value statistics:")
print(alpha_df['alpha_value'].describe())

print("\nAnalysis complete!")
print(f"Results saved with prefix: bmi_eur_chr1")
