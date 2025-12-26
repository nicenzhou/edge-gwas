.. _futureupdates:

Advanced Topics for Further Updates
===================

This guide covers potential troubleshooting walkthroughs in edge-gwas.

Extensions and Future Directions
---------------------------------

Planned Features
~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Gene-based testing**: Aggregate SNPs within genes
2. **Mixed models**: Account for relatedness and population structure
3. **X chromosome**: Proper handling of X-linked inheritance
4. **Rare variant analysis**: Optimized methods for MAF < 0.01
5. **Meta-analysis tools**: Combine α estimates across studies
6. **GPU acceleration**: Faster computation for biobank-scale data
7. **Additional transformations**: Box-Cox, Yeo-Johnson
8. **Adaptive method selection**: Automatically choose best optimization method
9. **Interactive optimization**: Real-time convergence monitoring
10. **Batch effect correction**: Built-in adjustments for technical covariates

Research Directions
~~~~~~~~~~~~~~~~~~~

* **Optimal splitting ratio**: Beyond 50/50 for different scenarios
* **Cross-validation**: Stability of α across folds
* **Bayesian EDGE**: Prior distributions on α
* **Multi-trait EDGE**: Joint analysis of related phenotypes
* **Gene-environment interaction**: EDGE with interaction terms
* **Optimization method benchmarking**: Systematic comparison across data types
* **Transformation selection**: Automated selection of optimal transformation
* **Non-linear relationships**: Spline-based EDGE encoding

Advanced Topics
---------------

Cross-Validation for Alpha Estimation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For extra stability, use k-fold cross-validation:

.. code-block:: python

   from sklearn.model_selection import KFold
   
   kf = KFold(n_splits=5, shuffle=True, random_state=42)
   alpha_estimates = []
   
   for train_idx, val_idx in kf.split(genotype.index):
       edge = EDGEAnalysis(outcome_type='continuous')
       alpha_df = edge.calculate_alpha(
           genotype.iloc[train_idx],
           phenotype.iloc[train_idx],
           outcome, covariates
       )
       alpha_estimates.append(alpha_df)
   
   # Average alpha across folds
   alpha_mean = pd.concat(alpha_estimates).groupby('variant_id')['alpha_value'].mean()

Meta-Analysis of EDGE Results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Combining α estimates across studies:

.. code-block:: python

   def meta_analyze_alpha(alpha_list, n_list):
       """
       Fixed-effects meta-analysis of alpha values.
       
       Args:
           alpha_list: List of alpha estimates from different studies
           n_list: List of sample sizes
       
       Returns:
           Meta-analyzed alpha value
       """
       weights = np.array(n_list) / np.sum(n_list)
       alpha_meta = np.sum(np.array(alpha_list) * weights)
       return alpha_meta
   
   # Example
   study1_alpha = 0.25
   study2_alpha = 0.30
   study3_alpha = 0.28
   
   study_alphas = [study1_alpha, study2_alpha, study3_alpha]
   study_ns = [10000, 15000, 12000]
   
   meta_alpha = meta_analyze_alpha(study_alphas, study_ns)
   print(f"Meta-analyzed alpha: {meta_alpha:.3f}")

Conditional Analysis
~~~~~~~~~~~~~~~~~~~~

Testing for association conditional on a known locus:

.. code-block:: python

   # Add lead SNP as covariate
   lead_snp = 'rs123456'
   phenotype_df['lead_snp_genotype'] = genotype[lead_snp]
   
   # Run EDGE with lead SNP as covariate
   edge = EDGEAnalysis(outcome_type='continuous')
   alpha_df, gwas_results = edge.run_full_analysis(
       train_g, train_p, test_g, test_p,
       outcome='trait',
       covariates=['age', 'sex', 'pc1', 'pc2', 'lead_snp_genotype']
   )

Gene-Based Testing (Future Feature Preview)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Aggregate EDGE effects across genes:

.. code-block:: python

   # Planned for v0.2.0
   from edge_gwas.gene_based import EdgeGeneTest
   
   gene_test = EdgeGeneTest(
       gene_annotation='gencode_v38.gtf',
       window_kb=50  # Include SNPs within 50kb of gene
   )
   
   gene_results = gene_test.run_analysis(
       genotype_data=test_g,
       phenotype_data=test_p,
       alpha_values=alpha_df,
       outcome='trait',
       covariates=['age', 'sex']
   )

Interaction Effects
~~~~~~~~~~~~~~~~~~~

Testing gene-environment interactions with EDGE encoding:

.. code-block:: python

   # Create interaction term
   phenotype_df['bmi_centered'] = phenotype_df['bmi'] - phenotype_df['bmi'].mean()
   
   # For each SNP, create interaction with environment
   def edge_interaction_test(edge_genotype, environment, outcome, covariates):
       """
       Test EDGE genotype × environment interaction.
       
       Args:
           edge_genotype: EDGE-encoded genotype
           environment: Environmental variable (centered)
           outcome: Outcome variable
           covariates: List of covariates
       
       Returns:
           Interaction p-value
       """
       import statsmodels.api as sm
       
       # Create interaction term
       interaction = edge_genotype * environment
       
       # Fit model with main effects and interaction
       X = sm.add_constant(pd.DataFrame({
           'edge_geno': edge_genotype,
           'environment': environment,
           'interaction': interaction,
           **{cov: phenotype_df[cov] for cov in covariates}
       }))
       
       model = sm.OLS(outcome, X)
       result = model.fit()
       
       return result.pvalues['interaction']
   
   # Example usage
   snp = 'rs123456'
   edge_encoded = test_g[snp].replace({0: 1.0, 1: alpha_dict[snp], 2: 0.0})
   
   interaction_p = edge_interaction_test(
       edge_encoded,
       phenotype_df['bmi_centered'],
       phenotype_df['outcome'],
       ['age', 'sex']
   )

Stratified Analysis
~~~~~~~~~~~~~~~~~~~

Running EDGE separately by subgroups:

.. code-block:: python

   # Stratify by sex
   for sex in ['male', 'female']:
       mask = phenotype_df['sex'] == sex
       
       train_g_strat = train_g[mask]
       train_p_strat = train_p[mask]
       test_g_strat = test_g[mask]
       test_p_strat = test_p[mask]
       
       edge = EDGEAnalysis(outcome_type='continuous')
       alpha_df_strat, gwas_strat = edge.run_full_analysis(
           train_g_strat, train_p_strat,
           test_g_strat, test_p_strat,
           outcome='trait',
           covariates=['age', 'pc1', 'pc2'],
           output_prefix=f'results/edge_{sex}'
       )
   
   # Compare alpha values between strata
   alpha_male = pd.read_csv('results/edge_male_alpha_values.csv')
   alpha_female = pd.read_csv('results/edge_female_alpha_values.csv')
   
   comparison = pd.merge(
       alpha_male[['variant_id', 'alpha_value']],
       alpha_female[['variant_id', 'alpha_value']],
       on='variant_id',
       suffixes=('_male', '_female')
   )
   
   comparison['alpha_diff'] = comparison['alpha_value_male'] - comparison['alpha_value_female']
   print(comparison.nlargest(10, 'alpha_diff'))

Polygenic Risk Scores with EDGE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Constructing PRS using EDGE-weighted genotypes:

.. code-block:: python

   def calculate_edge_prs(genotype, alpha_values, effect_sizes):
       """
       Calculate polygenic risk score using EDGE encoding.
       
       Args:
           genotype: Genotype matrix (samples × variants)
           alpha_values: Dict of {variant_id: alpha_value}
           effect_sizes: Dict of {variant_id: effect_size}
       
       Returns:
           PRS for each sample
       """
       prs = np.zeros(len(genotype))
       
       for variant_id in genotype.columns:
           if variant_id in alpha_values and variant_id in effect_sizes:
               # EDGE encode
               geno = genotype[variant_id]
               alpha = alpha_values[variant_id]
               edge_encoded = geno.replace({0: 1.0, 1: alpha, 2: 0.0})
               
               # Weight by effect size
               beta = effect_sizes[variant_id]
               prs += edge_encoded.values * beta
       
       return pd.Series(prs, index=genotype.index)
   
   # Example usage
   alpha_dict = dict(zip(alpha_df['variant_id'], alpha_df['alpha_value']))
   beta_dict = dict(zip(gwas_results['variant_id'], gwas_results['coef']))
   
   prs = calculate_edge_prs(test_g, alpha_dict, beta_dict)
   
   # Evaluate PRS performance
   from sklearn.metrics import roc_auc_score
   auc = roc_auc_score(test_p['case_status'], prs)
   print(f"PRS AUC: {auc:.3f}")

Sensitivity Analyses
~~~~~~~~~~~~~~~~~~~~

**Testing robustness to different splitting ratios:**

.. code-block:: python

   split_ratios = [0.3, 0.4, 0.5, 0.6, 0.7]
   results_by_split = {}
   
   for test_size in split_ratios:
       train_ids, test_ids = train_test_split(
           genotype.index, 
           test_size=test_size, 
           random_state=42
       )
       
       edge = EDGEAnalysis(outcome_type='continuous')
       alpha_df, gwas_df = edge.run_full_analysis(
           genotype.loc[train_ids], phenotype.loc[train_ids],
           genotype.loc[test_ids], phenotype.loc[test_ids],
           outcome='trait',
           covariates=['age', 'sex']
       )
       
       results_by_split[test_size] = {
           'n_sig': (gwas_df['pval'] < 5e-8).sum(),
           'lambda': calculate_lambda(gwas_df['pval']),
           'alpha_mean': alpha_df['alpha_value'].mean()
       }
   
   # Compare results
   sensitivity_df = pd.DataFrame(results_by_split).T
   print(sensitivity_df)

**Testing robustness to different transformations:**

.. code-block:: python

   transformations = [None, 'log', 'inverse_normal', 'rank_inverse_normal']
   results_by_transform = {}
   
   for transform in transformations:
       edge = EDGEAnalysis(
           outcome_type='continuous',
           outcome_transform=transform
       )
       
       alpha_df, gwas_df = edge.run_full_analysis(
           train_g, train_p, test_g, test_p,
           outcome='trait',
           covariates=['age', 'sex']
       )
       
       results_by_transform[str(transform)] = {
           'n_sig': (gwas_df['pval'] < 5e-8).sum(),
           'lambda': calculate_lambda(gwas_df['pval']),
           'min_pval': gwas_df['pval'].min()
       }
   
   transform_df = pd.DataFrame(results_by_transform).T
   print(transform_df)

Integration with Other Tools
-----------------------------

LD Score Regression
~~~~~~~~~~~~~~~~~~~~

Estimate heritability of EDGE results:

.. code-block:: bash

   # Convert EDGE results to LDSC format
   python convert_edge_to_ldsc.py \
       --edge-results edge_gwas_results.csv \
       --output edge_ldsc.sumstats.gz
   
   # Run LDSC
   ldsc.py \
       --h2 edge_ldsc.sumstats.gz \
       --ref-ld-chr eur_w_ld_chr/ \
       --w-ld-chr eur_w_ld_chr/ \
       --out edge_h2

Fine-Mapping
~~~~~~~~~~~~

Use EDGE p-values for fine-mapping:

.. code-block:: python

   # Prepare data for fine-mapping tools (e.g., FINEMAP, SuSiE)
   def prepare_finemap_input(gwas_results, ld_matrix, region):
       """
       Prepare EDGE results for fine-mapping.
       
       Args:
           gwas_results: EDGE GWAS results
           ld_matrix: LD matrix for region
           region: Genomic region (chr:start-end)
       
       Returns:
           Input files for FINEMAP
       """
       # Extract region
       region_results = gwas_results[
           gwas_results['region'] == region
       ].copy()
       
       # Calculate z-scores
       region_results['z_score'] = region_results['coef'] / region_results['std_err']
       
       # Write files
       region_results[['variant_id', 'z_score']].to_csv(
           f'{region}.z', 
           sep=' ', 
           index=False,
           header=False
       )
       
       np.savetxt(f'{region}.ld', ld_matrix)

Colocalization Analysis
~~~~~~~~~~~~~~~~~~~~~~~

Test for shared causal variants between EDGE GWAS and eQTL:

.. code-block:: python

   # Using coloc package (R)
   # Prepare EDGE summary statistics
   
   edge_coloc_input = {
       'beta': gwas_results['coef'],
       'varbeta': gwas_results['std_err']**2,
       'N': gwas_results['n_samples'],
       'type': 'quant',
       'snp': gwas_results['variant_id']
   }
   
   # Run in R:
   # library(coloc)
   # coloc_result <- coloc.abf(
   #     dataset1=edge_coloc_input,
   #     dataset2=eqtl_input
   # )

Pathway Enrichment
~~~~~~~~~~~~~~~~~~

Test for enrichment of EDGE signals in biological pathways:

.. code-block:: python

   from scipy.stats import hypergeom
   
   def pathway_enrichment(sig_snps, pathway_snps, total_snps):
       """
       Test for enrichment using hypergeometric test.
       
       Args:
           sig_snps: Set of significant SNPs
           pathway_snps: Set of SNPs in pathway
           total_snps: Total number of tested SNPs
       
       Returns:
           Enrichment p-value
       """
       overlap = len(sig_snps & pathway_snps)
       n_sig = len(sig_snps)
       n_pathway = len(pathway_snps)
       
       p_value = hypergeom.sf(
           overlap - 1,
           total_snps,
           n_pathway,
           n_sig
       )
       
       return p_value
   
   # Example
   sig_variants = set(gwas_results[gwas_results['pval'] < 5e-8]['variant_id'])
   lipid_pathway_variants = set(['rs123', 'rs456', 'rs789'])  # Example
   
   p_enrich = pathway_enrichment(
       sig_variants,
       lipid_pathway_variants,
       total_snps=len(gwas_results)
   )

Visualization Examples
----------------------

Alpha Distribution Plot
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import matplotlib.pyplot as plt
   import seaborn as sns
   
   plt.figure(figsize=(10, 6))
   sns.histplot(alpha_df['alpha_value'], bins=50, kde=True)
   plt.axvline(0.5, color='red', linestyle='--', label='Additive (α=0.5)')
   plt.axvline(0, color='blue', linestyle='--', label='Recessive (α=0)')
   plt.axvline(1, color='green', linestyle='--', label='Dominant (α=1)')
   plt.xlabel('Alpha Value')
   plt.ylabel('Frequency')
   plt.title('Distribution of EDGE Alpha Values')
   plt.legend()
   plt.savefig('alpha_distribution.png', dpi=300)

Manhattan Plot with Alpha Coloring
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import matplotlib.pyplot as plt
   import numpy as np
   
   def manhattan_with_alpha(gwas_results, alpha_df):
       """
       Manhattan plot colored by alpha values.
       """
       # Merge results with alpha
       plot_data = pd.merge(
           gwas_results,
           alpha_df[['variant_id', 'alpha_value']],
           on='variant_id'
       )
       
       # Classify inheritance
       plot_data['inheritance'] = pd.cut(
           plot_data['alpha_value'],
           bins=[-np.inf, 0.3, 0.7, 1.3, np.inf],
           labels=['Recessive', 'Additive', 'Dominant', 'Over-dominant']
       )
       
       # Plot
       fig, ax = plt.subplots(figsize=(16, 6))
       
       colors = {
           'Recessive': 'blue',
           'Additive': 'gray',
           'Dominant': 'green',
           'Over-dominant': 'red'
       }
       
       for inheritance, color in colors.items():
           subset = plot_data[plot_data['inheritance'] == inheritance]
           ax.scatter(
               range(len(subset)),
               -np.log10(subset['pval']),
               c=color,
               label=inheritance,
               alpha=0.6,
               s=10
           )
       
       ax.axhline(-np.log10(5e-8), color='red', linestyle='--', label='GW significance')
       ax.set_xlabel('Variant Index')
       ax.set_ylabel('-log10(p-value)')
       ax.set_title('EDGE GWAS Manhattan Plot (colored by inheritance pattern)')
       ax.legend()
       plt.tight_layout()
       plt.savefig('manhattan_alpha_colored.png', dpi=300)
   
   manhattan_with_alpha(gwas_results, alpha_df)

QQ Plot Comparison
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   def qq_plot_comparison(edge_pvals, additive_pvals):
       """
       Compare QQ plots for EDGE vs additive GWAS.
       """
       import matplotlib.pyplot as plt
       import numpy as np
       from scipy import stats
       
       fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
       
       def make_qq(pvals, ax, title):
           # Remove NaN
           pvals_clean = pvals[~np.isnan(pvals)]
           
           # Expected and observed
           n = len(pvals_clean)
           expected = -np.log10(np.arange(1, n+1) / (n+1))
           observed = -np.log10(np.sort(pvals_clean))
           
           # Plot
           ax.scatter(expected, observed, alpha=0.5, s=10)
           ax.plot([0, max(expected)], [0, max(expected)], 
                   'r--', label='y=x')
           
           # Calculate lambda
           chi2_stats = stats.chi2.ppf(1 - pvals_clean, df=1)
           lambda_gc = np.median(chi2_stats) / 0.455
           
           ax.set_xlabel('Expected -log10(p)')
           ax.set_ylabel('Observed -log10(p)')
           ax.set_title(f'{title}\nλ = {lambda_gc:.3f}')
           ax.legend()
           ax.grid(True, alpha=0.3)
       
       make_qq(edge_pvals, ax1, 'EDGE GWAS')
       make_qq(additive_pvals, ax2, 'Additive GWAS')
       
       plt.tight_layout()
       plt.savefig('qq_comparison.png', dpi=300)
   
   # Usage
   qq_plot_comparison(
       gwas_results['pval'].values,
       additive_results['pval'].values
   )

Regional Association Plot
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   def regional_plot_with_alpha(gwas_results, region_chr, region_start, region_end):
       """
       Regional association plot showing alpha values.
       """
       import matplotlib.pyplot as plt
       import numpy as np
       
       # Extract region
       region_data = gwas_results[
           (gwas_results['chr'] == region_chr) &
           (gwas_results['pos'] >= region_start) &
           (gwas_results['pos'] <= region_end)
       ].copy()
       
       # Create figure
       fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), 
                                       sharex=True, height_ratios=[2, 1])
       
       # Top panel: -log10(p) colored by alpha
       scatter = ax1.scatter(
           region_data['pos'] / 1e6,
           -np.log10(region_data['pval']),
           c=region_data['alpha_value'],
           cmap='RdYlBu_r',
           vmin=0, vmax=1,
           s=30,
           alpha=0.7
       )
       
       ax1.axhline(-np.log10(5e-8), color='red', linestyle='--', 
                   label='GW significance')
       ax1.set_ylabel('-log10(p-value)')
       ax1.set_title(f'Region: chr{region_chr}:{region_start}-{region_end}')
       ax1.legend()
       ax1.grid(True, alpha=0.3)
       
       # Add colorbar
       cbar = plt.colorbar(scatter, ax=ax1)
       cbar.set_label('Alpha Value')
       
       # Bottom panel: Alpha values
       ax2.scatter(
           region_data['pos'] / 1e6,
           region_data['alpha_value'],
           c='black',
           s=20,
           alpha=0.5
       )
       ax2.axhline(0.5, color='red', linestyle='--', label='Additive')
       ax2.axhline(0, color='blue', linestyle='--', label='Recessive')
       ax2.axhline(1, color='green', linestyle='--', label='Dominant')
       ax2.set_xlabel('Position (Mb)')
       ax2.set_ylabel('Alpha Value')
       ax2.legend(loc='upper right', fontsize=8)
       ax2.grid(True, alpha=0.3)
       ax2.set_ylim(-0.2, 1.5)
       
       plt.tight_layout()
       plt.savefig(f'regional_plot_chr{region_chr}_{region_start}_{region_end}.png', 
                   dpi=300)
   
   # Usage
   regional_plot_with_alpha(gwas_results, 1, 100000000, 102000000)

Effect Size Comparison Plot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   def compare_effect_sizes(edge_results, additive_results):
       """
       Compare effect sizes between EDGE and additive GWAS.
       """
       import matplotlib.pyplot as plt
       import numpy as np
       
       # Merge results
       comparison = pd.merge(
           edge_results[['variant_id', 'coef', 'alpha_value']],
           additive_results[['variant_id', 'coef']],
           on='variant_id',
           suffixes=('_edge', '_additive')
       )
       
       fig, axes = plt.subplots(1, 2, figsize=(14, 6))
       
       # Left panel: Effect size comparison
       axes[0].scatter(
           comparison['coef_additive'],
           comparison['coef_edge'],
           alpha=0.3,
           s=10
       )
       
       # Add diagonal line
       lims = [
           np.min([axes[0].get_xlim(), axes[0].get_ylim()]),
           np.max([axes[0].get_xlim(), axes[0].get_ylim()]),
       ]
       axes[0].plot(lims, lims, 'r--', alpha=0.75, zorder=0)
       axes[0].set_xlabel('Additive Effect Size')
       axes[0].set_ylabel('EDGE Effect Size')
       axes[0].set_title('Effect Size Comparison')
       axes[0].grid(True, alpha=0.3)
       
       # Right panel: Effect difference vs alpha
       comparison['effect_diff'] = comparison['coef_edge'] - comparison['coef_additive']
       
       scatter = axes[1].scatter(
           comparison['alpha_value'],
           comparison['effect_diff'],
           c=np.abs(comparison['effect_diff']),
           cmap='viridis',
           alpha=0.5,
           s=10
       )
       
       axes[1].axhline(0, color='red', linestyle='--')
       axes[1].axvline(0.5, color='blue', linestyle='--', alpha=0.5, 
                       label='Additive (α=0.5)')
       axes[1].set_xlabel('Alpha Value')
       axes[1].set_ylabel('Effect Size Difference (EDGE - Additive)')
       axes[1].set_title('Effect Difference vs Alpha')
       axes[1].legend()
       axes[1].grid(True, alpha=0.3)
       
       plt.colorbar(scatter, ax=axes[1], label='|Effect Difference|')
       plt.tight_layout()
       plt.savefig('effect_comparison.png', dpi=300)
   
   # Usage
   compare_effect_sizes(gwas_results, additive_results)

Convergence Diagnostic Plot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   def plot_convergence_diagnostics(alpha_df, genotype_maf):
       """
       Visualize relationship between convergence and variant properties.
       """
       import matplotlib.pyplot as plt
       import numpy as np
       
       # Merge with MAF
       plot_data = pd.merge(
           alpha_df,
           genotype_maf,
           left_on='variant_id',
           right_index=True
       )
       
       fig, axes = plt.subplots(2, 2, figsize=(12, 10))
       
       # Panel 1: Alpha vs MAF
       axes[0, 0].scatter(plot_data['maf'], plot_data['alpha_value'], 
                          alpha=0.3, s=5)
       axes[0, 0].axhline(0.5, color='red', linestyle='--', label='Additive')
       axes[0, 0].set_xlabel('Minor Allele Frequency')
       axes[0, 0].set_ylabel('Alpha Value')
       axes[0, 0].set_title('Alpha vs MAF')
       axes[0, 0].legend()
       axes[0, 0].grid(True, alpha=0.3)
       
       # Panel 2: Standard error vs MAF
       axes[0, 1].scatter(plot_data['maf'], plot_data['std_err_het'], 
                          alpha=0.3, s=5, label='Het')
       axes[0, 1].scatter(plot_data['maf'], plot_data['std_err_hom'], 
                          alpha=0.3, s=5, label='Hom')
       axes[0, 1].set_xlabel('Minor Allele Frequency')
       axes[0, 1].set_ylabel('Standard Error')
       axes[0, 1].set_title('Standard Error vs MAF')
       axes[0, 1].legend()
       axes[0, 1].grid(True, alpha=0.3)
       axes[0, 1].set_yscale('log')
       
       # Panel 3: Distribution of alpha by MAF bins
       maf_bins = pd.cut(plot_data['maf'], bins=[0, 0.01, 0.05, 0.1, 0.5])
       plot_data['maf_bin'] = maf_bins
       
       axes[1, 0].violinplot(
           [plot_data[plot_data['maf_bin'] == bin]['alpha_value'].dropna() 
            for bin in maf_bins.cat.categories],
           positions=range(len(maf_bins.cat.categories)),
           showmeans=True
       )
       axes[1, 0].axhline(0.5, color='red', linestyle='--')
       axes[1, 0].set_xticks(range(len(maf_bins.cat.categories)))
       axes[1, 0].set_xticklabels(maf_bins.cat.categories, rotation=45)
       axes[1, 0].set_xlabel('MAF Bin')
       axes[1, 0].set_ylabel('Alpha Value')
       axes[1, 0].set_title('Alpha Distribution by MAF')
       axes[1, 0].grid(True, alpha=0.3)
       
       # Panel 4: Sample size vs convergence
       axes[1, 1].hist(plot_data['n_samples'], bins=50, alpha=0.7)
       axes[1, 1].set_xlabel('Sample Size')
       axes[1, 1].set_ylabel('Frequency')
       axes[1, 1].set_title('Sample Size Distribution')
       axes[1, 1].grid(True, alpha=0.3)
       
       plt.tight_layout()
       plt.savefig('convergence_diagnostics.png', dpi=300)

Power Analysis Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   def plot_power_curves():
       """
       Visualize EDGE power for different alpha values.
       """
       import matplotlib.pyplot as plt
       import numpy as np
       from scipy.stats import norm
       
       # Simulation parameters
       n = 10000
       maf = 0.3
       alpha_true_values = [0, 0.25, 0.5, 0.75, 1.0]
       effect_sizes = np.linspace(0, 0.5, 50)
       
       fig, ax = plt.subplots(figsize=(10, 6))
       
       for alpha_true in alpha_true_values:
           power_edge = []
           power_additive = []
           
           for beta in effect_sizes:
               # Simplified power calculation
               # EDGE: uses correct alpha
               ncp_edge = beta * np.sqrt(n * 2 * maf * (1 - maf))
               power_e = 1 - norm.cdf(norm.ppf(1 - 5e-8/2) - ncp_edge)
               power_edge.append(power_e)
               
               # Additive: assumes alpha=0.5
               scaling = 1 - abs(alpha_true - 0.5)  # Power loss when alpha != 0.5
               ncp_add = beta * scaling * np.sqrt(n * 2 * maf * (1 - maf))
               power_a = 1 - norm.cdf(norm.ppf(1 - 5e-8/2) - ncp_add)
               power_additive.append(power_a)
           
           ax.plot(effect_sizes, power_edge, label=f'EDGE (α={alpha_true})', 
                   linewidth=2)
           ax.plot(effect_sizes, power_additive, '--', alpha=0.5,
                   label=f'Additive (α={alpha_true})')
       
       ax.axhline(0.8, color='gray', linestyle=':', label='80% power')
       ax.set_xlabel('Effect Size')
       ax.set_ylabel('Statistical Power')
       ax.set_title('EDGE vs Additive GWAS Power\n(N=10,000, MAF=0.3)')
       ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
       ax.grid(True, alpha=0.3)
       ax.set_ylim([0, 1])
       
       plt.tight_layout()
       plt.savefig('power_curves.png', dpi=300, bbox_inches='tight')

Simulation Code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Simulate GWAS data with known alpha:**

.. code-block:: python

   import numpy as np
   import pandas as pd
   from scipy.stats import bernoulli, norm
   
   def simulate_gwas_data(n_samples, n_variants, maf, alpha_true, beta_true, 
                          h2=0.3, outcome_type='continuous'):
       """
       Simulate GWAS data with specified inheritance pattern.
       
       Args:
           n_samples: Number of samples
           n_variants: Number of variants
           maf: Minor allele frequency
           alpha_true: True alpha value
           beta_true: True effect size
           h2: Heritability (for continuous outcomes)
           outcome_type: 'continuous' or 'binary'
       
       Returns:
           genotype_df, phenotype_df
       """
       # Simulate genotypes under HWE
       genotypes = np.zeros((n_samples, n_variants), dtype=int)
       
       for j in range(n_variants):
           # Allele frequencies
           p = maf  # Frequency of alt allele
           
           # Genotype probabilities under HWE
           p_00 = (1 - p) ** 2
           p_01 = 2 * p * (1 - p)
           p_11 = p ** 2
           
           # Sample genotypes
           geno_probs = np.random.choice([0, 1, 2], size=n_samples, 
                                         p=[p_00, p_01, p_11])
           genotypes[:, j] = geno_probs
       
       # Simulate phenotypes
       # Assume first variant is causal
       causal_geno = genotypes[:, 0]
       
       # Create EDGE encoding for causal variant
       edge_encoded = np.where(causal_geno == 0, 1.0,
                              np.where(causal_geno == 1, alpha_true, 0.0))
       
       # Genetic component
       genetic_component = beta_true * edge_encoded
       
       if outcome_type == 'continuous':
           # Add environmental noise
           var_genetic = np.var(genetic_component)
           var_environmental = var_genetic * (1 - h2) / h2
           
           environmental = np.random.normal(0, np.sqrt(var_environmental), n_samples)
           outcome = genetic_component + environmental
           
       else:  # binary
           # Logistic model
           prob_case = 1 / (1 + np.exp(-genetic_component))
           outcome = bernoulli.rvs(prob_case)
       
       # Create DataFrames
       genotype_df = pd.DataFrame(
           genotypes,
           columns=[f'rs{i}' for i in range(n_variants)]
       )
       
       # Simulate covariates
       age = np.random.normal(50, 10, n_samples)
       sex = np.random.binomial(1, 0.5, n_samples)
       
       phenotype_df = pd.DataFrame({
           'outcome': outcome,
           'age': age,
           'sex': sex
       })
       
       return genotype_df, phenotype_df
   
   # Example usage
   geno, pheno = simulate_gwas_data(
       n_samples=10000,
       n_variants=1000,
       maf=0.3,
       alpha_true=0.2,  # Recessive
       beta_true=0.5,
       outcome_type='continuous'
   )

**Validate EDGE performance:**

.. code-block:: python

   from edge_gwas import EDGEAnalysis
   from sklearn.model_selection import train_test_split
   
   # Split simulated data
   train_ids, test_ids = train_test_split(
       geno.index, test_size=0.5, random_state=42
   )
   
   # Run EDGE
   edge = EDGEAnalysis(outcome_type='continuous')
   alpha_df, gwas_results = edge.run_full_analysis(
       geno.loc[train_ids], pheno.loc[train_ids],
       geno.loc[test_ids], pheno.loc[test_ids],
       outcome='outcome',
       covariates=['age', 'sex']
   )
   
   # Check alpha recovery
   estimated_alpha = alpha_df[alpha_df['variant_id'] == 'rs0']['alpha_value'].values[0]
   print(f"True alpha: {0.2:.3f}")
   print(f"Estimated alpha: {estimated_alpha:.3f}")
   print(f"Error: {abs(estimated_alpha - 0.2):.3f}")
   
   # Check p-value
   causal_pval = gwas_results[gwas_results['variant_id'] == 'rs0']['pval'].values[0]
   print(f"P-value for causal variant: {causal_pval:.2e}")

Sample Size Calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Power calculation for EDGE:**

.. code-block:: python

   from scipy.stats import norm, nct
   import numpy as np
   
   def edge_power_calculation(n, maf, alpha_true, beta, sig_level=5e-8):
       """
       Calculate statistical power for EDGE analysis.
       
       Args:
           n: Sample size
           maf: Minor allele frequency
           alpha_true: True alpha value
           beta: Effect size
           sig_level: Significance threshold
       
       Returns:
           Statistical power (0-1)
       """
       # Genotype frequencies under HWE
       p_00 = (1 - maf) ** 2
       p_01 = 2 * maf * (1 - maf)
       p_11 = maf ** 2
       
       # EDGE encoding variance
       edge_00 = 1.0
       edge_01 = alpha_true
       edge_11 = 0.0
       
       mean_edge = p_00 * edge_00 + p_01 * edge_01 + p_11 * edge_11
       var_edge = (p_00 * edge_00**2 + p_01 * edge_01**2 + 
                   p_11 * edge_11**2 - mean_edge**2)
       
       # Non-centrality parameter
       ncp = beta * np.sqrt(n * var_edge)
       
       # Critical value
       z_crit = norm.ppf(1 - sig_level / 2)
       
       # Power
       power = 1 - norm.cdf(z_crit - ncp) + norm.cdf(-z_crit - ncp)
       
       return power
   
   # Example: Power curves
   import matplotlib.pyplot as plt
   
   sample_sizes = np.arange(1000, 50000, 1000)
   alpha_values = [0, 0.25, 0.5, 0.75, 1.0]
   
   fig, ax = plt.subplots(figsize=(10, 6))
   
   for alpha in alpha_values:
       powers = [edge_power_calculation(n, maf=0.3, alpha_true=alpha, beta=0.1)
                 for n in sample_sizes]
       ax.plot(sample_sizes, powers, label=f'α={alpha}')
   
   ax.axhline(0.8, color='gray', linestyle='--', label='80% power')
   ax.set_xlabel('Sample Size')
   ax.set_ylabel('Statistical Power')
   ax.set_title('EDGE Power by Inheritance Pattern\n(MAF=0.3, β=0.1)')
   ax.legend()
   ax.grid(True, alpha=0.3)
   plt.savefig('power_analysis.png', dpi=300)

**Sample size requirement:**

For 80% power at genome-wide significance (p < 5×10⁻⁸):

.. code-block:: python

   def required_sample_size(maf, alpha_true, beta, power=0.8, sig_level=5e-8):
       """
       Calculate required sample size for desired power.
       
       Returns:
           Required sample size
       """
       from scipy.optimize import fsolve
       
       def power_equation(n):
           return edge_power_calculation(n, maf, alpha_true, beta, sig_level) - power
       
       n_required = fsolve(power_equation, x0=10000)[0]
       return int(np.ceil(n_required))
   
   # Examples
   scenarios = [
       {'maf': 0.3, 'alpha': 0.2, 'beta': 0.1, 'label': 'Recessive'},
       {'maf': 0.3, 'alpha': 0.5, 'beta': 0.1, 'label': 'Additive'},
       {'maf': 0.3, 'alpha': 0.8, 'beta': 0.1, 'label': 'Dominant'},
   ]
   
   print("Sample size requirements for 80% power:")
   print("=" * 60)
   for scenario in scenarios:
       n_req = required_sample_size(
           scenario['maf'], 
           scenario['alpha'], 
           scenario['beta']
       )
       print(f"{scenario['label']:12s}: N = {n_req:,}")

Useful Helper Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Calculate Minor Allele Frequency:**

.. code-block:: python

   def calculate_maf(genotype_series):
       """
       Calculate minor allele frequency from genotype data.
       
       Args:
           genotype_series: Pandas Series with genotypes (0, 1, 2)
       
       Returns:
           MAF value
       """
       valid_geno = genotype_series.dropna()
       
       # Count alleles
       n_alt = (2 * (valid_geno == 2).sum() + 
                (valid_geno == 1).sum())
       n_total = 2 * len(valid_geno)
       
       # Frequency
       freq_alt = n_alt / n_total
       
       # MAF is minimum of alt and ref frequency
       maf = min(freq_alt, 1 - freq_alt)
       
       return maf

**Calculate Hardy-Weinberg Equilibrium p-value:**

.. code-block:: python

   def calculate_hwe(genotype_series, case_control=None):
       """
       Test for Hardy-Weinberg Equilibrium.
       
       Args:
           genotype_series: Pandas Series with genotypes
           case_control: Optional Series indicating cases (1) and controls (0)
                        If provided, tests HWE in controls only
       
       Returns:
           HWE p-value
       """
       from scipy.stats import chi2
       
       # Use controls only if case/control provided
       if case_control is not None:
           geno = genotype_series[case_control == 0]
       else:
           geno = genotype_series
       
       # Count genotypes
       n_00 = (geno == 0).sum()
       n_01 = (geno == 1).sum()
       n_11 = (geno == 2).sum()
       n_total = n_00 + n_01 + n_11
       
       # Allele frequency
       p_alt = (2 * n_11 + n_01) / (2 * n_total)
       p_ref = 1 - p_alt
       
       # Expected counts under HWE
       exp_00 = n_total * p_ref ** 2
       exp_01 = n_total * 2 * p_ref * p_alt
       exp_11 = n_total * p_alt ** 2
       
       # Chi-square test (1 df)
       chi2_stat = ((n_00 - exp_00) ** 2 / exp_00 +
                    (n_01 - exp_01) ** 2 / exp_01 +
                    (n_11 - exp_11) ** 2 / exp_11)
       
       p_value = 1 - chi2.cdf(chi2_stat, df=1)
       
       return p_value

**Calculate Genomic Inflation Factor:**

.. code-block:: python

   def calculate_lambda(pvalues):
       """
       Calculate genomic inflation factor (lambda).
       
       Args:
           pvalues: Array of p-values
       
       Returns:
           Lambda value
       """
       from scipy.stats import chi2
       
       # Remove NaN
       pvals_clean = pvalues[~np.isnan(pvalues)]
       
       # Convert to chi-square statistics
       chi2_stats = chi2.ppf(1 - pvals_clean, df=1)
       
       # Lambda is ratio of observed to expected median
       lambda_gc = np.median(chi2_stats) / chi2.ppf(0.5, df=1)
       
       return lambda_gc

**Apply Genomic Control:**

.. code-block:: python

   def apply_genomic_control(pvalues, lambda_gc):
       """
       Adjust p-values using genomic control.
       
       Args:
           pvalues: Original p-values
           lambda_gc: Genomic inflation factor
       
       Returns:
           Adjusted p-values
       """
       from scipy.stats import chi2
       
       # Convert to chi-square
       chi2_stats = chi2.ppf(1 - pvalues, df=1)
       
       # Adjust by lambda
       chi2_adjusted = chi2_stats / lambda_gc
       
       # Convert back to p-values
       pvalues_adjusted = 1 - chi2.cdf(chi2_adjusted, df=1)
       
       return pvalues_adjusted

**Identify Unrelated Samples:**

.. code-block:: python

   def identify_unrelated_samples(kinship_matrix, threshold=0.0884):
       """
       Identify maximal set of unrelated individuals.
       
       Args:
           kinship_matrix: N×N kinship matrix
           threshold: Kinship threshold (0.0884 = 3rd degree)
       
       Returns:
           List of unrelated sample indices
       """
       import networkx as nx
       
       n = kinship_matrix.shape[0]
       
       # Create graph of related pairs
       G = nx.Graph()
       G.add_nodes_from(range(n))
       
       for i in range(n):
           for j in range(i + 1, n):
               if kinship_matrix[i, j] > threshold:
                   G.add_edge(i, j)
       
       # Find maximum independent set (NP-hard, use greedy approximation)
       unrelated = nx.maximal_independent_set(G)
       
       return sorted(unrelated)

**Convert Between Genotype Formats:**

.. code-block:: python

   def convert_genotype_format(genotype, from_format, to_format):
       """
       Convert between different genotype encoding formats.
       
       Args:
           genotype: Genotype data
           from_format: '012', 'nucleotide', 'dosage'
           to_format: '012', 'nucleotide', 'dosage'
       
       Returns:
           Converted genotype
       """
       if from_format == '012' and to_format == 'nucleotide':
           # Requires allele information
           # Example: 0->AA, 1->AG, 2->GG
           raise NotImplementedError("Requires allele information")
       
       elif  from_format == 'dosage' and to_format == '012':
           # Round dosages to nearest integer
           return np.round(genotype).astype(int)
       
       elif from_format == '012' and to_format == 'dosage':
           # Already compatible (0, 1, 2 are valid dosages)
           return genotype.astype(float)
       
       else:
           raise ValueError(f"Conversion from {from_format} to {to_format} not implemented")

**Annotate Variants with Gene Information:**

.. code-block:: python

   def annotate_variants_with_genes(variants_df, gtf_file, window_kb=0):
       """
       Annotate variants with nearest gene information.
       
       Args:
           variants_df: DataFrame with 'chr', 'pos' columns
           gtf_file: Path to GTF annotation file
           window_kb: Window around gene (in kb)
       
       Returns:
           DataFrame with gene annotations
       """
       import pyranges as pr
       
       # Load GTF
       genes = pr.read_gtf(gtf_file)
       genes = genes[genes.Feature == "gene"]
       
       # Convert variants to ranges
       variants_gr = pr.PyRanges(
           chromosomes=variants_df['chr'].values,
           starts=variants_df['pos'].values,
           ends=variants_df['pos'].values + 1
       )
       
       # Find nearest genes
       nearest = variants_gr.nearest(
           genes,
           suffix="_gene",
           how="upstream"
       )
       
       # Add gene information to variants
       variants_annotated = variants_df.copy()
       variants_annotated['gene_name'] = nearest.Name.values
       variants_annotated['gene_distance'] = nearest.Distance.values
       
       return variants_annotated

See Also
--------

**Documentation:**

* :ref:`index` - Index page of the documentation
* :ref:`installation` - Installation instructions and requirements
* :ref:`quickstart` - Getting started guide with simple examples
* :ref:`user_guide` - Comprehensive user guide and tutorials
* :ref:`api_reference` - Complete API documentation
* :ref:`examples` - Example analyses and case studies
* :ref:`visualization` - Plotting and visualization guide
* :ref:`statistical_model` - Statistical methods and mathematical background
* :ref:`troubleshooting` - Troubleshooting guide and common issues
* :ref:`faq` - Frequently asked questions
* :ref:`citation` - How to cite EDGE in publications
* :ref:`changelog` - Version history and release notes
* :ref:`futureupdates` - Planned features and roadmap

---

*Last updated: 2025-12-25 for edge-gwas v0.1.1*

*For questions or issues, visit:* https://github.com/nicenzhou/edge-gwas/issues
