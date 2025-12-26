"""Example: Creating visualizations from EDGE GWAS results"""

import pandas as pd
from edge_gwas import manhattan_plot, qq_plot, plot_alpha_distribution

# Load results
gwas_df = pd.read_csv('edge_results_gwas.txt', sep='\t')
alpha_df = pd.read_csv('edge_results_alpha.txt', sep='\t')

# Create Manhattan plot
manhattan_plot(
    gwas_df, 
    output='manhattan.png',
    title='My EDGE GWAS Study',
    sig_threshold=5e-8
)

# Create QQ plot and get lambda
lambda_gc = qq_plot(gwas_df, output='qq_plot.png')

# Plot alpha distribution
plot_alpha_distribution(alpha_df, output='alpha_dist.png')

print(f"Genomic inflation factor: {lambda_gc:.3f}")
