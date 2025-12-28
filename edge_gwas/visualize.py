"""
Visualization functions for EDGE GWAS results
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from typing import Union, List


def manhattan_plot(
    gwas_df: Union[pd.DataFrame, List[pd.DataFrame]], 
    output: str = 'manhattan.png',
    title: str = 'EDGE GWAS Manhattan Plot',
    sig_threshold: float = 5e-8,
    figsize: tuple = (14, 6),
    colors: List[str] = None
) -> None:
    """
    Create Manhattan plot from EDGE GWAS results.
    
    Args:
        gwas_df: DataFrame or list of DataFrames with columns 'chrom', 'pos', 'pval'
        output: Output filename for the plot
        title: Plot title
        sig_threshold: Genome-wide significance threshold (default: 5e-8)
        figsize: Figure size as (width, height)
        colors: List of two colors for alternating chromosomes
        
    Returns:
        None (saves plot to file)
    """
    # Handle list of DataFrames
    if isinstance(gwas_df, list):
        gwas_df = pd.concat(gwas_df, ignore_index=True)
    else:
        gwas_df = gwas_df.copy()
    
    # Check for required columns
    if 'chrom' not in gwas_df.columns or 'pval' not in gwas_df.columns:
        raise ValueError("gwas_df must contain 'chrom' and 'pval' columns")
    
    # Remove rows with missing values
    gwas_df = gwas_df.dropna(subset=['chrom', 'pval'])
    
    # Calculate -log10(p)
    gwas_df['-log10p'] = -np.log10(gwas_df['pval'])
    
    # Set default colors
    if colors is None:
        colors = ['#1f77b4', '#ff7f0e']
    
    # Create plot
    fig, ax = plt.subplots(figsize=figsize)
    x_pos = 0
    x_labels = []
    x_labels_pos = []
    
    # Sort chromosomes
    unique_chroms = sorted(gwas_df['chrom'].unique(), key=lambda x: (int(x) if str(x).isdigit() else float('inf'), str(x)))
    
    for chrom in unique_chroms:
        data = gwas_df[gwas_df['chrom'] == chrom]
        
        # Sort by position if available
        if 'pos' in gwas_df.columns:
            data = data.sort_values('pos')
        
        # Determine color index
        if str(chrom).isdigit():
            color_idx = int(chrom) % 2
        else:
            color_idx = 0
        
        ax.scatter(
            x_pos + np.arange(len(data)), 
            data['-log10p'], 
            c=colors[color_idx], 
            s=2, 
            alpha=0.7
        )
        
        x_labels.append(str(chrom))
        x_labels_pos.append(x_pos + len(data) / 2)
        x_pos += len(data)
    
    # Add significance line
    ax.axhline(-np.log10(sig_threshold), color='red', linestyle='--', 
               linewidth=1, label=f'p={sig_threshold}')
    
    # Format plot
    ax.set_xlabel('Chromosome', fontsize=12)
    ax.set_ylabel('-log₁₀(p-value)', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Manhattan plot saved to {output}")


def qq_plot(
    gwas_df: Union[pd.DataFrame, List[pd.DataFrame]], 
    output: str = 'qq_plot.png',
    title: str = 'EDGE GWAS QQ Plot',
    figsize: tuple = (8, 8)
) -> float:
    """
    Create QQ plot from EDGE GWAS results and calculate genomic inflation factor.
    
    Args:
        gwas_df: DataFrame or list of DataFrames with column 'pval'
        output: Output filename for the plot
        title: Plot title
        figsize: Figure size as (width, height)
        
    Returns:
        lambda_gc: Genomic inflation factor
    """
    # Handle list of DataFrames
    if isinstance(gwas_df, list):
        gwas_df = pd.concat(gwas_df, ignore_index=True)
    
    # Get p-values
    pvals = gwas_df['pval'].dropna()
    pvals = pvals[pvals > 0]
    
    if len(pvals) == 0:
        raise ValueError("No valid p-values found")
    
    # Calculate expected vs observed
    n = len(pvals)
    observed = -np.log10(np.sort(pvals))
    expected = -np.log10(np.arange(1, n + 1) / (n + 1))
    
    # Genomic inflation factor
    chisq = stats.chi2.ppf(1 - pvals, df=1)
    lambda_gc = np.median(chisq) / stats.chi2.ppf(0.5, df=1)
    
    # Create plot
    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(expected, observed, s=10, alpha=0.5, color='blue')
    
    # Add diagonal line
    max_val = max(expected.max(), observed.max())
    ax.plot([0, max_val], [0, max_val], 'r--', linewidth=2, label='Expected (null)')
    
    # Add lambda annotation
    ax.text(0.05, 0.95, f'λ = {lambda_gc:.3f}', 
            transform=ax.transAxes, 
            fontsize=12, 
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Format plot
    ax.set_xlabel('Expected -log₁₀(p)', fontsize=12)
    ax.set_ylabel('Observed -log₁₀(p)', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"QQ plot saved to {output}")
    print(f"Genomic inflation factor (λ): {lambda_gc:.3f}")
    
    return lambda_gc


def plot_alpha_distribution(
    alpha_df: pd.DataFrame,
    output: str = 'alpha_distribution.png',
    bins: int = 50,
    figsize: tuple = (10, 6),
    xlim: tuple = None
) -> None:
    """
    Plot distribution of alpha values.
    
    Args:
        alpha_df: DataFrame with 'alpha_value' column
        output: Output filename
        bins: Number of histogram bins
        figsize: Figure size
        xlim: Optional tuple (min, max) for x-axis limits, e.g., (-3, 3)
              If None, uses full range of data
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    alpha_values = alpha_df['alpha_value'].dropna()
    
    # Filter values if xlim is specified
    if xlim is not None:
        alpha_values_plot = alpha_values[(alpha_values >= xlim[0]) & (alpha_values <= xlim[1])]
        n_filtered = len(alpha_values) - len(alpha_values_plot)
        if n_filtered > 0:
            print(f"Note: {n_filtered} alpha values outside range [{xlim[0]}, {xlim[1]}] excluded from plot")
    else:
        alpha_values_plot = alpha_values
    
    ax.hist(alpha_values_plot, bins=bins, edgecolor='black', alpha=0.7)
    ax.axvline(0.5, color='red', linestyle='--', linewidth=2, label='Additive (α=0.5)')
    ax.axvline(0, color='orange', linestyle='--', linewidth=2, label='Recessive (α=0)')
    ax.axvline(1, color='green', linestyle='--', linewidth=2, label='Dominant (α=1)')
    
    # Set x-axis limits if specified
    if xlim is not None:
        ax.set_xlim(xlim)
    
    ax.set_xlabel('Alpha Value', fontsize=12)
    ax.set_ylabel('Frequency', fontsize=12)
    ax.set_title('Distribution of EDGE Alpha Values', fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Alpha distribution plot saved to {output}")
    print(f"Total alpha values: {len(alpha_values)}")
    if xlim is not None:
        print(f"Values plotted (within [{xlim[0]}, {xlim[1]}]): {len(alpha_values_plot)}")
    print(f"Alpha statistics:")
    print(f"  Mean: {alpha_values.mean():.3f}")
    print(f"  Median: {alpha_values.median():.3f}")
    print(f"  Std: {alpha_values.std():.3f}")
    print(f"  Min: {alpha_values.min():.3f}")
    print(f"  Max: {alpha_values.max():.3f}")
