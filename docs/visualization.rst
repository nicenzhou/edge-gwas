.. _visualization:

Visualization Guide
===================

This guide covers all visualization options in edge-gwas.

Built-in Visualization Functions
---------------------------------

Manhattan Plot
~~~~~~~~~~~~~~

Create Manhattan plots to visualize genome-wide association results.

Basic Usage
^^^^^^^^^^^

.. code-block:: python

   from edge_gwas.visualize import manhattan_plot
   
   manhattan_plot(
       gwas_df,
       output='manhattan.png'
   )

Advanced Options
^^^^^^^^^^^^^^^^

.. code-block:: python

   manhattan_plot(
       gwas_df,
       output='manhattan.png',
       title='EDGE GWAS Manhattan Plot',
       sig_threshold=5e-8,           # Genome-wide significance
       suggestive_threshold=1e-5,    # Suggestive threshold
       figsize=(14, 6),              # Figure size
       colors=['#1f77b4', '#ff7f0e'], # Chromosome colors
       point_size=3,                 # Marker size
       annotate_top=5                # Annotate top N SNPs
   )

Multiple Chromosomes
^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   # List of DataFrames for each chromosome
   chr_results = [gwas_chr1, gwas_chr2, gwas_chr3]
   
   manhattan_plot(
       chr_results,
       output='multi_chr_manhattan.png',
       title='Multi-Chromosome GWAS'
   )

QQ Plot
~~~~~~~

Create QQ plots to assess genomic inflation.

Basic Usage
^^^^^^^^^^^

.. code-block:: python

   from edge_gwas.visualize import qq_plot
   
   lambda_gc = qq_plot(
       gwas_df,
       output='qq_plot.png'
   )
   
   print(f"Genomic inflation factor: {lambda_gc:.3f}")

Advanced Options
^^^^^^^^^^^^^^^^

.. code-block:: python

   lambda_gc = qq_plot(
       gwas_df,
       output='qq_plot.png',
       title='QQ Plot - EDGE GWAS',
       figsize=(8, 8),
       conf_interval=True,    # Show confidence interval
       color='blue',
       alpha=0.6
   )

Alpha Distribution Plot
~~~~~~~~~~~~~~~~~~~~~~~

Visualize the distribution of alpha values.

Basic Usage
^^^^^^^^^^^

.. code-block:: python

   from edge_gwas.visualize import plot_alpha_distribution
   
   plot_alpha_distribution(
       alpha_df,
       output='alpha_distribution.png'
   )

Advanced Options
^^^^^^^^^^^^^^^^

.. code-block:: python

   plot_alpha_distribution(
       alpha_df,
       output='alpha_distribution.png',
       bins=50,                       # Number of histogram bins
       figsize=(10, 6),               # Figure size
       title='Alpha Value Distribution',
       show_models=True,              # Show inheritance model regions
       kde=True,                      # Add kernel density estimate
       color='steelblue'
   )

Custom Visualizations
---------------------

Regional Association Plot
~~~~~~~~~~~~~~~~~~~~~~~~~

Zoom into a specific genomic region.

.. code-block:: python

   import matplotlib.pyplot as plt
   import numpy as np
   
   def regional_plot(gwas_df, chromosome, start_pos, end_pos, output='regional.png'):
       """
       Create regional association plot for a specific locus.
       """
       # Filter to region
       region = gwas_df[
           (gwas_df['chr'] == chromosome) &
           (gwas_df['pos'] >= start_pos) &
           (gwas_df['pos'] <= end_pos)
       ].copy()
       
       # Calculate -log10(p)
       region['-log10p'] = -np.log10(region['pval'])
       
       # Create plot
       fig, ax = plt.subplots(figsize=(12, 6))
       
       # Scatter plot
       scatter = ax.scatter(
           region['pos'] / 1e6,  # Convert to Mb
           region['-log10p'],
           c=region['-log10p'],
           cmap='YlOrRd',
           s=50,
           alpha=0.7,
           edgecolors='black',
           linewidths=0.5
       )
       
       # Add significance threshold
       ax.axhline(-np.log10(5e-8), color='red', linestyle='--', 
                  label='p = 5e-8', linewidth=2)
       
       # Labels
       ax.set_xlabel(f'Position on Chromosome {chromosome} (Mb)', fontsize=12)
       ax.set_ylabel('-log₁₀(p-value)', fontsize=12)
       ax.set_title(f'Regional Plot: Chr{chromosome}:{start_pos/1e6:.2f}-{end_pos/1e6:.2f} Mb', 
                    fontsize=14)
       
       # Colorbar
       cbar = plt.colorbar(scatter, ax=ax)
       cbar.set_label('-log₁₀(p-value)', fontsize=10)
       
       # Legend
       ax.legend(fontsize=10)
       ax.grid(True, alpha=0.3)
       
       plt.tight_layout()
       plt.savefig(output, dpi=300, bbox_inches='tight')
       plt.close()
       
       print(f"Regional plot saved: {output}")
   
   # Usage
   regional_plot(
       gwas_df,
       chromosome=1,
       start_pos=100e6,  # 100 Mb
       end_pos=110e6,    # 110 Mb
       output='chr1_region.png'
   )

Alpha vs Effect Size Plot
~~~~~~~~~~~~~~~~~~~~~~~~~

Visualize relationship between alpha values and effect sizes.

.. code-block:: python

   import matplotlib.pyplot as plt
   import numpy as np
   from edge_gwas.utils import merge_alpha_with_gwas
   
   def alpha_effect_plot(gwas_df, alpha_df, output='alpha_vs_effect.png'):
       """
       Plot alpha values vs effect sizes, colored by significance.
       """
       # Merge data
       merged = merge_alpha_with_gwas(gwas_df, alpha_df)
       
       # Calculate -log10(p)
       merged['-log10p'] = -np.log10(merged['pval'])
       
       # Create plot
       fig, ax = plt.subplots(figsize=(10, 8))
       
       # Color by significance
       colors = merged['-log10p'].values
       scatter = ax.scatter(
           merged['alpha_value'],
           merged['coef'],
           c=colors,
           cmap='viridis',
           s=30,
           alpha=0.6,
           edgecolors='none'
       )
       
       # Add reference lines
       ax.axvline(0.5, color='red', linestyle='--', alpha=0.5, label='Additive (α=0.5)')
       ax.axhline(0, color='gray', linestyle='-', alpha=0.3)
       
       # Inheritance model regions
       ax.axvspan(-0.5, 0.3, alpha=0.1, color='blue', label='Recessive')
       ax.axvspan(0.3, 0.7, alpha=0.1, color='green', label='Additive')
       ax.axvspan(0.7, 1.5, alpha=0.1, color='orange', label='Dominant')
       
       # Labels
       ax.set_xlabel('Alpha Value', fontsize=12)
       ax.set_ylabel('Effect Size (β)', fontsize=12)
       ax.set_title('Alpha Values vs Effect Sizes', fontsize=14)
       
       # Colorbar
       cbar = plt.colorbar(scatter, ax=ax)
       cbar.set_label('-log₁₀(p-value)', fontsize=10)
       
       # Legend
       ax.legend(loc='upper right', fontsize=9)
       ax.grid(True, alpha=0.3)
       
       plt.tight_layout()
       plt.savefig(output, dpi=300, bbox_inches='tight')
       plt.close()
       
       print(f"Alpha vs effect plot saved: {output}")
   
   # Usage
   alpha_effect_plot(gwas_df, alpha_df, 'alpha_vs_effect.png')

Forest Plot for Top Hits
~~~~~~~~~~~~~~~~~~~~~~~~

Create forest plot for significant associations.

.. code-block:: python

   import matplotlib.pyplot as plt
   import numpy as np
   
   def forest_plot(gwas_df, top_n=20, output='forest_plot.png'):
       """
       Create forest plot for top N associations.
       """
       # Get top hits
       top_hits = gwas_df.nsmallest(top_n, 'pval').copy()
       top_hits = top_hits.sort_values('coef')
       
       # Calculate 95% CI
       top_hits['ci_lower'] = top_hits['coef'] - 1.96 * top_hits['std_err']
       top_hits['ci_upper'] = top_hits['coef'] + 1.96 * top_hits['std_err']
       
       # Create plot
       fig, ax = plt.subplots(figsize=(10, top_n * 0.4 + 2))
       
       y_positions = np.arange(len(top_hits))
       
       # Plot error bars
       ax.errorbar(
           top_hits['coef'],
           y_positions,
           xerr=[
               top_hits['coef'] - top_hits['ci_lower'],
               top_hits['ci_upper'] - top_hits['coef']
           ],
           fmt='o',
           markersize=8,
           capsize=5,
           capthick=2,
           color='steelblue',
           ecolor='gray',
           alpha=0.7
       )
       
       # Add vertical line at 0
       ax.axvline(0, color='red', linestyle='--', linewidth=2, alpha=0.5)
       
       # Labels
       ax.set_yticks(y_positions)
       ax.set_yticklabels(top_hits['variant_id'], fontsize=9)
       ax.set_xlabel('Effect Size (95% CI)', fontsize=12)
       ax.set_title(f'Top {top_n} Associations - Forest Plot', fontsize=14)
       
       # Add p-values as annotations
       for i, (idx, row) in enumerate(top_hits.iterrows()):
           ax.text(
               ax.get_xlim()[1] * 0.95,
               i,
               f"p={row['pval']:.2e}",
               fontsize=8,
               ha='right',
               va='center'
           )
       
       ax.grid(True, alpha=0.3, axis='x')
       plt.tight_layout()
       plt.savefig(output, dpi=300, bbox_inches='tight')
       plt.close()
       
       print(f"Forest plot saved: {output}")
   
   # Usage
   forest_plot(gwas_df, top_n=20, output='top_hits_forest.png')

Comparing Multiple Analyses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Side-by-side comparison plots.

.. code-block:: python

   import matplotlib.pyplot as plt
   import numpy as np
   
   def compare_analyses(gwas_list, labels, output='comparison.png'):
       """
       Create side-by-side Manhattan and QQ plots for multiple analyses.
       """
       n_analyses = len(gwas_list)
       
       fig, axes = plt.subplots(n_analyses, 2, figsize=(16, 4 * n_analyses))
       
       if n_analyses == 1:
           axes = axes.reshape(1, -1)
       
       for i, (gwas_df, label) in enumerate(zip(gwas_list, labels)):
           # Manhattan plot
           gwas_df['-log10p'] = -np.log10(gwas_df['pval'])
           
           colors = ['#1f77b4', '#ff7f0e']
           x_pos = 0
           
           for chrom in sorted(gwas_df['chr'].unique()):
               data = gwas_df[gwas_df['chr'] == chrom]
               axes[i, 0].scatter(
                   x_pos + np.arange(len(data)),
                   data['-log10p'],
                   c=colors[int(chrom) % 2],
                   s=2,
                   alpha=0.7
               )
               x_pos += len(data)
           
           axes[i, 0].axhline(-np.log10(5e-8), color='red', linestyle='--', linewidth=1)
           axes[i, 0].set_ylabel('-log₁₀(p-value)', fontsize=10)
           axes[i, 0].set_title(f'{label} - Manhattan Plot', fontsize=12)
           axes[i, 0].grid(True, alpha=0.3, axis='y')
           
           # QQ plot
           pvals = gwas_df['pval'].dropna()
           pvals = pvals[pvals > 0]
           n = len(pvals)
           
           observed = -np.log10(np.sort(pvals))
           expected = -np.log10(np.arange(1, n + 1) / (n + 1))
           
           axes[i, 1].scatter(expected, observed, s=10, alpha=0.5, color='blue')
           
           max_val = max(expected.max(), observed.max())
           axes[i, 1].plot([0, max_val], [0, max_val], 'r--', linewidth=2)
           
           # Calculate lambda
           from scipy import stats
           chisq = stats.chi2.ppf(1 - pvals, df=1)
           lambda_gc = np.median(chisq) / stats.chi2.ppf(0.5, df=1)
           
           axes[i, 1].text(
               0.05, 0.95,
               f'λ = {lambda_gc:.3f}',
               transform=axes[i, 1].transAxes,
               fontsize=10,
               verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5)
           )
           
           axes[i, 1].set_xlabel('Expected -log₁₀(p)', fontsize=10)
           axes[i, 1].set_ylabel('Observed -log₁₀(p)', fontsize=10)
           axes[i, 1].set_title(f'{label} - QQ Plot', fontsize=12)
           axes[i, 1].grid(True, alpha=0.3)
       
       plt.tight_layout()
       plt.savefig(output, dpi=300, bbox_inches='tight')
       plt.close()
       
       print(f"Comparison plot saved: {output}")
   
   # Usage
   compare_analyses(
       [edge_gwas, additive_gwas, another_gwas],
       labels=['EDGE', 'Additive', 'Alternative'],
       output='method_comparison.png'
   )

Interactive Plots with Plotly
------------------------------

For interactive exploration of results.

.. code-block:: python

   import plotly.graph_objects as go
   import plotly.express as px
   import numpy as np
   
   def interactive_manhattan(gwas_df, output='manhattan_interactive.html'):
       """
       Create interactive Manhattan plot using Plotly.
       """
       gwas_df['-log10p'] = -np.log10(gwas_df['pval'])
       
       # Assign colors by chromosome
       colors = px.colors.qualitative.Plotly
       gwas_df['color'] = gwas_df['chr'].apply(lambda x: colors[int(x) % len(colors)])
       
       # Create cumulative position
       gwas_df = gwas_df.sort_values(['chr', 'pos'])
       gwas_df['x_pos'] = 0
       
       offset = 0
       for chrom in sorted(gwas_df['chr'].unique()):
           mask = gwas_df['chr'] == chrom
           gwas_df.loc[mask, 'x_pos'] = np.arange(mask.sum()) + offset
           offset += mask.sum()
       
       # Create figure
       fig = go.Figure()
       
       for chrom in sorted(gwas_df['chr'].unique()):
           data = gwas_df[gwas_df['chr'] == chrom]
           
           fig.add_trace(go.Scatter(
               x=data['x_pos'],
               y=data['-log10p'],
               mode='markers',
               name=f'Chr {chrom}',
               marker=dict(
                   size=4,
                   color=data['color'].iloc[0],
                   opacity=0.7
               ),
               text=[f"Variant: {v}<br>Position: {p}<br>P-value: {pv:.2e}" 
                     for v, p, pv in zip(data['variant_id'], data['pos'], data['pval'])],
               hoverinfo='text'
           ))
       
       # Add significance line
       fig.add_hline(
           y=-np.log10(5e-8),
           line_dash="dash",
           line_color="red",
           annotation_text="p = 5e-8"
       )
       
       fig.update_layout(
           title='Interactive EDGE GWAS Manhattan Plot',
           xaxis_title='Chromosome',
           yaxis_title='-log₁₀(p-value)',
           height=600,
           hovermode='closest',
           showlegend=True
       )
       
       fig.write_html(output)
       print(f"Interactive Manhattan plot saved: {output}")
   
   # Usage
   interactive_manhattan(gwas_df, 'manhattan_interactive.html')

Volcano Plot
~~~~~~~~~~~~

Visualize effect size vs significance.

.. code-block:: python

   import matplotlib.pyplot as plt
   import numpy as np
   
   def volcano_plot(gwas_df, output='volcano.png', effect_threshold=0.3, sig_threshold=5e-8):
       """
       Create volcano plot of effect size vs -log10(p-value).
       """
       # Calculate -log10(p)
       gwas_df['-log10p'] = -np.log10(gwas_df['pval'])
       
       # Categorize points
       gwas_df['category'] = 'Not significant'
       gwas_df.loc[
           (gwas_df['pval'] < sig_threshold) & (abs(gwas_df['coef']) >= effect_threshold),
           'category'
       ] = 'Significant'
       gwas_df.loc[
           (gwas_df['pval'] < sig_threshold) & (abs(gwas_df['coef']) < effect_threshold),
           'category'
       ] = 'Significant, small effect'
       
       # Create plot
       fig, ax = plt.subplots(figsize=(10, 8))
       
       # Plot by category
       colors = {
           'Not significant': 'gray',
           'Significant, small effect': 'orange',
           'Significant': 'red'
       }
       
       for category, color in colors.items():
           data = gwas_df[gwas_df['category'] == category]
           ax.scatter(
               data['coef'],
               data['-log10p'],
               c=color,
               label=category,
               s=20,
               alpha=0.6,
               edgecolors='none'
           )
       
       # Add threshold lines
       ax.axhline(-np.log10(sig_threshold), color='blue', linestyle='--', 
                  linewidth=2, alpha=0.7, label=f'p = {sig_threshold}')
       ax.axvline(-effect_threshold, color='green', linestyle='--', 
                  linewidth=1, alpha=0.5)
       ax.axvline(effect_threshold, color='green', linestyle='--', 
                  linewidth=1, alpha=0.5, label=f'Effect = ±{effect_threshold}')
       
       # Labels
       ax.set_xlabel('Effect Size (β)', fontsize=12)
       ax.set_ylabel('-log₁₀(p-value)', fontsize=12)
       ax.set_title('Volcano Plot - EDGE GWAS', fontsize=14)
       ax.legend(loc='upper right', fontsize=9)
       ax.grid(True, alpha=0.3)
       
       plt.tight_layout()
       plt.savefig(output, dpi=300, bbox_inches='tight')
       plt.close()
       
       print(f"Volcano plot saved: {output}")
   
   # Usage
   volcano_plot(gwas_df, 'volcano.png', effect_threshold=0.3, sig_threshold=5e-8)

Multi-Panel Summary Figure
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create comprehensive summary figure with multiple panels.

.. code-block:: python

   import matplotlib.pyplot as plt
   import matplotlib.gridspec as gridspec
   import numpy as np
   from scipy import stats
   
   def summary_figure(gwas_df, alpha_df, output='summary_figure.png'):
       """
       Create multi-panel summary figure with:
       - Manhattan plot
       - QQ plot
       - Alpha distribution
       - Top hits table
       """
       # Prepare data
       gwas_df['-log10p'] = -np.log10(gwas_df['pval'])
       
       # Create figure with custom layout
       fig = plt.figure(figsize=(16, 12))
       gs = gridspec.GridSpec(3, 2, height_ratios=[2, 1, 1], hspace=0.3, wspace=0.3)
       
       # Panel 1: Manhattan plot (top, spanning both columns)
       ax1 = fig.add_subplot(gs[0, :])
       
       colors = ['#1f77b4', '#ff7f0e']
       x_pos = 0
       chr_centers = []
       
       for chrom in sorted(gwas_df['chr'].unique()):
           data = gwas_df[gwas_df['chr'] == chrom]
           ax1.scatter(
               x_pos + np.arange(len(data)),
               data['-log10p'],
               c=colors[int(chrom) % 2],
               s=3,
               alpha=0.7
           )
           chr_centers.append(x_pos + len(data) / 2)
           x_pos += len(data)
       
       ax1.axhline(-np.log10(5e-8), color='red', linestyle='--', linewidth=2, label='p = 5e-8')
       ax1.axhline(-np.log10(1e-5), color='blue', linestyle=':', linewidth=1.5, label='p = 1e-5')
       ax1.set_xticks(chr_centers)
       ax1.set_xticklabels([str(c) for c in sorted(gwas_df['chr'].unique())])
       ax1.set_xlabel('Chromosome', fontsize=12)
       ax1.set_ylabel('-log₁₀(p-value)', fontsize=12)
       ax1.set_title('Manhattan Plot', fontsize=14, fontweight='bold')
       ax1.legend(loc='upper right')
       ax1.grid(True, alpha=0.3, axis='y')
       
       # Panel 2: QQ plot (middle left)
       ax2 = fig.add_subplot(gs[1, 0])
       
       pvals = gwas_df['pval'].dropna()
       pvals = pvals[pvals > 0]
       n = len(pvals)
       
       observed = -np.log10(np.sort(pvals))
       expected = -np.log10(np.arange(1, n + 1) / (n + 1))
       
       ax2.scatter(expected, observed, s=15, alpha=0.5, color='blue')
       max_val = max(expected.max(), observed.max())
       ax2.plot([0, max_val], [0, max_val], 'r--', linewidth=2)
       
       # Calculate lambda
       chisq = stats.chi2.ppf(1 - pvals, df=1)
       lambda_gc = np.median(chisq) / stats.chi2.ppf(0.5, df=1)
       
       ax2.text(
           0.05, 0.95,
           f'λ = {lambda_gc:.3f}',
           transform=ax2.transAxes,
           fontsize=11,
           verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7)
       )
       
       ax2.set_xlabel('Expected -log₁₀(p)', fontsize=11)
       ax2.set_ylabel('Observed -log₁₀(p)', fontsize=11)
       ax2.set_title('QQ Plot', fontsize=12, fontweight='bold')
       ax2.grid(True, alpha=0.3)
       
       # Panel 3: Alpha distribution (middle right)
       ax3 = fig.add_subplot(gs[1, 1])
       
       ax3.hist(alpha_df['alpha_value'], bins=50, color='steelblue', 
                edgecolor='black', alpha=0.7)
       
       # Add vertical lines for inheritance models
       ax3.axvline(0, color='purple', linestyle='--', linewidth=2, alpha=0.5, label='Recessive')
       ax3.axvline(0.5, color='green', linestyle='--', linewidth=2, alpha=0.5, label='Additive')
       ax3.axvline(1, color='orange', linestyle='--', linewidth=2, alpha=0.5, label='Dominant')
       
       ax3.set_xlabel('Alpha Value', fontsize=11)
       ax3.set_ylabel('Frequency', fontsize=11)
       ax3.set_title('Alpha Distribution', fontsize=12, fontweight='bold')
       ax3.legend(loc='upper right', fontsize=9)
       ax3.grid(True, alpha=0.3, axis='y')
       
       # Panel 4: Top hits table (bottom, spanning both columns)
       ax4 = fig.add_subplot(gs[2, :])
       ax4.axis('off')
       
       # Get top 10 hits
       top_hits = gwas_df.nsmallest(10, 'pval')[
           ['variant_id', 'chr', 'pos', 'pval', 'coef', 'std_err']
       ].copy()
       
       # Format p-values
       top_hits['pval'] = top_hits['pval'].apply(lambda x: f'{x:.2e}')
       top_hits['coef'] = top_hits['coef'].apply(lambda x: f'{x:.4f}')
       top_hits['std_err'] = top_hits['std_err'].apply(lambda x: f'{x:.4f}')
       
       # Create table
       table = ax4.table(
           cellText=top_hits.values,
           colLabels=['Variant ID', 'Chr', 'Position', 'P-value', 'Effect', 'Std Err'],
           cellLoc='center',
           loc='center',
           bbox=[0, 0, 1, 1]
       )
       
       table.auto_set_font_size(False)
       table.set_fontsize(9)
       table.scale(1, 2)
       
       # Style header
       for i in range(len(top_hits.columns)):
           table[(0, i)].set_facecolor('#4CAF50')
           table[(0, i)].set_text_props(weight='bold', color='white')
       
       # Alternate row colors
       for i in range(1, len(top_hits) + 1):
           for j in range(len(top_hits.columns)):
               if i % 2 == 0:
                   table[(i, j)].set_facecolor('#f0f0f0')
       
       ax4.set_title('Top 10 Significant Associations', fontsize=12, 
                     fontweight='bold', pad=20)
       
       # Add summary statistics text
       fig.text(
           0.5, 0.02,
           f'Total variants: {len(gwas_df):,} | '
           f'Significant (p<5e-8): {(gwas_df["pval"] < 5e-8).sum()} | '
           f'Suggestive (p<1e-5): {(gwas_df["pval"] < 1e-5).sum()} | '
           f'λ: {lambda_gc:.3f}',
           ha='center',
           fontsize=11,
           bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3)
       )
       
       plt.savefig(output, dpi=300, bbox_inches='tight')
       plt.close()
       
       print(f"Summary figure saved: {output}")
   
   # Usage
   summary_figure(gwas_df, alpha_df, 'edge_gwas_summary.png')

Plotting Tips
-------------

Color Schemes
~~~~~~~~~~~~~

Recommended color schemes for different plot types:

**Manhattan plots:**

.. code-block:: python

   # Classic
   colors = ['#1f77b4', '#ff7f0e']
   
   # Colorblind-friendly
   colors = ['#0173B2', '#DE8F05']
   
   # Publication quality
   colors = ['#2E3440', '#88C0D0']

**QQ plots:**

.. code-block:: python

   # Points
   point_color = '#3B4252'
   
   # Diagonal line
   line_color = '#BF616A'
   
   # Confidence interval
   ci_color = '#D8DEE9'

**Alpha distribution:**

.. code-block:: python

   # Histogram
   hist_color = '#5E81AC'
   
   # KDE line
   kde_color = '#BF616A'

Figure Resolution
~~~~~~~~~~~~~~~~~

For different output types:

.. code-block:: python

   # Screen display
   plt.savefig('plot.png', dpi=100)
   
   # Publication (print)
   plt.savefig('plot.png', dpi=300)
   
   # Poster
   plt.savefig('plot.png', dpi=600)
   
   # Vector format (infinite resolution)
   plt.savefig('plot.pdf')
   plt.savefig('plot.svg')

Font Sizes
~~~~~~~~~~

Recommended font sizes for readability:

.. code-block:: python

   # For standard figures (10x6 inches)
   title_size = 14
   label_size = 12
   tick_size = 10
   legend_size = 10
   
   # For presentations (larger)
   title_size = 18
   label_size = 16
   tick_size = 14
   legend_size = 14
   
   # For publications (smaller, more compact)
   title_size = 12
   label_size = 10
   tick_size = 8
   legend_size = 8

Exporting Plots
~~~~~~~~~~~~~~~

Multiple format export:

.. code-block:: python

   def save_plot_multiple_formats(fig, basename):
       """Save plot in multiple formats"""
       formats = ['png', 'pdf', 'svg', 'eps']
       for fmt in formats:
           filename = f'{basename}.{fmt}'
           fig.savefig(filename, dpi=300, bbox_inches='tight')
           print(f'Saved: {filename}')
   
   # Usage
   fig, ax = plt.subplots()
   # ... create your plot ...
   save_plot_multiple_formats(fig, 'my_plot')

See Also
--------

**Documentation:**

* `Documentation Home <index.html>`_ - Home
* :ref:`installation` - Installation instructions and requirements
* :ref:`quickstart` - Getting started guide with simple examples
* :ref:`statistical_model` - Statistical methods and mathematical background
* :ref:`examples` - Example analyses and case studies
* :ref:`visualization` - Plotting and visualization guide
* :ref:`api_reference` - Complete API documentation
* :ref:`troubleshooting` - Troubleshooting guide and common issues
* :ref:`faq` - Frequently asked questions
* :ref:`citation` - How to cite EDGE in publications
* :ref:`changelog` - Version history and release notes
* :ref:`futureupdates` - Planned features and roadmap

---

*Last updated: 2025-12-28 for edge-gwas v0.1.1*

*For questions or issues, visit:* https://github.com/nicenzhou/edge-gwas/issues
