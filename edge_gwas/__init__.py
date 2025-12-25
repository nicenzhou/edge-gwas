"""
EDGE-GWAS: Encoding Deviation Genotypic Effects GWAS

A Python package for performing EDGE GWAS analysis with flexible encoding
to identify nonadditive SNP effects.

References:
    Zhou, J. et al. Flexibly encoded GWAS identifies novel nonadditive SNPs 
    in individuals of African and European ancestry. medRxiv 2023.06.01.23290857
    
    Hall, M. A. et al. Novel EDGE encoding method enhances ability to identify 
    genetic interactions. PLoS Genetics 17, e1009534 (2021).
"""

__version__ = "0.0.0"
__author__ = "Jiayan Zhou"
__email__ = "jyzhou@stanford.edu"

from .core import EDGEAnalysis
from .utils import load_plink_data, prepare_phenotype_data
from .visualize import manhattan_plot, qq_plot, plot_alpha_distribution

__all__ = [
    "EDGEAnalysis",
    "load_plink_data",
    "prepare_phenotype_data",
    "manhattan_plot",
    "qq_plot",
    "plot_alpha_distribution",
]
