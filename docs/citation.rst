.. _citation:

Citation
========

If you use edge-gwas in your research, please cite:

Primary Citation
----------------

Zhou, J., Rico, A. L. G., Guare, L., et al. (2023). 
Flexibly encoded genome-wide association study identifies novel nonadditive 
genetic risk variants for cardiometabolic traits. 
*medRxiv*, 2023.06.01.23290857.

https://doi.org/10.1101/2023.06.01.23290857

BibTeX Entry
------------

.. code-block:: bibtex

   @article{zhou2023edgegwas,
     title={Flexibly encoded genome-wide association study identifies novel nonadditive genetic risk variants for cardiometabolic traits},
     author={Zhou, Jiayan and Rico, Andre Luis Garao and Guare, Lindsay and others},
     journal={medRxiv},
     year={2023},
     doi={10.1101/2023.06.01.23290857},
     url={https://doi.org/10.1101/2023.06.01.23290857}
   }

Software Citation
-----------------

For the edge-gwas Python package specifically:

**APA Style:**

Zhou, J., Hall, M. A., & Contributors. (2025). 
edge-gwas: Python package for Elastic Data-Driven Encoding GWAS (Version 0.1.0) [Computer software]. 
https://github.com/YOUR-USERNAME/edge-gwas

**BibTeX:**

.. code-block:: bibtex

   @software{edge_gwas_2025,
     title={edge-gwas: Python package for Elastic Data-Driven Encoding GWAS},
     author={Zhou, Jiayan and Hall, Molly Ann and {edge-gwas Contributors}},
     year={2025},
     version={0.1.0},
     url={https://github.com/YOUR-USERNAME/edge-gwas},
     license={GPL-3.0}
   }

Related Publications
--------------------

If you use EDGE methodology in your research, you may also want to cite:

**Original EDGE methodology papers** (add as they become available):

* Zhou, J., et al. (2023). Detailed methodology paper [In preparation]
* Hall, M. A., et al. (2023). Applications in cardiometabolic disease [In preparation]

Citing Specific Features
-------------------------

If you use specific features of edge-gwas:

**Two-stage analysis approach:**

Zhou, J., et al. (2023). Two-stage flexibly encoded GWAS methodology. 
*medRxiv*, 2023.06.01.23290857.

**Visualization tools:**

edge-gwas Contributors (2025). Visualization tools in edge-gwas Python package 
(Version 0.1.0). https://github.com/YOUR-USERNAME/edge-gwas

Example Citation in Methods Section
------------------------------------

**For a typical research paper:**

.. code-block:: text

   Genome-wide association analyses were performed using the EDGE 
   (Elastic Data-Driven Encoding) approach [1], implemented in the 
   edge-gwas Python package version 0.1.0 [2]. EDGE uses a two-stage 
   analysis framework that calculates flexible genetic encoding 
   parameters (α) on training data and applies them to test data, 
   allowing detection of nonadditive genetic effects without assuming 
   a specific inheritance model. Genotype data were filtered for minor 
   allele frequency (MAF > 0.01) and missingness (<5%). The dataset 
   was split 50/50 into training and test sets using stratified sampling 
   to maintain case-control balance. Alpha values were calculated on the 
   training set, and GWAS was performed on the test set using logistic 
   regression with age, sex, and the first 10 principal components as 
   covariates. Genome-wide significance was defined as p < 5×10⁻⁸.
   
   [1] Zhou, J., et al. (2023). medRxiv, 2023.06.01.23290857.
   [2] Zhou, J., Hall, M.A., & Contributors (2025). edge-gwas v0.1.0.

Acknowledgments
---------------

If you use edge-gwas, please consider acknowledging:

.. code-block:: text

   This research was conducted using edge-gwas (Zhou et al., 2025), 
   a Python package for flexibly encoded genome-wide association studies.

For Grant Applications
----------------------

When citing edge-gwas in grant applications:

.. code-block:: text

   We will perform GWAS using the EDGE (Elastic Data-Driven Encoding) 
   methodology (Zhou et al., 2023; medRxiv 2023.06.01.23290857), which 
   has been shown to identify novel nonadditive genetic associations 
   missed by traditional additive models. We will use the edge-gwas 
   Python package (v0.1.0), an open-source implementation providing 
   efficient two-stage analysis, quality control, and visualization tools.

Contact for Citation Questions
-------------------------------

For questions about how to cite edge-gwas:

* **Research collaboration**: molly.hall@pennmedicine.upenn.edu
* **Software questions**: jyzhou@stanford.edu
* **General inquiries**: Open an issue on GitHub

See Also
--------

* :ref:`changelog` - Version history
* :ref:`contributing` - How to contribute
* GitHub repository - https://github.com/YOUR-USERNAME/edge-gwas

*Last updated: 2025-12-25 for edge-gwas v0.1.1*

*For questions or issues, visit:* https://github.com/nicenzhou/edge-gwas/issues
