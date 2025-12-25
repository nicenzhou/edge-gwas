.. _statistical_model:

Statistical Model
=================

EDGE (Elastic Data-Driven Encoding) employs a flexible encoding approach to detect 
nonadditive genetic effects.

Overview
--------

EDGE uses a two-stage analysis:

1. **Training stage**: Calculate encoding parameters (α) from training data
2. **Test stage**: Apply α values to test data for GWAS

This approach allows detection of various inheritance patterns on a continuous scale.

Regression Model
----------------

EDGE separately estimates effects for heterozygous and homozygous alternate genotypes:

**Equation 1: Regression Model**

.. math::

   E(Y | SNP_{Het}, SNP_{HA}, COV_i) = \\beta_0 + \\beta_{Het} \\cdot SNP_{Het} + \\beta_{HA} \\cdot SNP_{HA} + \\sum_{i} \\beta_{cov_i} \\cdot COV_i

Where:

* :math:`Y` = phenotype/outcome
* :math:`SNP_{Het}` =
