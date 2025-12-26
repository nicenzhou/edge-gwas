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

   E(Y | SNP_{Het}, SNP_{HA}, COV_i) = \beta_0 + \beta_{Het} \cdot SNP_{Het} + \beta_{HA} \cdot SNP_{HA} + \sum_{i} \beta_{cov_i} \cdot COV_i

Where:

* :math:`Y` = phenotype/outcome
* :math:`SNP_{Het}` = indicator for heterozygous genotype
* :math:`SNP_{HA}` = indicator for homozygous alternate genotype
* :math:`COV_i` = covariates (age, sex, principal components, etc.)
* :math:`\beta_{Het}` = effect size for heterozygous genotype
* :math:`\beta_{HA}` = effect size for homozygous alternate genotype
* :math:`\beta_{cov_i}` = effect sizes for covariates

Encoding Parameter
------------------

**Equation 2: Alpha Calculation**

.. math::

   \alpha = \frac{\beta_{Het}}{\beta_{HA}}

Where:

* :math:`\alpha` = encoding parameter representing the ratio of heterozygous to homozygous alternate effects
* :math:`\beta_{Het}` = effect size for heterozygous genotype (from Equation 1)
* :math:`\beta_{HA}` = effect size for homozygous alternate genotype (from Equation 1)

Interpretation of Alpha
-----------------------

The α parameter indicates the inheritance pattern:

* **α ≈ 0**: Recessive model (heterozygotes similar to reference homozygotes)
* **α ≈ 0.5**: Additive model (heterozygotes have half the effect)
* **α ≈ 1**: Dominant model (heterozygotes similar to alternate homozygotes)
* **α < 0**: Under-recessive (heterozygotes have opposite effect)
* **α > 1**: Over-dominant (heterozygotes have stronger effect than homozygotes)

Two-Stage Analysis
------------------

Stage 1: Calculate Alpha (Training Data)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. For each SNP, fit Equation 1 on training data
2. Extract :math:`\beta_{Het}` and :math:`\beta_{HA}`
3. Calculate :math:`\alpha = \beta_{Het} / \beta_{HA}`
4. Store α values for each SNP

Stage 2: Apply Alpha (Test Data)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. For each SNP, use the α from training data
2. Encode genotypes: :math:`G_{encoded} = \alpha \cdot G_{Het} + G_{HA}`
3. Fit standard GWAS model:

.. math::

   E(Y | G_{encoded}, COV_i) = \beta_0 + \beta \cdot G_{encoded} + \sum_{i} \beta_{cov_i} \cdot COV_i

4. Test significance of :math:`\beta`

Outcome Types
-------------

Binary Outcomes (Logistic Regression)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For case-control studies:

.. math::

   logit(P(Y=1)) = \beta_0 + \beta_{Het} \cdot SNP_{Het} + \beta_{HA} \cdot SNP_{HA} + \sum_{i} \beta_{cov_i} \cdot COV_i

Quantitative Outcomes (Linear Regression)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For continuous traits:

.. math::

   Y = \beta_0 + \beta_{Het} \cdot SNP_{Het} + \beta_{HA} \cdot SNP_{HA} + \sum_{i} \beta_{cov_i} \cdot COV_i + \epsilon

Where :math:`\epsilon \sim N(0, \sigma^2)`

Advantages of EDGE
------------------

1. **Flexible inheritance**: Detects nonadditive effects without assuming a specific model
2. **Data-driven**: α is calculated from the data, not pre-specified
3. **Robust**: Two-stage approach reduces overfitting
4. **Interpretable**: α values have clear biological meaning
5. **Powerful**: Can detect associations missed by additive models

Statistical Testing
-------------------

Hypothesis Testing
~~~~~~~~~~~~~~~~~~

Null hypothesis: :math:`H_0: \beta = 0` (no association)

Alternative: :math:`H_a: \beta \neq 0` (SNP associated with outcome)

P-value Calculation
~~~~~~~~~~~~~~~~~~~

P-values are calculated using:

* Wald test for coefficients
* Chi-square distribution for test statistics
* Appropriate degrees of freedom based on model

Multiple Testing Correction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Genome-wide significance threshold:

* **Standard**: p < 5 × 10⁻⁸
* **Suggestive**: p < 1 × 10⁻⁵

Consider using:

* Bonferroni correction
* False Discovery Rate (FDR)
* Genomic control (λ adjustment)

Quality Control Considerations
-------------------------------

Minor Allele Frequency (MAF)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Recommend: MAF > 0.01 (1%)
* Rare variants may have unstable α estimates

Hardy-Weinberg Equilibrium (HWE)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Test in controls: HWE p > 1 × 10⁻⁶
* Deviation may indicate genotyping errors

Missingness
~~~~~~~~~~~

* Per-variant: < 5% missing
* Per-sample: < 5% missing

Population Structure
~~~~~~~~~~~~~~~~~~~~

* Include principal components (PCs) as covariates
* Typically 10 PCs recommended
* Check genomic inflation factor (λ)

References
----------

Zhou, J., et al. (2023). Flexibly encoded genome-wide association study identifies 
novel nonadditive genetic risk variants for cardiometabolic traits. 
*medRxiv*, 2023.06.01.23290857. 
https://doi.org/10.1101/2023.06.01.23290857

See Also
--------

* :ref:`quickstart` - Quick start guide
* :ref:`examples` - Example analyses
* :ref:`api_reference` - Function documentation
