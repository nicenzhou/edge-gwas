.. _statistical_model:

Statistical Model
=================

EDGE (Elastic Data-Driven Encoding) employs a flexible encoding approach to detect 
nonadditive genetic effects that are often missed by traditional additive GWAS models.

Overview
--------

EDGE uses a two-stage analysis framework:

1. **Training stage**: Calculate encoding parameters (α) from training data
2. **Test stage**: Apply α values to test data for association testing

This approach allows detection of various inheritance patterns on a continuous scale,
from recessive to dominant to over-dominant effects.

**Key Innovation:** Unlike standard GWAS that assumes additive inheritance (α = 0.5),
EDGE estimates α directly from the data, enabling detection of nonadditive effects.

Regression Model
----------------

Stage 1: Codominant Model (Training Data)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

EDGE separately estimates effects for heterozygous and homozygous alternate genotypes:

**Equation 1: Codominant Regression Model**

.. math::

   E(Y | SNP_{Het}, SNP_{HA}, COV_i) = \beta_0 + \beta_{Het} \cdot SNP_{Het} + \beta_{HA} \cdot SNP_{HA} + \sum_{i} \beta_{cov_i} \cdot COV_i

Where:

* :math:`Y` = phenotype/outcome
* :math:`SNP_{Het}` = indicator for heterozygous genotype (0 or 1)
* :math:`SNP_{HA}` = indicator for homozygous alternate genotype (0 or 1)
* :math:`COV_i` = covariates (age, sex, principal components, etc.)
* :math:`\beta_0` = intercept
* :math:`\beta_{Het}` = effect size for heterozygous genotype
* :math:`\beta_{HA}` = effect size for homozygous alternate genotype
* :math:`\beta_{cov_i}` = effect sizes for covariates

**Genotype Encoding for Equation 1:**

.. list-table::
   :header-rows: 1
   :widths: 25 25 25 25

   * - True Genotype
     - :math:`SNP_{Het}`
     - :math:`SNP_{HA}`
     - Interpretation
   * - 0/0 (Ref/Ref)
     - 0
     - 1
     - Reference homozygote
   * - 0/1 (Ref/Alt)
     - 1
     - 0
     - Heterozygote
   * - 1/1 (Alt/Alt)
     - 0
     - 0
     - Alternate homozygote

Encoding Parameter
------------------

**Equation 2: Alpha Calculation**

.. math::

   \alpha = \frac{\beta_{Het}}{\beta_{HA}}

Where:

* :math:`\alpha` = encoding parameter representing the ratio of heterozygous to homozygous alternate effects
* :math:`\beta_{Het}` = effect size for heterozygous genotype (from Equation 1)
* :math:`\beta_{HA}` = effect size for homozygous alternate genotype (from Equation 1)

**Special Cases:**

* If :math:`|\beta_{HA}| < 10^{-10}`, α is undefined (SNP skipped)
* If model fails to converge, SNP is skipped

Stage 2: EDGE-Encoded Model (Test Data)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Equation 3: EDGE-Encoded Genotypes**

.. math::

   G_{EDGE} = \begin{cases}
   1 & \text{if genotype is 0/0 (Ref/Ref)} \\
   \alpha & \text{if genotype is 0/1 (Het)} \\
   0 & \text{if genotype is 1/1 (Alt/Alt)}
   \end{cases}

**Equation 4: Association Testing Model**

.. math::

   E(Y | G_{EDGE}, COV_i) = \beta_0 + \beta \cdot G_{EDGE} + \sum_{i} \beta_{cov_i} \cdot COV_i

Where:

* :math:`G_{EDGE}` = EDGE-encoded genotype (from Equation 3)
* :math:`\beta` = effect coefficient being tested
* Null hypothesis: :math:`H_0: \beta = 0`

Interpretation of Alpha
-----------------------

The α parameter indicates the inheritance pattern:

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - Alpha Range
     - Inheritance Model
     - Biological Interpretation
   * - :math:`\alpha < 0`
     - Under-recessive
     - Heterozygotes have opposite effect direction
   * - :math:`\alpha \approx 0`
     - Recessive
     - Heterozygotes similar to reference homozygotes
   * - :math:`0 < \alpha < 0.3`
     - Partially recessive
     - Heterozygotes closer to reference
   * - :math:`\alpha \approx 0.5`
     - Additive
     - Heterozygotes have half the effect (standard GWAS)
   * - :math:`0.7 < \alpha < 1.3`
     - Partially dominant
     - Heterozygotes closer to alternate homozygotes
   * - :math:`\alpha \approx 1`
     - Dominant
     - Heterozygotes similar to alternate homozygotes
   * - :math:`\alpha > 1.3`
     - Over-dominant
     - Heterozygotes have stronger effect than both homozygotes

**Example Interpretations:**

* :math:`\alpha = 0`: Recessive disease allele (e.g., cystic fibrosis)
* :math:`\alpha = 0.5`: Additive effect (standard GWAS assumption)
* :math:`\alpha = 1`: Dominant disease allele (e.g., Huntington's)
* :math:`\alpha = 2`: Heterozygote advantage (e.g., sickle cell trait)

Two-Stage Analysis
------------------

Why Two Stages?
~~~~~~~~~~~~~~~

**Prevents Overfitting:**

* Using same data to estimate α and test association inflates false positives
* Two-stage approach provides unbiased p-values
* Independent test set ensures valid statistical inference

**Balances Power and Bias:**

* 50/50 split recommended for optimal balance
* Larger training set: more stable α estimates
* Larger test set: more powerful association testing

Stage 1: Calculate Alpha (Training Data)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For each SNP in the training dataset:

1. Fit codominant model (Equation 1) using:
   
   * Logistic regression for binary outcomes
   * Linear regression for continuous outcomes

2. Extract coefficients :math:`\beta_{Het}` and :math:`\beta_{HA}`

3. Calculate :math:`\alpha = \beta_{Het} / \beta_{HA}`

4. Store α value for this SNP

5. Check convergence:
   
   * If model converges: save α
   * If model fails: skip SNP and record in log

**Statistical Software:**

* Uses ``statsmodels`` package in Python
* BFGS optimization algorithm
* Maximum iterations: configurable (default: 1000)

Stage 2: Apply Alpha (Test Data)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For each SNP in the test dataset:

1. Retrieve α value calculated from training data

2. Encode genotypes using Equation 3:
   
   * 0/0 → 1
   * 0/1 → α
   * 1/1 → 0

3. Fit association model (Equation 4)

4. Test :math:`H_0: \beta = 0` using Wald test

5. Extract p-value, effect size, standard error

**Output:**

* Variant ID
* Chromosome and position
* P-value
* Effect coefficient (:math:`\beta`)
* Standard error
* Test statistic
* Applied α value

Outcome Types
-------------

Binary Outcomes (Logistic Regression)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For case-control studies (e.g., disease status):

**Training Stage:**

.. math::

   logit(P(Y=1)) = \beta_0 + \beta_{Het} \cdot SNP_{Het} + \beta_{HA} \cdot SNP_{HA} + \sum_{i} \beta_{cov_i} \cdot COV_i

**Test Stage:**

.. math::

   logit(P(Y=1)) = \beta_0 + \beta \cdot G_{EDGE} + \sum_{i} \beta_{cov_i} \cdot COV_i

Where :math:`P(Y=1)` is the probability of being a case.

**Effect Interpretation:**

* :math:`\beta > 0`: Risk allele (increases disease odds)
* :math:`\beta < 0`: Protective allele (decreases disease odds)
* Odds ratio: :math:`OR = e^\beta`

Quantitative Outcomes (Linear Regression)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For continuous traits (e.g., BMI, blood pressure):

**Training Stage:**

.. math::

   Y = \beta_0 + \beta_{Het} \cdot SNP_{Het} + \beta_{HA} \cdot SNP_{HA} + \sum_{i} \beta_{cov_i} \cdot COV_i + \epsilon

**Test Stage:**

.. math::

   Y = \beta_0 + \beta \cdot G_{EDGE} + \sum_{i} \beta_{cov_i} \cdot COV_i + \epsilon

Where :math:`\epsilon \sim N(0, \sigma^2)` is the error term.

**Effect Interpretation:**

* :math:`\beta` = change in trait per unit change in :math:`G_{EDGE}`
* Units depend on trait (e.g., kg for weight, mmHg for blood pressure)

Advantages of EDGE
------------------

Compared to Standard Additive GWAS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Flexible inheritance**: Detects recessive, dominant, and over-dominant effects
2. **Data-driven**: α estimated from data, not assumed to be 0.5
3. **More powerful**: Can detect associations missed by additive models
4. **Interpretable**: α values provide biological insights
5. **Robust**: Two-stage approach prevents overfitting

Compared to Other Nonadditive Methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

6. **Continuous spectrum**: Tests all inheritance patterns, not discrete models
7. **Single parameter**: α summarizes inheritance in one value
8. **Computational efficiency**: Faster than testing multiple genetic models
9. **No model selection**: Avoids multiple testing penalty from model selection

Statistical Properties
~~~~~~~~~~~~~~~~~~~~~~~

10. **Unbiased p-values**: Two-stage design ensures valid inference
11. **Consistent estimates**: α converges to true value with sufficient sample size
12. **Type I error control**: Maintains nominal false positive rate

Statistical Testing
-------------------

Hypothesis Testing
~~~~~~~~~~~~~~~~~~

**Null Hypothesis:**

.. math::

   H_0: \beta = 0 \text{ (no association between SNP and outcome)}

**Alternative Hypothesis:**

.. math::

   H_a: \beta \neq 0 \text{ (SNP is associated with outcome)}

**Test Statistic:**

Wald test statistic:

.. math::

   Z = \frac{\hat{\beta}}{SE(\hat{\beta})}

Where:

* :math:`\hat{\beta}` = estimated coefficient
* :math:`SE(\hat{\beta})` = standard error of the coefficient

P-value Calculation
~~~~~~~~~~~~~~~~~~~

**Two-tailed test:**

.. math::

   p = 2 \times P(|Z| > |z_{obs}|)

Where :math:`z_{obs}` is the observed test statistic.

**For logistic regression:**

* Test statistic follows standard normal distribution asymptotically
* P-value calculated from normal distribution

**For linear regression:**

* Can use t-distribution with :math:`n - p - 1` degrees of freedom
* Large samples: t-distribution ≈ normal distribution

Multiple Testing Correction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Genome-wide Significance:**

.. math::

   p < 5 \times 10^{-8}

This threshold corresponds to Bonferroni correction for ~1 million independent tests.

**Suggestive Threshold:**

.. math::

   p < 1 \times 10^{-5}

Used for exploratory analysis or hypothesis generation.

**Other Methods:**

* **Bonferroni correction**: :math:`\alpha_{corrected} = \frac{0.05}{n_{tests}}`
* **False Discovery Rate (FDR)**: Controls proportion of false discoveries
* **Genomic control (λ)**: Adjusts for population stratification

Genomic Inflation Factor
~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

   \lambda = \frac{median(\chi^2_{obs})}{median(\chi^2_{expected})} = \frac{median(\chi^2_{obs})}{0.455}

**Interpretation:**

* :math:`\lambda \approx 1.0`: No inflation (good)
* :math:`\lambda > 1.05`: Possible population stratification or cryptic relatedness
* :math:`\lambda < 0.95`: Possible over-correction

**Remedies for inflation:**

* Add more principal components as covariates
* Check for batch effects
* Verify sample quality control
* Consider mixed models for related individuals

Quality Control Considerations
-------------------------------

Minor Allele Frequency (MAF)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Recommendation:** MAF > 0.01 (1%)

**Rationale:**

* Rare variants have unstable α estimates
* Insufficient heterozygotes for reliable :math:`\beta_{Het}` estimation
* Higher chance of convergence failure

**For rare variant analysis:**

* Use gene-based or region-based methods
* Increase MAF threshold to 0.05 for small samples (<5,000)
* Consider burden tests or SKAT

Hardy-Weinberg Equilibrium (HWE)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Recommendation:** HWE p > 1 × 10⁻⁶ in controls

**Calculation:**

.. math::

   \chi^2_{HWE} = \frac{(O_{Het} - E_{Het})^2}{E_{Het}}

Where:

* :math:`O_{Het}` = observed heterozygote count
* :math:`E_{Het} = 2 \times n \times p \times (1-p)`
* :math:`n` = sample size
* :math:`p` = allele frequency

**Deviation indicates:**

* Genotyping errors
* Population stratification
* Selection pressure
* Copy number variation

Missingness
~~~~~~~~~~~

**Per-variant recommendation:** < 5% missing

**Per-sample recommendation:** < 5% missing

**Impact on EDGE:**

* Missing genotypes excluded from model fitting
* Reduces effective sample size
* May bias α estimates if missingness is non-random

Population Structure
~~~~~~~~~~~~~~~~~~~~

**Recommendation:** Include 10+ principal components (PCs) as covariates

**Why PCs matter:**

* Control for population stratification
* Reduce confounding
* Lower genomic inflation
* Improve power by reducing variance

**How many PCs:**

* Start with 10 PCs
* Check :math:`\lambda`; if > 1.05, add more (up to 20)
* Diminishing returns beyond 20 PCs
* Too many PCs can reduce power

Sample Size Considerations
---------------------------

Minimum Sample Size
~~~~~~~~~~~~~~~~~~~

**For binary outcomes:**

* Minimum: 500 cases + 500 controls
* Recommended: 5,000+ cases and controls
* Power depends on effect size, MAF, and α

**For continuous outcomes:**

* Minimum: 1,000 samples
* Recommended: 10,000+ samples

Power Calculations
~~~~~~~~~~~~~~~~~~

Power to detect association depends on:

1. Sample size (larger = more power)
2. Effect size (larger = more power)
3. MAF (intermediate MAF = most power)
4. α value (closer to 0.5 = less advantage over additive)
5. Significance threshold (more stringent = less power)

**EDGE has most advantage when:**

* True α ≠ 0.5 (nonadditive inheritance)
* Sufficient heterozygotes (MAF 0.05-0.50)
* Large sample size (N > 5,000)
* True effect size is moderate to large

Comparison with Other Methods
------------------------------

EDGE vs. Standard Additive GWAS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Standard Additive Model:**

.. math::

   E(Y | G) = \beta_0 + \beta \cdot G_{additive} + \sum_{i} \beta_{cov_i} \cdot COV_i

Where :math:`G_{additive}` is the allele count (0, 1, or 2).

**Key Differences:**

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Aspect
     - Additive GWAS
     - EDGE
   * - Genetic model
     - Fixed (α = 0.5)
     - Data-driven (α estimated)
   * - Inheritance patterns
     - Additive only
     - All patterns (recessive to dominant)
   * - Power for additive effects
     - Optimal
     - Comparable
   * - Power for nonadditive effects
     - Low
     - High
   * - Computational cost
     - Lower
     - Moderate (2-stage)
   * - Interpretation
     - Simple
     - More informative (α values)

EDGE vs. Genotypic Model
~~~~~~~~~~~~~~~~~~~~~~~~~

**Genotypic Model:**

Tests both :math:`\beta_{Het}` and :math:`\beta_{HA}` simultaneously using 2 df test.

**Advantages of EDGE:**

* 1 df test (more powerful with correct α)
* Provides interpretable α parameter
* Better for replication (apply α to new data)

**Disadvantages of EDGE:**

* Requires data splitting
* Less powerful if α is unstable

EDGE vs. Cochran-Armitage Trend Test
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Trend Test with weights:**

.. math::

   \text{weights} = (0, w, 1) \text{ for genotypes } (0/0, 0/1, 1/1)

**Relationship to EDGE:**

* EDGE α is analogous to weight w
* Standard additive: w = 0.5
* Recessive: w = 0
* Dominant: w = 1

**EDGE advantages:**

* α estimated from data, not pre-specified
* Two-stage design prevents overfitting

Practical Considerations
-------------------------

When to Use EDGE
~~~~~~~~~~~~~~~~

**EDGE is most beneficial when:**

1. Suspected nonadditive inheritance (family studies suggest recessive/dominant)
2. Previous additive GWAS found few/no associations
3. Biological mechanism suggests specific inheritance pattern
4. Follow-up of candidate genes with known inheritance
5. Large sample size available (>5,000 samples)

**Additive GWAS may be sufficient when:**

1. Initial exploratory analysis
2. Small sample size (<1,000 samples)
3. Known additive effects from literature
4. Computational resources limited

Recommended Workflow
~~~~~~~~~~~~~~~~~~~~

**Step 1: Run both EDGE and additive GWAS**

.. code-block:: python

   from edge_gwas import EDGEAnalysis
   from edge_gwas.utils import additive_gwas
   
   # EDGE analysis
   edge = EDGEAnalysis(outcome_type='binary')
   alpha_df, edge_results = edge.run_full_analysis(
       train_g, train_p, test_g, test_p, outcome, covariates
   )
   
   # Additive GWAS for comparison
   additive_results = additive_gwas(
       test_g, test_p, outcome, covariates, outcome_type='binary'
   )

**Step 2: Compare results**

* Identify SNPs significant in EDGE but not additive
* Examine α values for significant SNPs
* Check for biological plausibility

**Step 3: Validate findings**

* Replicate in independent cohort
* Check α consistency across cohorts
* Perform functional follow-up

Computational Complexity
~~~~~~~~~~~~~~~~~~~~~~~~

**Time Complexity:**

* Training stage: :math:`O(n_{train} \times m \times p^2 \times k)`
* Test stage: :math:`O(n_{test} \times m \times p^2)`

Where:

* :math:`n` = sample size
* :math:`m` = number of variants
* :math:`p` = number of covariates
* :math:`k` = iterations for convergence

**Space Complexity:**

* :math:`O(n \times m)` for genotype data
* Additional :math:`O(m)` for α values storage

**Optimization Strategies:**

* Process chromosomes separately
* Use parallel processing (multiple cores)
* Filter variants before analysis
* Use sparse matrix representations for rare variants

Implementation Details
----------------------

Software Implementation
~~~~~~~~~~~~~~~~~~~~~~~

**edge-gwas v0.1.1** uses:

* **statsmodels**: Logistic and linear regression
* **scikit-learn**: Data splitting, preprocessing
* **pandas**: Data management
* **numpy/scipy**: Numerical computations
* **joblib**: Parallel processing

**Numerical Stability:**

* Check for :math:`|\beta_{HA}| > 10^{-10}` before computing α
* Use BFGS optimization with line search
* Employ convergence tolerance of 1e-8
* Maximum iterations: configurable (default: 1000)

Convergence Issues
~~~~~~~~~~~~~~~~~~

**SNPs may fail to converge due to:**

1. **Rare variants**: Too few observations for stable estimates
2. **Perfect separation**: In logistic regression when outcome completely separated by genotype
3. **Multicollinearity**: Highly correlated covariates
4. **Numerical instability**: Extreme coefficient values

**Solutions:**

* Increase ``max_iter`` parameter
* Filter rare variants (MAF > 0.01)
* Check covariate correlations
* Use ``get_skipped_snps()`` to identify problematic variants

Extensions and Future Directions
---------------------------------

Planned Features (v0.2.0+)
~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Gene-based testing**: Aggregate SNPs within genes
2. **Mixed models**: Account for relatedness and population structure
3. **X chromosome**: Proper handling of X-linked inheritance
4. **Rare variant analysis**: Optimized methods for MAF < 0.01
5. **Meta-analysis tools**: Combine α estimates across studies
6. **GPU acceleration**: Faster computation for biobank-scale data

Research Directions
~~~~~~~~~~~~~~~~~~~

* **Optimal splitting ratio**: Beyond 50/50 for different scenarios
* **Cross-validation**: Stability of α across folds
* **Bayesian EDGE**: Prior distributions on α
* **Multi-trait EDGE**: Joint analysis of related phenotypes
* **Gene-environment interaction**: EDGE with interaction terms

Limitations
-----------

Current Limitations
~~~~~~~~~~~~~~~~~~~

1. **Data splitting required**: Reduces effective sample size
2. **Assumes independence**: Not designed for related individuals
3. **Single SNP analysis**: Does not model joint effects
4. **Autosomal variants only**: X/Y chromosomes require special handling
5. **European-ancestry optimized**: May need adjustment for diverse populations

Statistical Assumptions
~~~~~~~~~~~~~~~~~~~~~~~

EDGE assumes:

1. **Independent samples**: No cryptic relatedness
2. **No genotyping errors**: QC filters applied
3. **Hardy-Weinberg Equilibrium**: In controls for binary outcomes
4. **Correct model specification**: Covariates adequately control confounding
5. **Large sample size**: Asymptotic properties of tests hold

Validation Studies
------------------

Simulation Studies
~~~~~~~~~~~~~~~~~~

EDGE has been validated through simulations showing:

* **Type I error control**: Maintains nominal false positive rate (5%)
* **Increased power**: 10-50% power gain for nonadditive effects
* **Robustness**: Stable performance across MAF ranges
* **α recovery**: Accurately estimates true inheritance pattern

Real Data Applications
~~~~~~~~~~~~~~~~~~~~~~

EDGE has been applied to:

* **Cardiometabolic traits**: BMI, diabetes, blood pressure
* **UK Biobank data**: >500,000 individuals
* **Novel discoveries**: Identified nonadditive effects missed by additive GWAS

See **Zhou et al. (2023)** for detailed results.

Mathematical Proofs
-------------------

Unbiased p-values (Sketch)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Theorem**: Two-stage EDGE produces valid p-values under the null hypothesis.

**Proof sketch**:

1. α estimated on training data :math:`D_{train}`
2. Association tested on independent test data :math:`D_{test}`
3. Under :math:`H_0`: :math:`\beta = 0`, test statistic in :math:`D_{test}` follows null distribution
4. Independence of :math:`D_{train}` and :math:`D_{test}` ensures no inflation
5. Therefore, Type I error rate is controlled at nominal level

Consistency of α (Sketch)
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Theorem**: :math:`\hat{\alpha} \xrightarrow{p} \alpha_{true}` as :math:`n \rightarrow \infty`

**Conditions**:

* :math:`|\beta_{HA}| > c > 0` for some constant c
* Standard regularity conditions for MLE

**Implication**: With sufficient training data, α estimates converge to true value.

Glossary of Terms
-----------------

.. glossary::

   Alpha (α)
      Encoding parameter representing the ratio of heterozygous to homozygous effects
   
   Additive model
      Genetic model assuming heterozygote effect is exactly half of homozygous effect (α = 0.5)
   
   BFGS
      Broyden-Fletcher-Goldfarb-Shanno algorithm for optimization
   
   Codominant model
      Model that separately estimates effects for heterozygotes and homozygotes
   
   EDGE encoding
      Genotype encoding using estimated α values
   
   Genomic inflation (λ)
      Measure of systematic p-value inflation due to population structure
   
   HWE
      Hardy-Weinberg Equilibrium
   
   MAF
      Minor Allele Frequency
   
   Over-dominant
      Inheritance pattern where heterozygotes have stronger effect than either homozygote (α > 1)
   
   Recessive
      Inheritance pattern where heterozygotes behave like reference homozygotes (α ≈ 0)
   
   Wald test
      Statistical test using ratio of estimate to standard error
   
   Two-stage analysis
      Analysis approach separating parameter estimation (stage 1) from hypothesis testing (stage 2)

References
----------

Primary Citation
~~~~~~~~~~~~~~~~

**Zhou J, Rico ALG, Guare L, et al.** (2025). Flexibly encoded genome-wide association study identifies 
novel nonadditive genetic risk variants for cardiometabolic traits. 
*medRxiv*, 2023.06.01.23290857. 
https://doi.org/10.1101/2023.06.01.23290857

See Also
--------

* :ref:`installation` - Installation guide
* :ref:`quickstart` - Quick start tutorial
* :ref:`user_guide` - Comprehensive user guide
* :ref:`examples` - Example analyses with code
* :ref:`api_reference` - Complete function documentation
* :ref:`visualization` - Plotting and visualization

---

*Last updated: 2025-12-25 for edge-gwas v0.1.1*

*For questions or issues, visit:* https://github.com/nicenzhou/edge-gwas/issues
