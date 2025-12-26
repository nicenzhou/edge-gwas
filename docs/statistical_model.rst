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
   
   * Logistic regression for binary outcomes (BFGS optimization)
   * Linear regression for continuous outcomes (configurable optimization method)

2. Extract coefficients :math:`\beta_{Het}` and :math:`\beta_{HA}`

3. Calculate :math:`\alpha = \beta_{Het} / \beta_{HA}`

4. Store α value for this SNP

5. Check convergence:
   
   * If model converges: save α
   * If model fails: skip SNP and record in log

**Statistical Software:**

* Uses ``statsmodels`` package in Python
* Optimization algorithms (see Optimization Methods section)
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

**Optimization:**

* Uses BFGS algorithm by default
* Suitable for non-convex optimization problems
* Handles perfect separation robustly

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

**Optimization:**

* Configurable optimization method (default: BFGS)
* See Optimization Methods section for available algorithms

Optimization Methods
--------------------

For continuous outcomes (linear regression), users can select from multiple optimization algorithms:

Available Algorithms
~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Method
     - Description
     - Best Use Case
   * - ``'bfgs'``
     - **Broyden-Fletcher-Goldfarb-Shannon** (default)
     - General purpose, good balance of speed and accuracy
   * - ``'newton'``
     - **Newton-Raphson**
     - Fast convergence when close to optimum
   * - ``'lbfgs'``
     - **Limited-memory BFGS**
     - Large datasets, memory-constrained systems
   * - ``'nm'``
     - **Nelder-Mead**
     - Derivative-free, robust to noise
   * - ``'cg'``
     - **Conjugate Gradient**
     - Large-scale problems, sparse matrices
   * - ``'ncg'``
     - **Newton Conjugate Gradient**
     - Large problems requiring second-order information
   * - ``'powell'``
     - **Powell's method**
     - Derivative-free, quadratic convergence
   * - ``'basinhopping'``
     - **Basin-hopping**
     - Global optimization, multiple local minima

Algorithm Details
~~~~~~~~~~~~~~~~~

**BFGS (Default):**

* Quasi-Newton method using approximation of Hessian matrix
* Good balance between computational cost and convergence speed
* **Pros:** Fast, memory-efficient, handles ill-conditioned problems
* **Cons:** May struggle with very large datasets
* **Recommended for:** Most GWAS applications

**Newton-Raphson:**

* Uses exact Hessian matrix (second derivatives)
* **Pros:** Quadratic convergence near optimum
* **Cons:** Expensive to compute Hessian, requires good starting values
* **Recommended for:** Small datasets, when high precision needed

**L-BFGS:**

* Limited-memory version of BFGS
* **Pros:** Handles large-scale problems, low memory footprint
* **Cons:** Slightly slower convergence than BFGS
* **Recommended for:** Biobank-scale data (>100,000 samples)

**Nelder-Mead:**

* Simplex-based, derivative-free method
* **Pros:** Robust, doesn't require gradients
* **Cons:** Slower convergence, not suitable for high dimensions
* **Recommended for:** Debugging, noisy data

**Conjugate Gradient:**

* Iterative method using conjugate directions
* **Pros:** Memory-efficient, good for sparse problems
* **Cons:** Slower than Newton methods
* **Recommended for:** Very large covariate sets

**Newton Conjugate Gradient:**

* Combines Newton method with CG for Hessian inversion
* **Pros:** Handles large-scale problems efficiently
* **Cons:** Requires gradient and Hessian-vector products
* **Recommended for:** Problems with >1000 covariates

**Powell's Method:**

* Derivative-free, uses conjugate directions
* **Pros:** Robust, no gradient calculation
* **Cons:** Can be slow, may not converge
* **Recommended for:** Non-smooth objective functions

**Basin-hopping:**

* Global optimization combining local search with random jumps
* **Pros:** Finds global optimum, robust to starting values
* **Cons:** Very slow, many function evaluations
* **Recommended for:** Suspected multiple local minima (rare in GWAS)

Usage Example
~~~~~~~~~~~~~

.. code-block:: python

   from edge_gwas import EDGEAnalysis
   
   # Default: BFGS optimization
   edge = EDGEAnalysis(
       outcome_type='continuous',
       ols_method='bfgs'  # Default, can be omitted
   )
   
   # Use L-BFGS for large dataset
   edge_large = EDGEAnalysis(
       outcome_type='continuous',
       ols_method='lbfgs'
   )
   
   # Use Newton for high precision
   edge_precise = EDGEAnalysis(
       outcome_type='continuous',
       ols_method='newton',
       max_iter=2000  # May need more iterations
   )
   
   # Calculate alpha values
   alpha_df = edge.calculate_alpha(
       train_genotype,
       train_phenotype,
       outcome='bmi',
       covariates=['age', 'sex', 'pc1', 'pc2']
   )

Choosing an Optimization Method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Decision Tree:**

1. **Start with default (BFGS)**
   
   * Works well for >95% of cases
   * Good balance of speed and accuracy

2. **If convergence issues occur:**
   
   * Try ``'newton'`` for better precision
   * Try ``'nm'`` or ``'powell'`` for robustness
   * Increase ``max_iter`` parameter

3. **For very large datasets (N > 100,000):**
   
   * Use ``'lbfgs'`` to reduce memory usage
   * Consider ``'cg'`` if memory is still an issue

4. **For many covariates (P > 100):**
   
   * Use ``'ncg'`` for efficiency
   * Consider ``'cg'`` as alternative

5. **If suspect multiple solutions:**
   
   * Use ``'basinhopping'`` (very slow!)
   * Check results for biological plausibility

**Performance Comparison (approximate):**

.. list-table::
   :header-rows: 1
   :widths: 20 20 20 20 20

   * - Method
     - Speed
     - Memory
     - Accuracy
     - Robustness
   * - BFGS
     - Fast
     - Medium
     - High
     - High
   * - Newton
     - Very Fast
     - High
     - Very High
     - Medium
   * - L-BFGS
     - Fast
     - Low
     - High
     - High
   * - Nelder-Mead
     - Slow
     - Low
     - Medium
     - Very High
   * - CG
     - Medium
     - Low
     - Medium
     - Medium
   * - NCG
     - Medium
     - Medium
     - High
     - Medium
   * - Powell
     - Slow
     - Low
     - Medium
     - High
   * - Basin-hopping
     - Very Slow
     - Medium
     - Very High
     - Very High

Convergence Diagnostics
~~~~~~~~~~~~~~~~~~~~~~~~

When a model fails to converge, the variant is automatically skipped and logged. You can:

1. **Check skipped variants:**

   .. code-block:: python

      skipped = edge.get_skipped_snps()
      print(f"Number of skipped SNPs: {len(skipped)}")

2. **Try different optimization method:**

   .. code-block:: python

      edge_robust = EDGEAnalysis(
          outcome_type='continuous',
          ols_method='nm',  # More robust
          max_iter=5000     # More iterations
      )

3. **Examine specific failing variants:**

   * Check MAF (may be too low)
   * Check for perfect collinearity with covariates
   * Verify no data quality issues

**Warning Signs:**

* >10% of SNPs fail to converge → Check data quality
* Consistent failures for common variants → Try different method
* Only rare variants fail → Expected, may increase MAF threshold

Outcome Transformations
-----------------------

For continuous outcomes, EDGE supports several transformations to improve normality and stabilize variance.

Available Transformations
~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 35 40

   * - Transformation
     - Mathematical Form
     - Use Case
   * - None (default)
     - :math:`Y' = Y`
     - Normally distributed outcomes
   * - ``'log'``
     - :math:`Y' = \ln(Y)`
     - Right-skewed, positive values only
   * - ``'log10'``
     - :math:`Y' = \log_{10}(Y)`
     - Right-skewed, prefer base-10 scale
   * - ``'inverse_normal'``
     - Parametric INT
     - Approximately normal, reduces outliers
   * - ``'rank_inverse_normal'``
     - Rank-based INT (RINT)
     - Heavy-tailed, robust to outliers

Log Transformation
~~~~~~~~~~~~~~~~~~

**Natural Log:**

.. math::

   Y' = \ln(Y)

**Log Base 10:**

.. math::

   Y' = \log_{10}(Y)

**Requirements:**

* All values must be positive (:math:`Y > 0`)
* Raises error if any :math:`Y \leq 0`

**Effect on interpretation:**

* :math:`\beta` represents proportional change
* :math:`e^\beta` is multiplicative effect (for natural log)
* Common for: income, biomarkers, gene expression

**Example:**

.. code-block:: python

   edge = EDGEAnalysis(
       outcome_type='continuous',
       outcome_transform='log'
   )
   
   # For outcome = triglycerides (right-skewed)
   alpha_df = edge.calculate_alpha(
       train_genotype,
       train_phenotype,
       outcome='triglycerides',  # Must be > 0
       covariates=['age', 'sex', 'bmi']
   )

Inverse Normal Transformation (INT)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Parametric INT:**

.. math::

   Y' = \Phi^{-1}\left(\frac{rank(Y) - 0.5}{n}\right)

Where :math:`\Phi^{-1}` is the inverse normal CDF.

**Process:**

1. Standardize Y using sample mean and SD
2. Calculate ranks
3. Convert ranks to quantiles
4. Apply inverse normal transformation

**Properties:**

* Assumes underlying normal distribution
* Less robust to outliers than RINT
* Preserves relative ordering

**Example:**

.. code-block:: python

   edge = EDGEAnalysis(
       outcome_type='continuous',
       outcome_transform='inverse_normal'
   )

Rank-Based Inverse Normal Transformation (RINT)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**RINT with Blom's Formula:**

.. math::

   Y' = \Phi^{-1}\left(\frac{rank(Y) - 3/8}{n + 1/4}\right)

**Alternative formulas:**

* **Van der Waerden:** :math:`\frac{rank(Y)}{n + 1}`
* **Blom (default):** :math:`\frac{rank(Y) - 3/8}{n + 1/4}`
* **Tukey:** :math:`\frac{rank(Y) - 1/3}{n + 1/3}`

**Advantages:**

* Extremely robust to outliers
* No assumptions about original distribution
* Handles ties appropriately (average rank method)
* Forces exact normal distribution

**Disadvantages:**

* Loses information about original scale
* Interpretation in normalized units only
* May reduce power for truly normal traits

**When to use:**

* Heavy-tailed distributions
* Presence of outliers
* Unknown original distribution
* Standard practice in many GWAS studies

**Example:**

.. code-block:: python

   edge = EDGEAnalysis(
       outcome_type='continuous',
       outcome_transform='rank_inverse_normal'
   )
   
   # For outcome with outliers (e.g., C-reactive protein)
   alpha_df = edge.calculate_alpha(
       train_genotype,
       train_phenotype,
       outcome='crp',  # Often has outliers
       covariates=['age', 'sex', 'bmi', 'smoking']
   )

Choosing a Transformation
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Decision Guide:**

1. **Check outcome distribution:**

   .. code-block:: python

      import matplotlib.pyplot as plt
      
      plt.hist(phenotype_df['outcome'], bins=50)
      plt.title('Outcome Distribution')
      plt.show()

2. **Apply transformation based on distribution:**

   * **Normal:** No transformation needed
   * **Right-skewed (e.g., income, biomarkers):** ``'log'`` or ``'log10'``
   * **Heavy-tailed with outliers:** ``'rank_inverse_normal'``
   * **Moderate deviation from normal:** ``'inverse_normal'``

3. **Verify transformation effectiveness:**

   .. code-block:: python

      from scipy import stats
      
      # Original
      _, p_original = stats.shapiro(y_original)
      
      # After transformation  
      _, p_transformed = stats.shapiro(y_transformed)
      
      print(f"Normality p-value: {p_original:.4f} → {p_transformed:.4f}")

**Practical Recommendations:**

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Trait Type
     - Recommended Transformation
   * - Anthropometric (height, weight)
     - None or ``'inverse_normal'``
   * - Lipids (cholesterol, triglycerides)
     - ``'log'`` or ``'rank_inverse_normal'``
   * - Inflammatory markers (CRP, IL-6)
     - ``'log'`` or ``'rank_inverse_normal'``
   * - Blood pressure
     - None or ``'inverse_normal'``
   * - Glucose, insulin
     - ``'log'``
   * - Gene expression
     - ``'log'`` or ``'rank_inverse_normal'``
   * - Counts (cell counts)
     - ``'log'`` or Poisson model
   * - Unknown distribution
     - ``'rank_inverse_normal'`` (safest)

Effect Size Interpretation After Transformation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Log transformation:**

* :math:`\beta` = log-fold change per EDGE unit
* Original scale: multiply by :math:`e^\beta` (natural log) or :math:`10^\beta` (log10)

**Inverse normal transformations:**

* :math:`\beta` in standard deviation units
* Not directly interpretable on original scale
* Focus on p-values and relative rankings

**Example interpretation:**

.. code-block:: python

   # Log-transformed triglycerides
   # beta = 0.05, alpha = 0.3
   # Interpretation: 
   # - Each EDGE-encoded allele associated with exp(0.05) = 1.05x
   # - 5% increase in triglycerides
   # - Partially recessive (alpha = 0.3)

**Reporting transformed results:**

* Always state transformation method used
* Report effects in transformed units
* Optionally back-transform for key findings
* Note limitations of interpretation

Transformation Diagnostics
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Check for invalid values:**

.. code-block:: python

   # After transformation
   print(f"NaN values: {y_transformed.isna().sum()}")
   print(f"Inf values: {np.isinf(y_transformed).sum()}")
   
   # Summary statistics
   print(y_transformed.describe())

**Common issues:**

* **Log of negative/zero:** Add small constant or use log1p
* **Perfect ties in RINT:** Rare, but can occur
* **Extreme outliers:** May still influence log transformation

**Solutions:**

.. code-block:: python

   # For log transformation with zeros
   y_transformed = np.log1p(y)  # log(1 + y)
   
   # For negative values, shift to positive
   y_shifted = y - y.min() + 1
   y_transformed = np.log(y_shifted)
   
   # Winsorize extreme outliers before transformation
   from scipy.stats import mstats
   y_winsorized = mstats.winsorize(y, limits=[0.01, 0.01])

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
   * - Optimization
     - Standard OLS/GLM
     - Configurable algorithms
   * - Outcome transformations
     - User must handle
     - Built-in support

EDGE vs. Genotypic Model
~~~~~~~~~~~~~~~~~~~~~~~~~

**Genotypic Model:**

Tests both :math:`\beta_{Het}` and :math:`\beta_{HA}` simultaneously using 2 df test.

**Advantages of EDGE:**

* 1 df test (more powerful with correct α)
* Provides interpretable α parameter
* Better for replication (apply α to new data)
* Optimized computational methods

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
* Flexible optimization for stability
* Built-in outcome transformations

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
6. Outcome is non-normal (use transformations)
7. Need interpretable inheritance patterns

**Additive GWAS may be sufficient when:**

1. Initial exploratory analysis
2. Small sample size (<1,000 samples)
3. Known additive effects from literature
4. Computational resources limited
5. Outcome is normally distributed

Recommended Workflow
~~~~~~~~~~~~~~~~~~~~

**Step 1: Prepare data and select parameters**

.. code-block:: python

   from edge_gwas import EDGEAnalysis
   from edge_gwas.utils import additive_gwas
   
   # For continuous outcome with right-skewed distribution
   edge = EDGEAnalysis(
       outcome_type='continuous',
       outcome_transform='log',  # or 'rank_inverse_normal'
       ols_method='bfgs',  # default, suitable for most cases
       max_iter=1000,
       verbose=True
   )

**Step 2: Run both EDGE and additive GWAS**

.. code-block:: python

   # EDGE analysis
   alpha_df, edge_results = edge.run_full_analysis(
       train_g, train_p, test_g, test_p, 
       outcome='triglycerides', 
       covariates=['age', 'sex', 'bmi', 'pc1', 'pc2', 'pc3'],
       output_prefix='results/edge_triglycerides'
   )
   
   # Additive GWAS for comparison
   additive_results = additive_gwas(
       test_g, test_p, 
       outcome='triglycerides',
       covariates=['age', 'sex', 'bmi', 'pc1', 'pc2', 'pc3'],
       outcome_type='continuous',
       outcome_transform='log'
   )

**Step 3: Compare results**

.. code-block:: python

   # Identify SNPs significant in EDGE but not additive
   edge_sig = edge_results[edge_results['pval'] < 5e-8]
   add_sig = additive_results[additive_results['pval'] < 5e-8]
   
   edge_only = edge_sig[~edge_sig['variant_id'].isin(add_sig['variant_id'])]
   print(f"SNPs significant only in EDGE: {len(edge_only)}")
   
   # Examine alpha values for significant SNPs
   print("\nAlpha distribution for EDGE-specific findings:")
   print(edge_only['alpha_value'].describe())

**Step 4: Validate findings**

* Replicate in independent cohort
* Check α consistency across cohorts
* Perform functional follow-up
* Check biological plausibility

**Step 5: Handle convergence issues**

.. code-block:: python

   # Check skipped SNPs
   skipped = edge.get_skipped_snps()
   print(f"Skipped {len(skipped)} SNPs")
   
   # If many SNPs fail, try different method
   if len(skipped) > 0.1 * len(test_g.columns):
       edge_robust = EDGEAnalysis(
           outcome_type='continuous',
           outcome_transform='log',
           ols_method='nm',  # More robust
           max_iter=5000
       )
       # Re-run analysis
       alpha_df2, edge_results2 = edge_robust.run_full_analysis(
           train_g, train_p, test_g, test_p,
           outcome='triglycerides',
           covariates=['age', 'sex', 'bmi', 'pc1', 'pc2', 'pc3']
       )

Computational Complexity
~~~~~~~~~~~~~~~~~~~~~~~~

**Time Complexity:**

* Training stage: :math:`O(n_{train} \times m \times p^2 \times k)`
* Test stage: :math:`O(n_{test} \times m \times p^2)`

Where:

* :math:`n` = sample size
* :math:`m` = number of variants
* :math:`p` = number of covariates
* :math:`k` = iterations for convergence (method-dependent)

**Space Complexity:**

* :math:`O(n \times m)` for genotype data
* Additional :math:`O(m)` for α values storage

**Optimization Strategies:**

* Process chromosomes separately
* Use parallel processing (multiple cores)
* Filter variants before analysis
* Use sparse matrix representations for rare variants
* Choose efficient optimization method (BFGS or L-BFGS)

**Method-specific performance:**

.. list-table::
   :header-rows: 1
   :widths: 25 25 50

   * - Method
     - Relative Speed
     - Notes
   * - BFGS
     - 1.0x (baseline)
     - Good default choice
   * - Newton
     - 0.8x (faster)
     - When close to optimum
   * - L-BFGS
     - 1.1x
     - Slightly slower but memory-efficient
   * - Nelder-Mead
     - 3-5x (slower)
     - Robust but slow
   * - CG
     - 1.5x
     - Good for sparse problems
   * - NCG
     - 1.3x
     - Efficient for large covariate sets
   * - Powell
     - 4-6x (slower)
     - Derivative-free
   * - Basin-hopping
     - 50-100x (much slower)
     - Global optimization, rarely needed

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
* Multiple optimization algorithms for robustness
* Employ convergence tolerance of 1e-6 (configurable)
* Maximum iterations: configurable (default: 1000)
* Automatic handling of invalid transformation values

Convergence Issues
~~~~~~~~~~~~~~~~~~

**SNPs may fail to converge due to:**

1. **Rare variants**: Too few observations for stable estimates
2. **Perfect separation**: In logistic regression when outcome completely separated by genotype
3. **Multicollinearity**: Highly correlated covariates
4. **Numerical instability**: Extreme coefficient values
5. **Poor starting values**: Especially with Newton method
6. **Transformation issues**: Invalid values after transformation

**Solutions:**

* Increase ``max_iter`` parameter
* Try different optimization method:
  
  .. code-block:: python
  
     edge = EDGEAnalysis(
         outcome_type='continuous',
         ols_method='nm',  # More robust
         max_iter=5000
     )

* Filter rare variants (MAF > 0.01)
* Check covariate correlations
* Use outcome transformation
* Use ``get_skipped_snps()`` to identify problematic variants

**Systematic convergence failures:**

.. code-block:: python

   # If >10% of SNPs fail
   skipped = edge.get_skipped_snps()
   
   if len(skipped) / total_snps > 0.1:
       # Possible data quality issues
       print("Warning: High failure rate. Check:")
       print("1. MAF distribution")
       print("2. Covariate multicollinearity")
       print("3. Outcome distribution")
       print("4. Sample size per genotype group")

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

Limitations
-----------

Current Limitations
~~~~~~~~~~~~~~~~~~~

1. **Data splitting required**: Reduces effective sample size
2. **Assumes independence**: Not designed for related individuals (mixed models coming)
3. **Single SNP analysis**: Does not model joint effects
4. **Autosomal variants only**: X/Y chromosomes require special handling
5. **European-ancestry optimized**: May need adjustment for diverse populations
6. **Linear transformations**: Non-linear relationships not fully captured
7. **Computational cost**: Higher than standard GWAS (but manageable)

Statistical Assumptions
~~~~~~~~~~~~~~~~~~~~~~~

EDGE assumes:

1. **Independent samples**: No cryptic relatedness (unless using mixed models)
2. **No genotyping errors**: QC filters applied
3. **Hardy-Weinberg Equilibrium**: In controls for binary outcomes
4. **Correct model specification**: Covariates adequately control confounding
5. **Large sample size**: Asymptotic properties of tests hold
6. **Proper transformation**: Normality after transformation (for linear regression)
7. **Convergence**: Optimization reaches true optimum

Validation Studies
------------------

Simulation Studies
~~~~~~~~~~~~~~~~~~

EDGE has been validated through simulations showing:

* **Type I error control**: Maintains nominal false positive rate (5%)
* **Increased power**: 10-50% power gain for nonadditive effects
* **Robustness**: Stable performance across MAF ranges
* **α recovery**: Accurately estimates true inheritance pattern
* **Transformation effectiveness**: Proper error control with non-normal outcomes
* **Optimization stability**: Consistent results across optimization methods

Real Data Applications
~~~~~~~~~~~~~~~~~~~~~~

EDGE has been applied to:

* **Cardiometabolic traits**: BMI, diabetes, blood pressure, lipids
* **UK Biobank data**: >500,000 individuals
* **Novel discoveries**: Identified nonadditive effects missed by additive GWAS
* **Diverse outcomes**: Binary and continuous with various transformations

See **Zhou et al. (2023)** for detailed results.

Method Comparison Studies
~~~~~~~~~~~~~~~~~~~~~~~~~

Benchmarking shows:

* **BFGS vs Newton**: BFGS more robust, Newton slightly faster when converging
* **Transformation impact**: RINT provides best Type I error control for heavy-tailed traits
* **Power curves**: EDGE superior for α < 0.3 or α > 0.7
* **Replication rates**: Higher for EDGE-discovered loci with extreme α values

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
5. Transformation applied identically to both datasets maintains independence
6. Therefore, Type I error rate is controlled at nominal level

Consistency of α (Sketch)
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Theorem**: :math:`\hat{\alpha} \xrightarrow{p} \alpha_{true}` as :math:`n \rightarrow \infty`

**Conditions**:

* :math:`|\beta_{HA}| > c > 0` for some constant c
* Standard regularity conditions for MLE
* Convergence of optimization algorithm to true optimum

**Implication**: With sufficient training data and proper optimization, α estimates converge to true value.

Convergence of Optimization Methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Theorem**: Under standard regularity conditions, quasi-Newton methods (BFGS, L-BFGS) achieve superlinear convergence.

**Rate**: :math:`\|\theta_{k+1} - \theta^*\| \leq C \|\theta_k - \theta^*\|^{1+\epsilon}`

Where:

* :math:`\theta^*` = true parameter value
* :math:`\epsilon > 0` for superlinear convergence
* :math:`C` = positive constant

Glossary of Terms
-----------------

.. glossary::

   Alpha (α)
      Encoding parameter representing the ratio of heterozygous to homozygous effects
   
   Additive model
      Genetic model assuming heterozygote effect is exactly half of homozygous effect (α = 0.5)
   
   BFGS
      Broyden-Fletcher-Goldfarb-Shannon algorithm for optimization (default for EDGE)
   
   Codominant model
      Model that separately estimates effects for heterozygotes and homozygotes
   
   EDGE encoding
      Genotype encoding using estimated α values
   
   Genomic inflation (λ)
      Measure of systematic p-value inflation due to population structure
   
   HWE
      Hardy-Weinberg Equilibrium
   
   INT
      Inverse Normal Transformation (parametric)
   
   L-BFGS
      Limited-memory BFGS, memory-efficient variant for large datasets
   
   MAF
      Minor Allele Frequency
   
   Nelder-Mead
      Derivative-free simplex optimization method
   
   Newton-Raphson
      Second-order optimization using exact Hessian matrix
   
   Over-dominant
      Inheritance pattern where heterozygotes have stronger effect than either homozygote (α > 1)
   
   Quasi-Newton methods
      Optimization methods using approximation of Hessian (e.g., BFGS, L-BFGS)
   
   Recessive
      Inheritance pattern where heterozygotes behave like reference homozygotes (α ≈ 0)
   
   RINT
      Rank-based Inverse Normal Transformation, robust to outliers
   
   Wald test
      Statistical test using ratio of estimate to standard error
   
   Two-stage analysis
      Analysis approach separating parameter estimation (stage 1) from hypothesis testing (stage 2)

Best Practices Summary
----------------------

Data Preparation
~~~~~~~~~~~~~~~~

1. **Quality Control:**
   
   * MAF > 0.01 (or 0.05 for small samples)
   * HWE p > 1×10⁻⁶ in controls
   * Genotype call rate > 95%
   * Sample call rate > 95%
   * Remove related individuals (kinship > 0.0884)

2. **Covariates:**
   
   * Include age, sex, assessment center (if applicable)
   * Add 10-20 principal components
   * Check for multicollinearity (VIF < 10)
   * Consider batch effects

3. **Outcome Preparation:**
   
   * Check distribution (histogram, Q-Q plot)
   * Remove extreme outliers (>5 SD) before transformation
   * Choose appropriate transformation:
     
     * Normal → None
     * Right-skewed → ``'log'`` or ``'rank_inverse_normal'``
     * Heavy-tailed → ``'rank_inverse_normal'``
     * Unknown → ``'rank_inverse_normal'`` (safest)
   
   * Verify transformation effectiveness

Parameter Selection
~~~~~~~~~~~~~~~~~~~

1. **Optimization Method:**
   
   * Default: ``ols_method='bfgs'`` (suitable for >95% of cases)
   * Large datasets (N>100k): ``ols_method='lbfgs'``
   * Convergence issues: ``ols_method='nm'`` or ``ols_method='newton'``
   * Many covariates (P>100): ``ols_method='ncg'``

2. **Maximum Iterations:**
   
   * Default: ``max_iter=1000`` (usually sufficient)
   * Increase to 2000-5000 if convergence issues
   * Monitor skipped SNPs

3. **Data Splitting:**
   
   * 50/50 split recommended (balanced power and stability)
   * Larger training set (60/40) for unstable traits
   * Larger test set (40/60) for well-behaved traits


Result Interpretation
~~~~~~~~~~~~~~~~~~~~~

1. **Genome-wide significant hits (p < 5×10⁻⁸):**
   
   * Report α value for each hit
   * Classify inheritance pattern
   * Compare with additive GWAS results
   * Check for biological plausibility

2. **Alpha interpretation:**
   
   .. code-block:: python
   
      def classify_inheritance(alpha):
          if alpha < 0:
              return "Under-recessive (heterozygote disadvantage)"
          elif alpha < 0.3:
              return "Recessive to partially recessive"
          elif 0.45 <= alpha <= 0.55:
              return "Additive"
          elif 0.7 < alpha < 1.3:
              return "Partially dominant to dominant"
          elif alpha >= 1.3:
              return "Over-dominant (heterozygote advantage)"
          else:
              return "Intermediate"
      
      sig_hits = gwas_results[gwas_results['pval'] < 5e-8]
      sig_hits['inheritance'] = sig_hits['alpha_value'].apply(classify_inheritance)
      print(sig_hits[['variant_id', 'alpha_value', 'inheritance', 'pval']])

3. **Quality checks:**
   
   * λ should be ≈1.0 (acceptable: 0.95-1.05)
   * Q-Q plot should show no systematic deviation
   * Check consistency of α values for nearby SNPs (LD)

Appendix A: Optimization Method Details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**BFGS Algorithm:**

The Broyden-Fletcher-Goldfarb-Shannon algorithm approximates the Hessian matrix using gradient information:

.. math::

   H_{k+1} = H_k + \frac{y_k y_k^T}{y_k^T s_k} - \frac{H_k s_k s_k^T H_k}{s_k^T H_k s_k}

where:

- :math:`s_k = x_{k+1} - x_k` (step)
- :math:`y_k = \nabla f(x_{k+1}) - \nabla f(x_k)` (gradient change)
- :math:`H_k` is the approximate Hessian at iteration k

**Convergence criteria:**

.. math::

   \|\nabla f(x_k)\| < \epsilon \quad \text{or} \quad |f(x_k) - f(x_{k-1})| < \epsilon

Default tolerance: :math:`\epsilon = 10^{-6}`

**L-BFGS Algorithm:**

Limited-memory BFGS stores only m recent vector pairs (typically m=10):

Memory requirement: :math:`O(m \cdot n)` instead of :math:`O(n^2)` for full BFGS

Trade-off: Slightly slower convergence but much lower memory usage

**Newton-Raphson:**

Uses exact Hessian matrix:

.. math::

   x_{k+1} = x_k - H^{-1}(x_k) \nabla f(x_k)

where :math:`H(x_k)` is the true Hessian matrix.

**Convergence rate:** Quadratic (fastest when close to optimum)

**Cost:** Expensive Hessian computation at each iteration

**Nelder-Mead Simplex:**

Derivative-free method using simplex of n+1 points:

Operations: reflection, expansion, contraction, shrinkage

**Pros:** Robust, no gradient needed
**Cons:** Slow, may not converge for high dimensions

See Also
--------

**Documentation:**

* :ref:`installation` - Installation instructions and requirements
* :ref:`quickstart` - Getting started guide with simple examples
* :ref:`user_guide` - Comprehensive user guide and tutorials
* :ref:`api_reference` - Complete API documentation
* :ref:`examples` - Example analyses and case studies
* :ref:`visualization` - Plotting and visualization guide
* :ref:`faq` - Frequently asked questions

---

*Last updated: 2025-12-25 for edge-gwas v0.1.1*

*For questions or issues, visit:* https://github.com/nicenzhou/edge-gwas/issues
