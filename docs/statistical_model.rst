.. _statistical_model:

Statistical Model
=================

Overview
--------

EDGE (Elastic Data-Driven Encoding) detects nonadditive genetic effects through a two-stage framework:

1. **Training**: Calculate encoding parameters (α)
2. **Testing**: Apply α for association testing

Regression Models
-----------------

**Stage 1: Codominant Model (Training)**

.. math::

   E(Y) = \beta_0 + \beta_{Het} \cdot SNP_{Het} + \beta_{HA} \cdot SNP_{HA} + \sum \beta_i \cdot COV_i

**Stage 2: EDGE-Encoded Model (Testing)**

.. math::

   E(Y) = \beta_0 + \beta \cdot G_{EDGE} + \sum \beta_i \cdot COV_i

Where:

.. math::

   G_{EDGE} = \begin{cases}
   1 & \text{if 0/0} \\
   \alpha & \text{if 0/1} \\
   0 & \text{if 1/1}
   \end{cases}

Encoding Parameter
------------------

.. math::

   \alpha = \frac{\beta_{Het}}{\beta_{HA}}

**Interpretation:**

* α < 0: Under-recessive
* α ≈ 0: Recessive
* α ≈ 0.5: Additive
* α ≈ 1: Dominant
* α > 1: Over-dominant

Outcome Types
-------------

**Binary (Logistic Regression):**

.. math::

   logit(P(Y=1)) = \beta_0 + \beta \cdot G_{EDGE} + \sum \beta_i \cdot COV_i

**Continuous (Linear Regression):**

.. math::

   Y = \beta_0 + \beta \cdot G_{EDGE} + \sum \beta_i \cdot COV_i + \epsilon

Optimization Methods
--------------------

Available algorithms for continuous outcomes:

* ``'bfgs'`` (default): Quasi-Newton, general purpose
* ``'lbfgs'``: Memory-efficient for large datasets
* ``'newton'``: Fast when close to optimum
* ``'nm'``: Robust, derivative-free
* ``'cg'``: Efficient for sparse problems

**Usage:**

.. code-block:: python

   edge = EDGEAnalysis(
       outcome_type='continuous',
       ols_method='bfgs'  # or 'lbfgs', 'newton', etc.
   )

Outcome Transformations
-----------------------

For continuous outcomes:

* ``'log'``: Natural log for right-skewed data (Y > 0)
* ``'log10'``: Base-10 log
* ``'inverse_normal'``: Parametric INT
* ``'rank_inverse_normal'``: Rank-based INT (robust)

**Usage:**

.. code-block:: python

   edge = EDGEAnalysis(
       outcome_type='continuous',
       outcome_transform='rank_inverse_normal'
   )

Statistical Testing
-------------------

**Hypothesis:**

* H₀: β = 0
* Hₐ: β ≠ 0

**Test statistic:**

.. math::

   Z = \frac{\hat{\beta}}{SE(\hat{\beta})}

**Significance thresholds:**

* Genome-wide: p < 5×10⁻⁸
* Suggestive: p < 1×10⁻⁵

Quality Control
---------------

**Recommendations:**

* MAF > 0.01
* HWE p > 1×10⁻⁶
* Missing rate < 5%
* Include 10+ PCs as covariates

**Minimum sample sizes:**

* Binary: 500+ cases and controls
* Continuous: 1,000+ samples

Best Practices
--------------

1. **Check outcome distribution, apply transformation if needed**
2. **Use default BFGS optimization (try lbfgs for N>100k)**
3. **Include adequate covariates (age, sex, 10+ PCs)**
4. **Compare with additive GWAS**
5. **Verify α consistency for significant hits**


See Also
--------

**Documentation:**

* `Documentation Home <index.html>`_ - Home
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
