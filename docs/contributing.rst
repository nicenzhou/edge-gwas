.. _contributing:

Contributing Guide
==================

Thank you for your interest in contributing to edge-gwas!

This guide will help you get started with contributing to the project.

Ways to Contribute
------------------

You can contribute to edge-gwas in several ways:

* **Report bugs** - Submit bug reports via GitHub Issues
* **Suggest features** - Propose new features or improvements
* **Improve documentation** - Fix typos, add examples, clarify explanations
* **Submit code** - Fix bugs or implement new features
* **Share use cases** - Tell us how you're using edge-gwas
* **Help others** - Answer questions in GitHub Discussions

Getting Started
---------------

1. Fork the Repository
~~~~~~~~~~~~~~~~~~~~~~

Visit https://github.com/nicenzhou/Visit https://github.com/ni

2. Clone Your Fork
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   git clone https://github.com/nicenzhou/edge-gwas.git
   cd edge-gwas

3. Set Up Development Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Create virtual environment
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   
   # Install development dependencies
   pip install -r requirements.txt
   pip install -e ".[dev]"

4. Create a Branch
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   git checkout -b feature/my-new-feature
   # or
   git checkout -b fix/bug-description

Development Workflow
--------------------

Code Style
~~~~~~~~~~

* Follow PEP 8 style guide
* Use descriptive variable names
* Add docstrings to all functions and classes
* Keep functions focused and concise

**Example:**

.. code-block:: python

   def calculate_maf(genotypes):
       """
       Calculate minor allele frequency for each variant.
       
       Parameters
       ----------
       genotypes : pandas.DataFrame
           Genotype matrix (samples × variants)
       
       Returns
       -------
       pandas.Series
           MAF for each variant
       
       Examples
       --------
       >>> maf = calculate_maf(genotype_df)
       >>> print(maf.head())
       """
       # Implementation here
       pass

Documentation
~~~~~~~~~~~~~

* Add docstrings using NumPy style
* Update relevant .rst files in docs/
* Add examples for new features
* Update changelog

**Docstring template:**

.. code-block:: python

   def my_function(param1, param2, param3=None):
       """
       Brief description of function.
       
       Longer description with more details if needed.
       
       Parameters
       ----------
       param1 : type
           Description of param1
       param2 : type
           Description of param2
       param3 : type, optional
           Description of param3 (default: None)
       
       Returns
       -------
       return_type
           Description of return value
       
       Raises
       ------
       ExceptionType
           When this exception is raised
       
       Examples
       --------
       >>> result = my_function(arg1, arg2)
       >>> print(result)
       """
       pass

Testing
~~~~~~~

* Write unit tests for new functions
* Ensure all tests pass before submitting
* Aim for high code coverage

.. code-block:: bash

   # Run tests
   pytest tests/
   
   # Run with coverage
   pytest --cov=edge_gwas tests/

**Example test:**

.. code-block:: python

   import pytest
   from edge_gwas import EDGEAnalysis
   
   def test_edge_analysis_initialization():
       """Test EDGEAnalysis initialization"""
       edge = EDGEAnalysis(outcome_type='binary')
       assert edge.outcome_type == 'binary'
   
   def test_invalid_outcome_type():
       """Test that invalid outcome type raises error"""
       with pytest.raises(ValueError):
           EDGEAnalysis(outcome_type='invalid')

Commit Messages
~~~~~~~~~~~~~~~

Use conventional commit format:

.. code-block:: text

   type(scope): brief description
   
   Longer description if needed.
   
   Fixes #issue_number

**Types:**

* ``feat``: New feature
* ``fix``: Bug fix
* ``docs``: Documentation changes
* ``style``: Code style changes (formatting)
* ``refactor``: Code refactoring
* ``test``: Adding or updating tests
* ``chore``: Maintenance tasks

**Examples:**

.. code-block:: text

   feat(gwas): add support for VCF file format
   
   fix(viz): correct lambda calculation in QQ plot
   
   docs(api): add examples for calculate_alpha function
   
   test(utils): add tests for data loading functions

Submitting Changes
------------------

1. Push to Your Fork
~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   git add .
   git commit -m "feat(module): description of changes"
   git push origin feature/my-new-feature

2. Create Pull Request
~~~~~~~~~~~~~~~~~~~~~~

1. Go to https://github.com/nicenzhou/edge-gwas
2. Click "Pull requests" → "New pull request"
3. Select your branch
4. Fill out the PR template:

   * **Title**: Brief description
   * **Description**: Detailed explanation of changes
   * **Related Issues**: Link to related issues
   * **Testing**: Describe how you tested changes
   * **Checklist**: Complete the checklist

**PR Template:**

.. code-block:: markdown

   ## Description
   Brief description of changes
   
   ## Related Issues
   Fixes #123
   
   ## Changes Made
   - Added feature X
   - Fixed bug in Y
   - Updated documentation for Z
   
   ## Testing
   - [ ] All tests pass
   - [ ] Added new tests
   - [ ] Manual testing performed
   
   ## Checklist
   - [ ] Code follows style guidelines
   - [ ] Documentation updated
   - [ ] Tests added/updated
   - [ ] Changelog updated

3. Code Review
~~~~~~~~~~~~~~

* Respond to reviewer comments
* Make requested changes
* Update your branch if needed

.. code-block:: bash

   # Update your branch with latest main
   git checkout main
   git pull upstream main
   git checkout feature/my-new-feature
   git rebase main

Reporting Bugs
--------------

When reporting bugs, please include:

1. **Description**: Clear description of the bug
2. **Steps to Reproduce**: Minimal example to reproduce
3. **Expected Behavior**: What you expected to happen
4. **Actual Behavior**: What actually happened
5. **Environment**: Python version, OS, package versions
6. **Error Messages**: Full error traceback

**Bug Report Template:**

.. code-block:: markdown

   **Describe the bug**
   A clear description of what the bug is.
   
   **To Reproduce**
   ```python
   # Minimal code to reproduce
   from edge_gwas import EDGEAnalysis
   edge = EDGEAnalysis()
   # ... steps to reproduce
