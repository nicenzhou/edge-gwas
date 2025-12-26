from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Core dependencies (required)
with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="edge-gwas",
    version="0.1.1",  # Updated version
    author="Jiayan Zhou", 
    author_email="jyzhou@stanford.edu",
    description="EDGE (Elastic Data-Driven Encoding) GWAS: Flexibly encoded GWAS for identifying nonadditive SNP effects",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/nicenzhou/edge-gwas", 
    project_urls={
        "Documentation": "https://edge-gwas.readthedocs.io",
        "Source": "https://github.com/nicenzhou/edge-gwas",
        "Bug Tracker": "https://github.com/nicenzhou/edge-gwas/issues",
    },
    packages=find_packages(exclude=["tests", "tests.*", "docs", "docs.*"]),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",  # Changed from MIT
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    install_requires=requirements,
    extras_require={
        # File format support
        "pgen": ["pgenlib>=0.81.0"],
        "bgen": ["bgen-reader>=4.0.8"],
        "vcf": ["cyvcf2>=0.30.0"],
        "all": [
            "pgenlib>=0.81.0",
            "bgen-reader>=4.0.8",
            "cyvcf2>=0.30.0",
        ],
        # Development tools
        "dev": [
            "pytest>=6.2.0",
            "pytest-cov>=2.12.0",
            "black>=21.0",
            "flake8>=3.9.0",
            "isort>=5.9.0",
            "mypy>=0.910",
        ],
        # Documentation
        "docs": [
            "sphinx>=4.0.0",
            "sphinx-rtd-theme>=1.0.0",
            "myst-parser>=0.15.0",
        ],
    },
    keywords=[
        "gwas", 
        "genomics", 
        "genetics", 
        "bioinformatics", 
        "nonadditive", 
        "flexible-encoding",
        "edge",
        "association-studies",
    ],
)
