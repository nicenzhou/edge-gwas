from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Read all dependencies from single requirements.txt
with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#") and not line.startswith("=")]

setup(
    name="edge-gwas",
    version="0.1.1",
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
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
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
