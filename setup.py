from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="edge-gwas",
    version="0.1.0",
    author="Jiayan Zhou",
    author_email="jyzhou@stanford.edu",
    description="EDGE (Encoding Deviation Genotypic Effects) GWAS: Flexibly encoded GWAS for identifying nonadditive SNPs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/edge-gwas",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    python_requires=">=3.7",
    install_requires=requirements,
    extras_require={
        "dev": ["pytest>=6.0", "pytest-cov>=2.0"],
    },
)
