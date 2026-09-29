from setuptools import setup, find_packages

setup(
    name="drug-discovery-pipeline",
    version="0.1.0",
    author="AI Drug Discovery Contributors",
    description="AI-driven drug discovery pipeline: from computational genetics to protein design",
    long_description=open("README.md", encoding="utf-8").read() if open("README.md", encoding="utf-8") else "",
    long_description_content_type="text/markdown",
    url="https://github.com/alok-108/drug-discovery-pipeline",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
    ],
    python_requires=">=3.10",
    install_requires=[
        "numpy>=1.26.0",
        "pandas>=2.1.0",
        "scipy>=1.11.0",
        "pyyaml>=6.0.1",
        "requests>=2.31.0",
        "networkx>=3.1",
        "biopython>=1.83",
        "matplotlib>=3.8.0",
    ],
    extras_require={
        "dev": [
            "pytest>=7.4.0",
            "pytest-cov>=4.1.0",
            "ruff>=0.1.0",
        ]
    },
    entry_points={
        "console_scripts": [
            "drug-pipeline=utils.io_utils:main_cli",
        ],
    },
)
