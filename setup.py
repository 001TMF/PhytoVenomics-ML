#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Phytovenomics

A machine learning platform for plant-based antivenom development.
This setup script installs the package and its dependencies.
"""

import os
from setuptools import setup, find_packages

# Read the long description from README.md
with open('../../Downloads/README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

# Package requirements
requirements = [
    'numpy>=1.20.0',
    'pandas>=1.3.0',
    'scikit-learn>=1.0.0',
    'matplotlib>=3.4.0',
    'seaborn>=0.11.0',
    'tensorflow>=2.8.0',
    'torch>=1.10.0',
    'pyyaml>=6.0',
    'click>=8.0.0',
    'tqdm>=4.62.0',
    'requests>=2.26.0',
    'scipy>=1.7.0',
    'biopython>=1.79',
    'rdkit>=2022.03.1',
    'joblib>=1.1.0',
]

# Development requirements
dev_requirements = [
    'pytest>=7.0.0',
    'pytest-cov>=3.0.0',
    'flake8>=4.0.0',
    'black>=22.1.0',
    'sphinx>=4.4.0',
    'sphinx-rtd-theme>=1.0.0',
    'pre-commit>=2.17.0',
]

# Testing requirements
test_requirements = [
    'pytest>=7.0.0',
    'pytest-cov>=3.0.0',
]

# Data packages (optional)
data_requirements = [
    'h5py>=3.6.0', 
    'pillow>=9.0.0',
    'pyarrow>=7.0.0',
]

# Version extraction
VERSION = {}
with open(os.path.join('phytovenomics', '_version.py'), 'r') as f:
    exec(f.read(), VERSION)

setup(
    name='phytovenomics',
    version=VERSION.get('__version__', '0.1.0'),
    description='Machine learning platform for plant-based antivenom development',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Phytovenomics Team',
    author_email='info@phytovenomics.org',
    url='https://github.com/phytovenomics/phytovenomics',
    packages=find_packages(include=['phytovenomics', 'phytovenomics.*']),
    package_data={
        'phytovenomics': [
            'config/*.yaml',
            'data/reference/*.csv',
            'data/schemas/*.json',
        ],
    },
    python_requires='>=3.8',
    install_requires=requirements,
    extras_require={
        'dev': dev_requirements,
        'test': test_requirements,
        'data': data_requirements,
        'all': dev_requirements + test_requirements + data_requirements,
    },
    entry_points={
        'console_scripts': [
            'phytovenomics-setup=phytovenomics.cli.setup_cli:main',
            'phytovenomics-run=phytovenomics.cli.run_cli:main',
            'phytovenomics-analyze=phytovenomics.cli.analyze_cli:main',
        ],
    },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
    ],
    license='MIT',
    zip_safe=False,
)