"""
Setup script for Antibody Design Pipeline
"""

from setuptools import setup, find_packages
import os

# Read README
def read_readme():
    readme_path = os.path.join(os.path.dirname(__file__), 'README_PIPELINE.md')
    if os.path.exists(readme_path):
        with open(readme_path, 'r', encoding='utf-8') as f:
            return f.read()
    return "Integrated Antibody Design Pipeline"

# Read version
def get_version():
    version_file = os.path.join(os.path.dirname(__file__), 'antibody_pipeline', '__init__.py')
    with open(version_file, 'r') as f:
        for line in f:
            if line.startswith('__version__'):
                return line.split('"')[1]
    return "1.0.0"

setup(
    name='antibody-design-pipeline',
    version=get_version(),
    description='Integrated antibody design pipeline combining DiffAb, Germinal, and RFdiffusion',
    long_description=read_readme(),
    long_description_content_type='text/markdown',
    author='Phytovenomics Team',
    author_email='',
    url='https://github.com/your-repo/antibody-pipeline',
    packages=find_packages(),
    python_requires='>=3.8',
    install_requires=[
        # Core dependencies
        'numpy>=1.21.0',
        'torch>=2.0.0',
        'biopython>=1.79',

        # Structure prediction
        'abnumber>=0.3.0',  # Chothia numbering

        # IgLM
        'iglm>=0.4.0',
        'transformers>=4.30.0',

        # Configuration
        'pyyaml>=6.0',

        # Optional but recommended
        'scipy>=1.7.0',
        'pandas>=1.3.0',
    ],
    extras_require={
        'dev': [
            'pytest>=7.0.0',
            'pytest-cov>=3.0.0',
            'black>=22.0.0',
            'flake8>=4.0.0',
            'mypy>=0.950',
        ],
        'chai': [
            'chai_lab>=0.1.0',  # Chai-1 structure prediction
        ],
        'pyrosetta': [
            # PyRosetta requires special installation
            # See: https://www.pyrosetta.org/downloads
        ],
    },
    entry_points={
        'console_scripts': [
            'antibody-design=antibody_pipeline.cli:main',
        ],
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    keywords='antibody design machine-learning structure-prediction',
    project_urls={
        'Documentation': 'https://github.com/your-repo/antibody-pipeline/docs',
        'Source': 'https://github.com/your-repo/antibody-pipeline',
        'Bug Reports': 'https://github.com/your-repo/antibody-pipeline/issues',
    },
)