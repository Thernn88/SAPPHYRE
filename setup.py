from setuptools import setup, find_packages

setup(
    name="sapphyre",
    version="0.1.3",
    description="An assembly-less solution for processing high-throughput sequencing reads for phylogenetics",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/Thernn88/SAPPHYRE",
    author="Kevin Moran",
    author_email="thernn8891@gmail.com",
    maintainer="Ashton Pooley",
    maintainer_email="admin@ashi.digital",
    license="GPLv3",
    packages=find_packages(where="."),
    zip_safe=False,
    entry_points={
        "console_scripts": [
            "sapphyre = sapphyre.__main__:main",
        ],
    },
    install_requires=[
        "beautifulsoup4>=4.12.2",
        "Bio>=1.6.0",
        "biopython>=1.81",
        "msgspec>=0.18.4",
        "needletail>=0.5.0",
        "numpy>=1.24.2",
        "openpyxl>=3.1.2",
        "pandas>=2.1.3",
        "parasail>=1.3.4",
        "phymmr_tools>=0.6.5",
        "pro2codon>=1.2.4",
        "pyfastx>=2.0.1",
        "Requests>=2.31.0",
        "scipy>=1.11.4",
        "setuptools>=67.6.1",
        "tqdm>=4.66.1",
        "wrap_rocks>=0.3.7",
        "xxhash>=3.4.1",
        "pr2codon>=1.1.14",
    ],
    dependency_links=[
        'https://pypi.org/simple/',
    ],
    keywords=[
        "phylogenetics",
        "bioinformatics",
        "genomics",
        "phylogenomics",
        "phylogeny",
        "bio",
        "assembly",
        "sequencing",
    ],
    platforms="any",
    classifiers=[
        "Natural Language :: English",
        "Environment :: Console",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ]
)
