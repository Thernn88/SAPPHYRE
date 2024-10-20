from setuptools import find_packages, setup

setup(
    name="sapphyre",
    version="0.3.8",
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
        "beautifulsoup4>=4.12.3",
        "biopython>=1.83",
        "blosum>=2.0.3",
        "msgspec>=0.18.6",
        "needletail>=0.5.0",
        "numpy>=2.0.1",
        "openpyxl>=3.1.2",
        "pandas>=2.2.2",
        "parasail>=1.3.4",
        "pr2codon>=1.1.18",
        "psutil>=5.9.8",
        "pyarrow>=17.0.0",
        "pyfamsa>=0.4.0",
        "requests>=2.32.3",
        "sapphyre_tools>=0.9.0",
        "setuptools>=69.5.1",
        "tqdm>=4.66.2",
        "wrap_rocks>=0.3.8",
        "xxhash>=3.4.1",
        "isal>=1.7.0",
        
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
