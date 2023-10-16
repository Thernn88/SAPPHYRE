from setuptools import setup

setup(name='sapphyre',
      version='0.1',
      description='An assembly-less solution for processing high-throughput sequencing reads for phylogenetics',
      url='https://github.com/Thernn88/SAPPHYRE',
      author='Sapphyre TEAM',
      license='GPL3',
      py_modules=['sapphyre'],
      zip_safe=False,
      entry_points={
        'console_scripts': [
            'sapphyre = sapphyre.__main__:main',
        ],
    },
    install_requires=[
        "phymmr_tools>=0.5.10",
        "wrap_rocks>=0.3.7",
        "biopython>=1.79",
        "numpy>=1.23.3",
        "blosum==2.0.2",
        "pro2codon>=1.2.4",
        "pyfastx>=2.0.1",
        "needletail>=0.5.0",
        "pandas>=2.1.1",
        "msgspec>=0.18.2",
        "parasail>=1.3.4",
        "xxhash>=3.3.0",
        "pyarrow>=13.0.0",
        "openpyxl>=3.1.2",
        "requests>=2.31.0",
        "beautifulsoup4>=4.12.2",
        "tqdm>=4.66.1"
    ]
)