from setuptools import setup

setup(
    name='PhyMMR',
    version='1.0.0',
    packages=['phymmr', 'phymmr.rocky'],
    url='https://github.com/Thernn88/PhyMMR',
    license='GPLv3+',
    author='Ashton Pooley',
    author_email='',
    description='',
    install_requires=[
        "biopython>=1.79",
        "numpy>=1.23.3",
        "tqdm>=4.64.1",
        "wrap_rocks>=0.3.1",
        "xxhash>=3.0.0",
        "phymmr_tools>=0.2.4",
        "pro2codon>=1.2.4",
    ],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        "Environment :: Console",
        "Operating System :: OS Independent",
    ],
    entry_points={
        'console_scripts': [
            'phymmr = phymmr.__main__:main',
        ]
    }
)
