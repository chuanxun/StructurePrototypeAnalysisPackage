from setuptools import setup, find_packages

setup(
    name='spap',
    version='1.0.0',
    packages=find_packages(),
    # package_dir={'':'spap'},
    url='https://github.com/chuanxun/StructurePrototypeAnalysisPackage',
    license='GPL',
    author='Chuanxun Su',
    author_email='suchuanxun@163.com',
    keywords=["crystal structure", 'prototype', 'spap', 'similarity comparison', "clustering", 'unsupervised learning',
              'pattern recognition', 'machine learning', 'AI', 'structure prediction', 'DFT', 'CALYPSO'],
    platforms="Independant",
    # install_requires=['numpy>=1.8.0', 'spglib>=1.10.0', 'ase>=3.13.0'],
    install_requires=['numpy', 'spglib', 'ase'],
    description='This program can analyze symmetry and compare similarity of atomic structures.',
    long_description='Structure Prototype Analysis Package (SPAP) can analyze symmetry and compare similarity of a '
                     'large number of atomic structures. Typically, SPAP can process structures predicted by CALYPSO '
                     '(www.calypso.cn).',
    classifiers=[
        "Programming Language :: Python",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
    ],
    entry_points={
        'console_scripts': [
            'spap = spap.spap:start_cli',
        ]
    }
)
