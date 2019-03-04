"""A setuptools based setup module.
See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path
import LDpred

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='LDpred',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=LDpred.__version__,

    description='A Python package that adjusts GWAS summary statistics for the effects of linkage disequilibrium (LD)',
    long_description=long_description,

    # The project's main homepage.
    url='https://github.com/ldpred',

    # Author details
    author='Bjarni J Vilhjalmsson',
    author_email='bjarni.vilhjalmsson@gmail.com',

    # Choose your license
    license='MIT',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 4 - Beta',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Software Development :: Build Tools',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: MIT License',

        #'Programming Language :: Python :: 2',
        #'Programming Language :: Python :: 2.6',
        #'Programming Language :: Python :: 2.7',
        #'Programming Language :: Python :: 3',
        #'Programming Language :: Python :: 3.3',
        #'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],

    # What does your project relate to?
    keywords='Polygenic Risk Scores, GWAS, Linkage Disequilibrium, Risk Prediction',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(),

    # Alternatively, if you want to distribute just a my_module.py, uncomment
    # this:
    #   py_modules=["my_module"],

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=['h5py','scipy','plinkio'],

    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # $ pip install -e .[dev,test]
    extras_require={
        #'dev': ['check-manifest'],
        #'test': ['coverage'],
    },

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    package_data={'test_data':
                ['./test_data/LDpred_data_p1.000_test_0.fam','./test_data/LDpred_data_p1.000_test_0.bim','./test_data/LDpred_data_p1.000_test_0.bed',
                './test_data/LDpred_data_p1.000_ss_0.txt','./test_data/LDpred_data_p1.000_train_0.fam','./test_data/LDpred_data_p1.000_train_0.bim',
                './test_data/LDpred_data_p1.000_train_0.bed','./test_data/LDpred_data_p0.001_test_0.fam','./test_data/LDpred_data_p0.001_test_0.bim',
                './test_data/LDpred_data_p0.001_test_0.bed','./test_data/LDpred_data_p0.001_ss_0.txt','./test_data/LDpred_data_p0.001_train_0.fam',
                './test_data/LDpred_data_p0.001_train_0.bim','./test_data/LDpred_data_p0.001_train_0.bed','./test_data/LDpred_cc_data_p0.001_test_0.fam',
                './test_data/LDpred_cc_data_p0.001_test_0.bim','./test_data/LDpred_cc_data_p0.001_test_0.bed','./test_data/LDpred_cc_data_p0.001_ss_0.txt',
                './test_data/LDpred_cc_data_p0.001_train_0.fam','./test_data/LDpred_cc_data_p0.001_train_0.bim','./test_data/LDpred_cc_data_p0.001_train_0.bed',
                './reference/hm3_sids.txt.gz']
                },
    include_package_data=True,
    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'

#     scripts=['ldpred/LDpred.py', 'ldpred/LDpred_inf.py', 'ldpred/LD_pruning_thres.py', 
#              'ldpred/validate.py', 'ldpred/coord_genotypes.py'],
    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        'console_scripts': [
            'ldpred=LDpred:main',
        ],
    },
)
