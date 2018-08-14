import os
import shutil
import sys
import glob
from setuptools import setup, find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

version = '0.0.1'
if os.path.exists('VERSION'):
  version = open('VERSION').read().strip()

setup(
    name='gffmunger',
    version=version,
    description='Munger of GFF files.  Proper description to follow.',
	 long_description=read('README.md'),
    packages = find_packages(),
    author='Tim Stickland',
    author_email='path-help@sanger.ac.uk',
    url='https://github.com/sanger-pathogens/gffmunger',
    scripts=glob.glob('scripts/*'),
    data_files=[('config', ['gffmunger-config.yml'])],
    test_suite='nose.collector',
    tests_require=['nose >= 1.3'],
    install_requires=[
         'biopython >= 1.68',
         #'pyfastaq >= 3.12.0'
         'gffutils'
       ],
    license='GPLv3',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience  :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)'
    ],
)
