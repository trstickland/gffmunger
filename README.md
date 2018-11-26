# GFF munger
Munges GFF3 files exported from Chado

[![Build Status](https://travis-ci.org/sanger-pathogens/gffmunger.svg?branch=master)](https://travis-ci.org/sanger-pathogens/gffmunger)   
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/sanger-pathogens/gffmunger/blob/master/LICENSE)   
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/gffmunger/README.html)   
[![Container ready](https://img.shields.io/badge/container-ready-brightgreen.svg)](https://quay.io/repository/biocontainers/gffmunger)  

## Contents
  * [Introduction](#introduction)
  * [Installation](#installation)
    * [Using pip3](#using-pip3)
    * [Conda](#conda)
    * [Debian/Ubuntu (Trusty/Xenial)](#debianubuntu-trustyxenial)
    * [Running the tests](#running-the-tests)
  * [Usage](#usage)
  * [License](#license)
  * [Feedback/Issues](#feedbackissues)

## Introduction
Munges GFF3 files exported from Chado (http://www.gmod.org/) database to make them suitable for loading into WebApollo. Currently this involved transferring annotations from polypeptide features to the feature (e.g. mRNA) from which the polypeptide derives.

## Installation
There are a number of ways to install GFF munger and details are provided below. If you encounter an issue when installing GFF munger please contact your local system administrator. If you encounter a bug please log it [here](https://github.com/sanger-pathogens/gffmunger/issues) or email us at path-help@sanger.ac.uk.

### Using pip3
```
pip3 install git+git://github.com/sanger-pathogens/gffmunger.git
```
### Conda
Set up bioconda channel:

`conda config --add channels bioconda`

Install GFF munger:

`conda install -c bioconda gffmunger`

### Debian/Ubuntu (Trusty/Xenial)
To install Python3 on Ubuntu, as root run:
```
apt-get update -qq
apt-get install -y git python3 python3-setuptools python3-biopython python3-pip
pip3 install git+git://github.com/sanger-pathogens/gffmunger.git
```
### Running the tests
The test can be run from the top level directory:  
```
./run_tests.sh
```
## Usage

```
gffmunger [--input chado_export.gff3.gz] [--fasta chado_export.fasta] [--output webapollo_compatible.gff3] [--quiet|--verbose]
```

Without `--input`, will read from standard input; without `--output`, will write new GFF3 to standard output.  If  `--fasta` is not used, then will attempt to read FASTA data from the input GFF3 file.

## License
GFF munger is free software, licensed under [GPLv3](https://github.com/sanger-pathogens/gffmunger/blob/master/LICENSE).

## Feedback/Issues
Please report any issues to the [issues page](https://github.com/sanger-pathogens/gffmunger/issues) or email path-help@sanger.ac.uk.
