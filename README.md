# GFF munger

Munges GFF3 files exported from Chado (http://www.gmod.org/) database to make them suitable for loading into WebApollo.

Currently this involved transferring annotations from polypeptide features to the feature (e.g. mRNA) from which the polypeptide derives.

## Synopsis

```
gffmunger [--input chado_export.gff3.gz] [--fasta chado_export.fasta] [--output webapollo_compatible.gff3] [--quiet|--verbose]
```

Without `--input`, will read from standard input; without `--output`, will write new GFF3 to standard output.  If  `--fasta` is not used, then will attempt to read FASTA data from the input GFF3 file.


[![Build Status](https://travis-ci.org/sanger-pathogens/gffmunger.svg?branch=master)](https://travis-ci.org/sanger-pathogens/gffmunger)

# Installation
The only dependancy is Python3. Assuming you have python 3.6+ and pip installed, just run:
```
pip3 install git+git://github.com/sanger-pathogens/gffmunger.git
```

## Debian/Ubuntu (Trusty/Xenial)
To install Python3 on Ubuntu, as root run:
```
apt-get update -qq
apt-get install -y git python3 python3-setuptools python3-biopython python3-pip
pip3 install git+git://github.com/sanger-pathogens/gffmunger.git
```
