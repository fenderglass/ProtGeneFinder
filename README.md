SpectroGene
===========

Version 0.1 beta

Description
-----------

SpectroGene is a tool for top-down protein identification using only bacterial genome sequence.
It can be used for annotation of new genes and correction of existing ones. The software is
currently under active development.


Availability
------------

Currently only linux 64-bit systems are supported. We are planning to add MacOS support soon.

Dependencies
------------

* Biopython [biopython.org]

Install
-------
SpectroGene is written on Python and requires no installation. However, it requires some extra Python
packages listed in "Dependencies" section.

Usage
-----
    usage: spectrogene.py genome [-h] [-p NUM_PROC] [-e E_VALUE]
                                 genome_file spectrum_file output_dir

    positional arguments:
      genome_file           path to genome file
      spectrum_file         path to spectrum file
      output_dir            output directory

    optional arguments:
      -h, --help            show this help message and exit
      -p NUM_PROC, --num_proc NUM_PROC
                            number of processes
      -e E_VALUE, --evalue E_VALUE
                            E-value threshold for PrSM [default 0.01]


Authors
-------
- Mikhail Kolmogorov (UCSD)


Contacts
--------
Please report any bugs directly to the issue tracker of this project.
Also, you can send your feedback at fenderglass@gmail.com
